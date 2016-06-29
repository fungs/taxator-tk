/*
taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

Copyright (C) 2010 Johannes Dröge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/alignmentrecord.hh"
#include "src/taxonpredictionmodelsequence.hh"
#include "src/taxonpredictionmodel.hh"
#include "src/constants.hh"
#include "src/sequencestorage.hh"
#include "src/predictionrecord.hh"
#include "src/profiling.hh"
#include "src/boundedbuffer.hh"
#include "src/concurrentoutstream.hh"
#include "src/exception.hh"

using namespace std;

typedef list< AlignmentRecordTaxonomy* > RecordSetType;

void doPredictionsSerial( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments,bool alignments_sorted, std::ostream& logsink ) {
    AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( seqid2taxid, tax );
    FileParser< AlignmentRecordFactory< AlignmentRecordTaxonomy > > parser( cin, fac );
    RecordSetGenerator<AlignmentRecordTaxonomy, RecordSetType>* recgen = NULL; // TODO: boost smpt??
    
    if (alignments_sorted) { // stupid nesting because template parameters must be const
        if (split_alignments) recgen = new RecordSetGeneratorSorted<AlignmentRecordTaxonomy, RecordSetType, true>( parser );
        else recgen = new RecordSetGeneratorSorted<AlignmentRecordTaxonomy, RecordSetType, false>( parser );
    }
    else {
        if (split_alignments) recgen = new RecordSetGeneratorUnsorted<AlignmentRecordTaxonomy, RecordSetType, true>( parser );
        else recgen = new RecordSetGeneratorUnsorted<AlignmentRecordTaxonomy, RecordSetType, false>( parser );
    }

    RecordSetType rset;
    
    PredictionRecord prec( tax );

    std::cout << GFF3Header();
    while( recgen->notEmpty() ) {
        recgen->getNext( rset );
        predictor->predict( rset, prec, logsink );
        deleteRecords( rset );
        std::cout << prec;
    }

//     delete recgen;
}


class BoostProducer {
public:
    BoostProducer( BoundedBuffer< RecordSetType >& buffer, AlignmentRecordFactory< AlignmentRecordTaxonomy >& fac, bool split_alignments, bool alignments_sorted ) :
        buffer_( buffer ),
        fac_( fac ),
        split_alignments_( split_alignments ),
        alignments_sorted_( alignments_sorted )
    {}

    void operator()() {
        produce();
    }

private:

    BoundedBuffer< RecordSetType >& buffer_;
    AlignmentRecordFactory< AlignmentRecordTaxonomy >& fac_;
    bool split_alignments_;
    bool alignments_sorted_;

    void produce() {  //TODO: use boost smart pointers for factory
        FileParser< AlignmentRecordFactory< AlignmentRecordTaxonomy > > parser( cin, fac_ );
        RecordSetGenerator<AlignmentRecordTaxonomy, RecordSetType>* recgen = NULL;

        if (alignments_sorted_) { // stupid nesting because template parameters must be const
            if (split_alignments_) recgen = new RecordSetGeneratorSorted<AlignmentRecordTaxonomy, RecordSetType, true>( parser );
            else recgen = new RecordSetGeneratorSorted<AlignmentRecordTaxonomy, RecordSetType, false>( parser );
        }
        else {
            if (split_alignments_) recgen = new RecordSetGeneratorUnsorted<AlignmentRecordTaxonomy, RecordSetType, true>( parser );
            else recgen = new RecordSetGeneratorUnsorted<AlignmentRecordTaxonomy, RecordSetType, false>( parser );
        }
        
        RecordSetType tmprset;

        while( recgen->notEmpty() ) {
            recgen->getNext( tmprset );
            buffer_.push( tmprset );
            tmprset.clear();  // ownership transferred, clear for next cycle
        }

        delete recgen;
    }

};


class BoostConsumer {
public:
    BoostConsumer( BoundedBuffer< RecordSetType >& buffer, TaxonPredictionModel< RecordSetType >* predictor, const Taxonomy* tax, ConcurrentOutStream& log, ConcurrentOutStream& output ) :
        buffer_( buffer ),
        predictor_( *predictor ),
        tax_( tax ),
        output_( output ),
        log_( log ),
        thread_count_( 0 )
    {}

    void operator()() {
        consume();
    }

private:
    BoundedBuffer< RecordSetType >& buffer_;
    TaxonPredictionModel< RecordSetType >& predictor_;
    const Taxonomy* tax_;
    ConcurrentOutStream& output_;
    ConcurrentOutStream& log_;
    boost::mutex count_mutex_; //needed for concurrent thread count
    uint thread_count_;

    void consume() {
        PredictionRecord prec( tax_ );

        // determine count of this thread to index concurrent stream
        boost::mutex::scoped_lock count_lock( count_mutex_ );
        const uint this_thread = thread_count_++;
        count_lock.unlock();

        while ( true ) {
            RecordSetType rset;
            try {
                rset = buffer_.pop();
            } catch ( boost::thread_interrupted ) {
                break;
            }

            // run prediction
            predictor_.predict( rset, prec, log_( this_thread ) );
            log_.flush( this_thread );

            // output to stdout
            output_( this_thread ) << prec;
            output_.flush( this_thread );

            deleteRecords( rset );
        }
    }
};


void doPredictionsParallel( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments, bool alignments_sorted , std::ostream& logsink, uint number_threads  ) {
    AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( seqid2taxid, tax );

    //print GFF3Header
    std::cout << GFF3Header();

    //adjust thread number
    uint procs = boost::thread::hardware_concurrency();
    if ( ! number_threads ) number_threads = procs;  // set number of threads to available (producer thread is really lightweight)
    else if ( procs ) number_threads = std::min( number_threads, procs );

    BoundedBuffer< RecordSetType > buffer( 10*number_threads );  // hold ten data chunks per consumer TODO: make option
    ConcurrentOutStream output( std::cout, number_threads, 1000 );  // TODO: analyse number and increase buffer size
    ConcurrentOutStream log( logsink, number_threads, 20000 );

    BoostProducer producer( buffer, fac, split_alignments, alignments_sorted );
    BoostConsumer consumer( buffer, predictor, tax, log, output );

    // start the consumers that wait for data in buffer
    boost::thread_group t_consumers;
    for( uint i = 0; i < number_threads; ++i ) t_consumers.create_thread( boost::ref( consumer ) );

    producer();  // main thread is the producer that fills the buffer (not counted separately)

    buffer.waitUntilEmpty();
    t_consumers.interrupt_all();  // tell waiting consumers to quit, there will be no more data coming
    t_consumers.join_all();

    assert( buffer.empty() );  // TODO: remove
}


// TODO: use function template?
void doPredictions( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments, bool alignments_sorted, std::ostream& logsink, uint number_threads ) {
    if ( number_threads > 1 ) return doPredictionsParallel( predictor, seqid2taxid, tax, split_alignments, alignments_sorted, logsink, number_threads );
    doPredictionsSerial( predictor, seqid2taxid, tax, split_alignments, alignments_sorted, logsink );
}

template<typename StringType>
void execute(string db_filename,string db_index_filename,string query_filename, string query_index_filename,
            boost::scoped_ptr< Taxonomy >& tax, float filterout, float toppercent, boost::scoped_ptr< StrIDConverter >& seqid2taxid,//
            bool split_alignments, bool alignments_sorted, std::ofstream& logsink, uint number_threads){

   //  &RPAPredictionModel< RecordSetType, RandomSeqStoreROInterface< StringType >, RandomSeqStoreROInterface< StringType>, StringType>( tax.get(), *query_storage, *db_storage, filterout, toppercent )
   //  , *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );  // T
    
    // load query sequences
    std::cerr << "load queries\n";
    boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > query_storage;
    
    if( query_index_filename.empty() ) query_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( query_filename ) );
    else query_storage.reset( new RandomIndexedSeqstoreRO< StringType >( query_filename, query_index_filename ) );
    
    std::cerr << "end load queries\n";
    
    // reference query sequences
    boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > db_storage;
    StopWatchCPUTime measure_db_loading( "loading reference db" );
    measure_db_loading.start();
    std::cerr << "db loading\n";
    if( db_index_filename.empty() ) {
        std::cerr << "db indexfile empty";
        db_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( db_filename ) );
    }
    else {
        std::cerr << "not empty";
        db_storage.reset( new RandomIndexedSeqstoreRO< StringType >( db_filename, db_index_filename ) );
    }
    measure_db_loading.stop();
    std::cerr << "db loaded\n";
    doPredictions( &RPAPredictionModel< RecordSetType, RandomSeqStoreROInterface< StringType >, RandomSeqStoreROInterface< StringType>, StringType>( tax.get(), *query_storage, *db_storage, filterout, toppercent ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );  // TODO: reuse toppercent param?
}

int main( int argc, char** argv ) {

    vector< string > ranks;
    string accessconverter_filename, algorithm, query_filename, query_index_filename, db_filename, db_index_filename, whitelist_filename, log_filename, data_format;
    bool delete_unmarked, split_alignments, alignments_sorted;
    uint nbest, minsupport, number_threads;
    float toppercent, minscore, filterout;
    double maxevalue;

    namespace po = boost::program_options;
    po::options_description visible_options ( "Allowed options" );
    visible_options.add_options()
    ( "help,h", "show help message")
    ( "citation", "show citation info" )
    ( "advanced-options", "show advanced program options" )
    ( "algorithm,a", po::value< string >( &algorithm )->default_value( "rpa" ), "set the algorithm that is used to predict taxonomic ids from alignments" )
    ( "seqid-taxid-mapping,g", po::value< string >( &accessconverter_filename ), "filename of seqid->taxid mapping for reference" )
    ( "query-sequences,q", po::value< string >( &query_filename ), "query sequences FASTA" )
    ( "query-sequences-index,v", po::value< string >( &query_index_filename ), "query sequences FASTA index, for out-of-memory operation; is created if not existing" )
    ( "ref-sequences,f", po::value< string >( &db_filename ), "reference sequences FASTA" )
    ( "ref-sequences-index,i", po::value< string >( &db_index_filename ), "FASTA file index, for out-of-memory operation; is created if not existing" )
    ( "processors,p", po::value< uint >( &number_threads )->default_value( 1 ), "sets number of threads, number > 2 will heavily profit from multi-core architectures, set to 0 for max. performance" )
    ( "logfile,l", po::value< std::string >( &log_filename )->default_value( "/dev/null" ), "specify name of file for logging (appending lines)" )
    ( "dataformat,b", po::value< std::string >( &data_format)->default_value("nucleotide"), "specify used data (nucleotide, protein)");
    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
    ( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
    ( "split-alignments,s", po::value< bool >( &split_alignments )->default_value( true ), "decompose alignments into disjunct segments and treat them separately (for algorithms where applicable)" )
    ( "alignments-sorted,o", po::value< bool>( &alignments_sorted )->default_value( false ), "avoid sorting if alignments are sorted")
    ( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks" )
    ( "heuristic-cutoff,x", po::value<float>(&filterout)->default_value(0.5), "filter out alignments, increase means faster run-time whereas 0 means no filtering at all")
    ( "toppercent,t", po::value< float >( &toppercent )->default_value( 0.05 ), "RPA re-evaluation band or top percent parameter for LCA methods" )
    ( "max-evalue,e", po::value< double >( &maxevalue )->default_value( 1000.0 ), "set maximum evalue for filtering" )
    ( "min-support,c", po::value< uint >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering) for MEGAN algorithm" )
    ( "minscore,m", po::value< float >( &minscore )->default_value( 0.0 ), "min score parameter for MEGAN classification" )
    ( "nbest,n", po::value< uint >( &nbest )->default_value( 1 ), "n-best LCA classification parameter" )
    ( "ignore-unclassified,u", "alignments for partly unclassified taxa will be ignored" )
    ( "db-whitelist,w", po::value< string >( &whitelist_filename ), "specifiy list of sequence identifiers in reference to be used to reduce memory footprint (RPA algorithm)" );

    po::options_description all_options;
    all_options.add( visible_options ).add( hidden_options );

    po::variables_map vm;
    po::store ( po::command_line_parser ( argc, argv ).options ( all_options ).run(), vm );

    if ( vm.count ( "help" ) ) {
        cout << visible_options << endl;
        return EXIT_SUCCESS;
    }

    if ( vm.count ( "citation" ) ) {
        cout << citation_note << endl;
        return EXIT_SUCCESS;
    }

    if ( vm.count ( "advanced-options" ) ) {
        cout << hidden_options << endl;
        return EXIT_SUCCESS;
    }

    po::notify ( vm );  // check required etc.

    if( ! vm.count( "ranks" ) ) {  // set to fallback if not given
        ranks = default_ranks;
    }

    if( ! vm.count( "seqid-taxid-mapping" ) ) {
        cout << "Specify a taxonomy mapping file for the reference sequence identifiers" << endl;
        cout << visible_options << endl;
        return EXIT_FAILURE;
    }

    bool ignore_unclassified = vm.count( "ignore-unclassified" );
    
    std::cerr << "load taxonomy\n";
    boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );  // create taxonomy
    if( ! tax ) return EXIT_FAILURE;
    std::cerr << "end load taxonomy\n";
    if( delete_unmarked ) tax->deleteUnmarkedNodes();  // do everything only with the major NCBI ranks given by "ranks"
    
    std::cerr << "load idtotax\n";
    boost::scoped_ptr< StrIDConverter > seqid2taxid( loadStrIDConverterFromFile( accessconverter_filename, 1000 ) );
    std::ofstream logsink( log_filename.c_str(), std::ios_base::app );
    std::cerr << "end load idtotax\n";
    
    try {
      // choose appropriate prediction model from command line parameters
      //TODO: "address of temporary warning" is annoying but life-time is guaranteed until function returns
      if( algorithm == "dummy" ) doPredictions( &DummyPredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );
      else if( algorithm == "simple-lca" ) doPredictions( &LCASimplePredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );
      else if( algorithm == "megan-lca" ) doPredictions( &MeganLCAPredictionModel< RecordSetType >( tax.get(), ignore_unclassified, toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );
      else if( algorithm == "ic-megan-lca" ) doPredictions( &MeganLCAPredictionModel< RecordSetType >( tax.get(), ignore_unclassified, toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );
      else if( algorithm == "n-best-lca" ) doPredictions( &NBestLCAPredictionModel< RecordSetType >( tax.get(), nbest ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );
      else if( algorithm == "rpa" ) {
          
          if(data_format == "nucleotide"){
              execute<seqan::String<seqan::Dna5>>(db_filename, db_index_filename, query_filename, query_index_filename,
                                                  tax , filterout, toppercent, seqid2taxid, 
                                                  split_alignments, alignments_sorted, logsink, number_threads);
//                typedef seqan::String<seqan::Dna5> StringType;
//                ///temp func
//          
//                // load query sequences
//                std::cerr << "load queries\n";
//                boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > query_storage;
//                if( query_index_filename.empty() ) query_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( query_filename ) );
//                else query_storage.reset( new RandomIndexedSeqstoreRO< StringType >( query_filename, query_index_filename ) );
//                std::cerr << "end load queries\n";
//                // reference query sequences
//                boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > db_storage;
//                StopWatchCPUTime measure_db_loading( "loading reference db" );
//                measure_db_loading.start();
//                std::cerr << "db loading\n";
//                if( db_index_filename.empty() ) {
//                    std::cerr << "db indexfile empty";
//                    db_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( db_filename ) );
//                }
//                else {
//                    std::cerr << "not empty";
//                    db_storage.reset( new RandomIndexedSeqstoreRO< StringType >( db_filename, db_index_filename ) );
//                }
//                measure_db_loading.stop();
//                std::cerr << "db loaded\n";
//                doPredictions( &RPAPredictionModel< RecordSetType, RandomSeqStoreROInterface< StringType >, RandomSeqStoreROInterface< StringType>, StringType>( tax.get(), *query_storage, *db_storage, filterout, toppercent ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );  // TODO: reuse toppercent param?
//                    
          }
          else if(data_format == "protein"){
               execute<seqan::String<seqan::AminoAcid>>(db_filename, db_index_filename, query_filename, query_index_filename,
                                                  tax , filterout, toppercent, seqid2taxid, 
                                                  split_alignments, alignments_sorted, logsink, number_threads);
          }
//                typedef seqan::String<seqan::AminoAcid> StringType;
//                ///temp func
//          
//                // load query sequences
//                std::cerr << "load queries\n";
//                boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > query_storage;
//                if( query_index_filename.empty() ) query_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( query_filename ) );
//                else query_storage.reset( new RandomIndexedSeqstoreRO< StringType >( query_filename, query_index_filename ) );
//                std::cerr << "end load queries\n";
//                // reference query sequences
//                boost::scoped_ptr< RandomSeqStoreROInterface< StringType > > db_storage;
//                StopWatchCPUTime measure_db_loading( "loading reference db" );
//                measure_db_loading.start();
//                std::cerr << "db loading\n";
//                if( db_index_filename.empty() ) {
//                    std::cerr << "db indexfile empty";
//                    db_storage.reset( new RandomInmemorySeqStoreRO< StringType, StringType >( db_filename ) );
//                }
//                else {
//                    std::cerr << "not empty";
//                    db_storage.reset( new RandomIndexedSeqstoreRO< StringType >( db_filename, db_index_filename ) );
//                }
//                measure_db_loading.stop();
//                std::cerr << "db loaded\n";
//                doPredictions( &RPAPredictionModel< RecordSetType, RandomSeqStoreROInterface< StringType >, RandomSeqStoreROInterface< StringType>, StringType>( tax.get(), *query_storage, *db_storage, filterout, toppercent ), *seqid2taxid, tax.get(), split_alignments, alignments_sorted, logsink, number_threads );  // TODO: reuse toppercent param?
//          }
          else{
          cout << "data format can either be nucleotide or protein" << endl;
          return EXIT_FAILURE;
          }
          } else {
          cout << "classification algorithm can either be: rpa (default), simple-lca, megan-lca, ic-megan-lca, n-best-lca" << endl;
          return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    } catch(Exception &e) {
       cerr << "An unrecoverable error occurred: " << e.what() << endl;
       cerr << endl << "Here is some debugging information to locate the problem:" << endl << boost::diagnostic_information(e) << endl;

       return EXIT_FAILURE;
    }
}
