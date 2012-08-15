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
// #include "src/ostreamtwrapperts.hh"



using namespace std;

typedef list< AlignmentRecordTaxonomy* > RecordSetType;

void doPredictionsSerial( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments = false ) {
	AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( seqid2taxid, tax );
	AlignmentFileParser< AlignmentRecordTaxonomy > parser( cin, fac );
	RecordSetGenerator< AlignmentRecordTaxonomy > recgen( parser );
	
	std::queue< RecordSetType > workload; //ready for parallelization

	//run sequential code instead of parallel
	RecordSetType tmprset;
	PredictionRecord prec( tax );
	
	while( recgen.notEmpty() ) {
		
		// parse
		recgen.getNext( tmprset );
		if ( split_alignments ) separateAlignmentsByRange( tmprset, workload ); //one query can have several aligning regions
		else workload.push( tmprset );
		tmprset.clear(); //ownership was transferred, clear for use in next cycle
		
		// run
		while ( ! workload.empty() ) {
			RecordSetType& rset = workload.front();
			predictor->predict( rset, prec );
			deleteRecords( rset ); //tidy up
			cout << prec << std::flush; //TODO: remove flush
			workload.pop();
		}
	}
}



class BoostProducer {
	public:
		BoostProducer( BoundedBuffer< RecordSetType >& buffer, AlignmentRecordFactory< AlignmentRecordTaxonomy >& fac, bool split_alignments ) : buffer_( buffer ), fac_( fac ), split_alignments_( split_alignments ) {}
		void operator()() { produce(); }
		
	private:
		BoundedBuffer< RecordSetType >& buffer_;
		AlignmentRecordFactory< AlignmentRecordTaxonomy >& fac_;
		bool split_alignments_;
		
		void produce() {
// 			std::cerr << "PRODUCER STARTED" << std::endl;
			AlignmentFileParser< AlignmentRecordTaxonomy > parser( cin, fac_ );
			RecordSetGenerator< AlignmentRecordTaxonomy > recgen( parser );
			RecordSetType tmprset;
			
			while( recgen.notEmpty() ) {
				recgen.getNext( tmprset );
				if ( split_alignments_ ) separateAlignmentsByRange( tmprset, buffer_ );
				else buffer_.push( tmprset );
				tmprset.clear(); //ownership was transferred, clear for use in next cycle
			}
			
// 			std::cerr << "PRODUCER FINISH!" << std::endl;			
		}
		
};



class BoostConsumer {
	public:
		BoostConsumer( BoundedBuffer< RecordSetType >& buffer, TaxonPredictionModel< RecordSetType >* predictor, const Taxonomy* tax, std::ostream& strm ) : buffer_( buffer ), predictor_( *predictor ), tax_( tax ), strm_( strm ) {}
		void operator()() { consume(); }

	private:
		BoundedBuffer< RecordSetType >& buffer_;
		TaxonPredictionModel< RecordSetType >& predictor_;
		const Taxonomy* tax_;
		std::ostream& strm_;
		
		void consume() {
			std::ostringstream obuf;
			PredictionRecord prec( tax_ );
			
			while ( true ) {
				// must be thread save!
				RecordSetType rset;
				try {
					rset = buffer_.pop();
				} catch ( boost::thread_interrupted ) {
					break;
				}
				
				// run prediction
				predictor_.predict( rset, prec );
				deleteRecords( rset ); //tidy up
				
				// thread-safe output to ostream
				obuf << prec;
				strm_ << obuf.str(); // should be atomic under POSIX
				obuf.str( "" ); //clear sstream
			}
		}
};



void doPredictionsParallel( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments, uint number_threads  ) {
	AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( seqid2taxid, tax );
	
	//adjust thread number
	uint procs = boost::thread::hardware_concurrency();
	if ( ! number_threads ) number_threads = procs;
	else if ( procs ) number_threads = std::min( number_threads, procs ); 
	
	BoundedBuffer< RecordSetType > buffer( 3*(number_threads - 1) ); //hold three data chunks per consumer
	
	BoostProducer producer( buffer, fac, split_alignments );
	BoostConsumer consumer( buffer, predictor, tax, std::cout );
	
	// start the consumers that wait for data in buffer
	boost::thread_group t_consumers;
	for( uint i = 0; i < number_threads - 1; ++i ) t_consumers.create_thread( boost::ref( consumer ) );
	
	producer(); //main thread is the producer that fills the buffer
	
	buffer.waitUntilEmpty();
	t_consumers.interrupt_all(); //tell waiting consumers to quit, there will be no more data coming
	t_consumers.join_all();
	
	assert( buffer.empty() ); //TODO: remove
}



void doPredictions( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments, uint number_threads ) {
	if ( number_threads > 1 ) return doPredictionsParallel( predictor, seqid2taxid, tax, split_alignments, number_threads );
	doPredictionsSerial( predictor, seqid2taxid, tax, split_alignments );
}



int main( int argc, char** argv ) {

	vector< string > ranks;
	string accessconverter_filename, algorithm, query_filename, db_filename, whitelist_filename;
	bool delete_unmarked, split_alignments;
	uint nbest, minsupport, number_threads;
	float toppercent, minscore;
	double maxevalue;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "seqid-conv-file,g", po::value< string >( &accessconverter_filename ), "filename of seqid->taxid mappings" )
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks" )
	( "algorithm,a", po::value< string >( &algorithm )->default_value( "lca" ), "set the algorithm that is used to predict taxonomic ids from alignments" )
	( "toppercent,t", po::value< float >( &toppercent )->default_value( 0.05 ), "RPA re-evaluation band or top percent parameter for LCA methods" )
	( "minscore,m", po::value< float >( &minscore )->default_value( 0.0 ), "min score parameter for MEGAN classification" )
	( "nbest,n", po::value< uint >( &nbest )->default_value( 1 ), "n-best LCA classification parameter" )
	( "max-evalue,e", po::value< double >( &maxevalue )->default_value( 1000.0 ), "set maximum evalue for filtering" )
	( "min-support,c", po::value< uint >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering) for MEGAN algorithm" )
	( "db-whitelist,w", po::value< string >( &whitelist_filename ), "specifiy list of sequence identifiers in reference to be used to reduce memory footprint (RPA algorithm)" )
	( "query-sequence-file,q", po::value< string >( &query_filename ), "fasta file to query sequences (respect order of alignments file!)" )
	( "ref-sequence-file,f", po::value< string >( &db_filename ), "fasta file to DB sequences" )
	( "ignore-unclassified,i", "alignments for partly unclassified taxa will be ignored" )
	( "split-alignments,s", po::value< bool >( &split_alignments )->default_value( true ), "decompose alignments into disjunct segments and treat them separately (for algorithms where applicable)" )
	( "processors,p", po::value< uint >( &number_threads )->default_value( 1 ), "sets number of threads, number > 2 will heavily profit from multi-core architectures, set to 0 for max. performance" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

	if( ! vm.count( "ranks" ) ) { //set to fallback if not given
		ranks = default_ranks;
	}

	if( ! vm.count( "seqid-conv-file" ) ) { //set to fallback if not given
		cout << "Specify seqid converter file" << endl;
		cout << desc << endl;
		return EXIT_FAILURE;
	}
	
	bool ignore_unclassified = vm.count( "ignore-unclassified" );
	
	// create taxonomy
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;
	
	if( delete_unmarked ) { //do everything only with the major NCBI ranks given by "ranks"
	tax->deleteUnmarkedNodes();
	}

	boost::scoped_ptr< StrIDConverter > seqid2taxid( loadStrIDConverterFromFile( accessconverter_filename, 1000 ) );
	
	// choose appropriate prediction model from command line parameters
	{
		if( algorithm == "lca" ) {
			doPredictions( &LCASimplePredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
		} else {
			if( algorithm == "megan-lca" ) {
				// check for parameters
				doPredictions( &MeganLCAPredictionModel< RecordSetType >( tax.get(), ignore_unclassified, toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments, number_threads );
			} else {
				if( algorithm == "ic-megan-lca" ) {
					doPredictions( &ICMeganLCAPredictionModel< RecordSetType >( tax.get(), toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments, number_threads );
				} else {
					if( algorithm == "n-best-lca" ) {
						doPredictions( &NBestLCAPredictionModel< RecordSetType >( tax.get(), nbest ), *seqid2taxid, tax.get(), split_alignments, number_threads );
					} else { //the rest are hidden algorithms that use a known query taxid from the identifier
						if( algorithm == "best-in-tree" ) {
							tax->setRankDistances( ranks );
							doPredictions( &ClosestNodePredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
						} else {
							if( algorithm == "correction" ) {
								tax->setRankDistances( ranks );
								doPredictions( &CorrectionPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
							} else {
								if( algorithm == "ic-correction" ) {
									tax->setRankDistances( ranks );
									doPredictions( &ICCorrectionPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
								} else {
									if( algorithm == "query-best-lca" ) {
										doPredictions( &QueryBestLCAPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
									} else {
										if( algorithm == "rpa" ) {
											typedef RandomSeqStorRO< seqan::String< seqan::Dna5 > > DBStorType;
											typedef RandomSeqStorRO< seqan::String< seqan::Dna5 > > QStorType;
											
											QStorType query_storage( query_filename );

											StopWatchCPUTime measure_db_loading( "loading reference db" );
											boost::scoped_ptr< DBStorType > db_storage;
											
											measure_db_loading.start();
											if ( vm.count( "db-whitelist" ) ) {
												set< string > whitelist;
												populateIdentSet( whitelist, whitelist_filename );
												db_storage.reset( new DBStorType( db_filename, whitelist ) );
											} else {
												db_storage.reset( new DBStorType( db_filename ) );
											}
											measure_db_loading.stop();
											
											std::ofstream prediction_debug_output( "prediction.log" );
											doPredictions( &DoubleAnchorRPAPredictionModel< RecordSetType, QStorType, DBStorType >( tax.get(), query_storage, *db_storage, toppercent, prediction_debug_output ), *seqid2taxid, tax.get(), split_alignments, number_threads ); //TODO: reuse toppercent param?
										} else {
											if( algorithm == "dummy" ) {
												doPredictions( &DummyPredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get(), split_alignments, number_threads );
											} else {
												cout << "classification algorithm can either be: rpa, lca, megan-lca, ic-megan-lca, n-best-lca" << endl;
												return EXIT_FAILURE;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return EXIT_SUCCESS;
}
