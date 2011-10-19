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

#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include <iostream>
#include <queue>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/alignmentrecord.hh"
#include "src/taxonpredictionmodelsequence.hh"
#include "src/taxonpredictionmodel.hh"
#include "src/constants.hh"
#include "src/sequencestorage.hh"
#include "src/predictionrecord.hh"



using namespace std;

typedef list< AlignmentRecordTaxonomy* > RecordSetType;
void doPredictions( TaxonPredictionModel< RecordSetType >* predictor, StrIDConverter& seqid2taxid, const Taxonomy* tax, bool split_alignments = false ) {
	AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( seqid2taxid, tax );
	AlignmentFileParser< AlignmentRecordTaxonomy > parser( cin, fac );
	RecordSetGenerator< AlignmentRecordTaxonomy > recgen( parser );
	PredictionRecord prec( tax );
	std::queue< RecordSetType > workload; //ready for parallelization
	RecordSetType tmprset;
	
	while( recgen.notEmpty() ) {
		recgen.getNext( tmprset );
		
		if ( split_alignments ) separateAlignmentsByRange( tmprset, workload ); //one query can have several aligning regions
		else workload.push( tmprset );
		tmprset.clear(); //ownership was transferred, clear for use in next cycle
		
		while ( ! workload.empty() ) { //TODO: parallelization
			RecordSetType& rset = workload.front();
			prec.query_identifier = rset.front()->getQueryIdentifier();
			predictor->predict( rset, prec );
			cout << prec << std::flush;
			deleteRecords( rset );
			workload.pop();
		}
	}
}



int main( int argc, char** argv ) {

	vector< string > ranks;
	string accessconverter_filename, algorithm, query_filename, db_filename;
	bool delete_unmarked, split_alignments;
	unsigned int nbest, minsupport;
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
	( "toppercent,t", po::value< float >( &toppercent )->default_value( 0.3 ), "top percent parameter for MEGAN classification or parameter for extended-lca method" )
	( "minscore,m", po::value< float >( &minscore )->default_value( 0.0 ), "min score parameter for MEGAN classification" )
	( "nbest,n", po::value< unsigned int >( &nbest )->default_value( 1 ), "parameter for n-best LCA classification" )
	( "max-evalue,e", po::value< double >( &maxevalue )->default_value( 1000.0 ), "set maximum evalue for filtering" )
	( "min-support,c", po::value< unsigned int >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering)" )
	( "query-sequence-file,q", po::value< string >( &query_filename ), "fasta file to query sequences (respect order of alignments file!)" )
	( "ref-sequence-file,b", po::value< string >( &db_filename ), "fasta file to DB sequences" )
	( "ignore-unclassified,i", "alignments for partly unclassified taxa are not considered" )
	( "split-alignments,s", po::value< bool >( &split_alignments )->default_value( true ), "decompose alignments into disjunct pieces and treat each region separately (for algorithms where applicable)" );

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
	boost::scoped_ptr< int > test( NULL );
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;
	
	if( delete_unmarked ) { //do everything only with the major NCBI ranks given by "ranks"
	tax->deleteUnmarkedNodes();
	}

	boost::scoped_ptr< StrIDConverter > seqid2taxid( loadStrIDConverterFromFile( accessconverter_filename, 1000 ) );
	
	// choose appropriate prediction model from command line parameters
	{
		if( algorithm == "lca" ) {
			doPredictions( &LCASimplePredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get(), split_alignments );
		} else {
			if( algorithm == "megan-lca" ) {
				// check for parameters
				doPredictions( &MeganLCAPredictionModel< RecordSetType >( tax.get(), ignore_unclassified, toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments );
			} else {
				if( algorithm == "ic-megan-lca" ) {
					doPredictions( &ICMeganLCAPredictionModel< RecordSetType >( tax.get(), toppercent, minscore, minsupport, maxevalue ), *seqid2taxid, tax.get(), split_alignments );
				} else {
					if( algorithm == "n-best-lca" ) {
						doPredictions( &NBestLCAPredictionModel< RecordSetType >( tax.get(), nbest ), *seqid2taxid, tax.get() );
					} else { //the rest are hidden algorithms that use a known query taxid from the identifier
						if( algorithm == "best-in-tree" ) {
							tax->setRankDistances( ranks );
							doPredictions( &ClosestNodePredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get() );
						} else {
							if( algorithm == "correction" ) {
								tax->setRankDistances( ranks );
								doPredictions( &CorrectionPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get() );
							} else {
								if( algorithm == "ic-correction" ) {
									tax->setRankDistances( ranks );
									doPredictions( &ICCorrectionPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get() );
								} else {
									if( algorithm == "query-best-lca" ) {
										doPredictions( &QueryBestLCAPredictionModel< RecordSetType >( tax.get(), seqid2taxid.get() ), *seqid2taxid, tax.get() );
									} else {
										if( algorithm == "test" ) {
											typedef RandomSeqStorRO< seqan::String< seqan::Dna5 > > DBStorType;
											typedef SequentialSeqStorRO< seqan::Dna5String > QStorType;
											QStorType query_storage( query_filename );
											DBStorType db_storage( db_filename );
											doPredictions( &TestPM< RecordSetType, QStorType, DBStorType >( tax.get(), query_storage, db_storage ), *seqid2taxid, tax.get() );
										} else {
											if( algorithm == "dummy" ) {
												doPredictions( &DummyPredictionModel< RecordSetType >( tax.get() ), *seqid2taxid, tax.get() );
											} else {
												cout << "classification algorithm can either be: lca, megan-lca, ic-megan-lca, n-best-lca" << endl;
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
