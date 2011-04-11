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

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/alignmentrecord.hh"
#include "src/taxonprediction.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {

	vector< string > ranks;
	string accessconverter_filename, algorithm;
	//bool delete_unmarked;
	unsigned int nbest, minsupport;
	float toppercent, minscore, weight_lower, balance_beta, t1, t2;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "seqid-conv-file,g", po::value< string >( &accessconverter_filename ), "filename of seqid->taxid mappings" )
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
// 	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks" )
	( "prediction-algorithm,a", po::value< string >( &algorithm )->default_value( "lca" ), "set the algorithm that is used to predict taxonomic ids from alignments" )
	( "toppercent,t", po::value< float >( &toppercent )->default_value( 0.3 ), "top percent parameter for MEGAN classification or parameter for extended-lca method" )
	( "minscore,m", po::value< float >( &minscore )->default_value( 0.0 ), "min score parameter for MEGAN classification" )
	( "nbest,n", po::value< unsigned int >( &nbest )->default_value( 1 ), "parameter for n-best LCA classification" )
	( "min-support,c", po::value< unsigned int >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering)" )
	( "beta-balance,b", po::value< float >( &balance_beta )->default_value( 0.5 ), "heuristic parameter for extended-lca feature balancing method" )
	( "weight-lower-bound,w", po::value< float >( &weight_lower )->default_value( 0.5 ), "heuristic parameter for extended-lca feature balancing method" )
	( "t-generic-1,1", po::value< float >( &t1 )->default_value( 0.0 ), "generic heuristic parameter for extended-lca feature balancing method" )
	( "t-geneneric-2,2", po::value< float >( &t2 )->default_value( 0.5 ), "generic heuristic parameter for extended-lca feature balancing method" )
	( "ignore-unclassified,i", "alignments for partly unclassified taxa are not considered" );

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
	Taxonomy* tax = loadTaxonomyFromEnvironment( &ranks );
	if( ! tax ) {
		return EXIT_FAILURE;
	}
// 	if( delete_unmarked ) { //do everything only with the major NCBI ranks given by "ranks"
	tax->deleteUnmarkedNodes();
// 	}

	TaxonomyInterface taxinter( tax );
	StrIDConverter* seqid2taxid = loadStrIDConverterFromFile( accessconverter_filename, 1000 );

	typedef list< AlignmentRecord* > RecordSetType;
	TaxonPredictionModel< RecordSetType >* predictor = NULL;

	// choose appropriate prediction model from command line parameters
	{
		if( algorithm == "lca" ) {
			predictor = new LCASimplePredictionModel< RecordSetType >( tax );
		} else {
			if( algorithm == "megan-lca" ) {
				// check for parameters
				predictor = new MeganLCAPredictionModel< RecordSetType >( tax, ignore_unclassified, toppercent, minscore, minsupport );
			} else {
				if( algorithm == "best-in-tree" ) {
					tax->setRankDistances( ranks );
					predictor = new ClosestNodePredictionModel< RecordSetType >( tax, seqid2taxid );
				} else {
					if( algorithm == "correction" ) {
						tax->setRankDistances( ranks );
						predictor = new CorrectionPredictionModel< RecordSetType >( tax, seqid2taxid );
					} else {
						if( algorithm == "ic-correction" ) {
							tax->setRankDistances( ranks );
							predictor = new ICCorrectionPredictionModel< RecordSetType >( tax, seqid2taxid );
						} else {
							if( algorithm == "query-best-lca" ) {
								predictor = new QueryBestLCAPredictionModel< RecordSetType >( tax, seqid2taxid );
							} else {
								if( algorithm == "n-best-lca" ) {
									predictor = new NBestLCAPredictionModel< RecordSetType >( tax, nbest );
								} else {
									if( algorithm == "extended-lca" ) {
										tax->setRankDistances( ranks );
										predictor = new ExtendedLCAPredictionModel< RecordSetType >( tax, t1, t2 );
									} else {
										if( algorithm == "top-percent-outlier-lca" ) {
											tax->setRankDistances( ranks );
											predictor = new TopPercentOutlierLCAPredictionModel< RecordSetType >( tax, toppercent, t1 );
										} else {
											if( algorithm == "ic-megan-lca" ) {
												predictor = new ICMeganLCAPredictionModel< RecordSetType >( tax, toppercent, minscore );
											} else {
												if( algorithm == "dummy" ) {
													predictor = new DummyPredictionModel< RecordSetType >( tax );
												} else {
													cout << "classification algorithm can either be: lca, megan-lca, best-in-tree, query-best-lca, correction, n-best-lca, extended-lca, dummy" << endl;
						// 								cout << desc << endl;
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
	}


	AlignmentFileParser parser( cin, seqid2taxid );
	RecordSetGenerator recgen( parser );

	RecordSetType recordset;
	while( recgen.notEmpty() ) {
		recgen.getNext( recordset );
		cout << recordset.front()->query_identifier << FSEP << predictor->predict( recordset )->data->taxid << endl;
		deleteRecords( recordset );
	}

	// tidy up and quit
	delete predictor;
	delete seqid2taxid;
	delete tax;
	return EXIT_SUCCESS;
}
