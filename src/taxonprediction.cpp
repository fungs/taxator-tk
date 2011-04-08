#include <list>
#include <string>
#include <stack>
#include <boost/lexical_cast.hpp>
#include "taxonprediction.hh"
#include "utils.hh"
#include "constants.hh"



void MyTaxonPredictionModel::loadModel( const std::string& filename ) {
	loadMapFromFile( filename, bitscore_cutoff, FSEP );
}



TaxonNode* MyTaxonPredictionModel::predict( AlignmentRecord* record ) {
		unsigned int taxid = record->reference_taxid;
		TaxonNode* node = taxinter.getNode( taxid );

		// generate path from root to leaf
		std::stack< TaxonNode* > working_list;
		const TaxonNode* root_node = taxinter.getRoot();
		while( node != root_node ) {
			if( node->data->mark_special ) {
				working_list.push( node );
			}
			node = node->parent;
		}
		// node = root_node here

		// traverse decision tree
		while( ! working_list.empty() ) {
			taxid = working_list.top()->data->taxid;
// 			std::cerr << "treshold at node " << working_list.top()->data->taxid << " is " << bitscore_cutoff[ taxid ] << std::endl;
			if( record->bitscore <= bitscore_cutoff[ taxid ]  ) {
// 				std::cout << "bitscore " << record->bitscore << " < threshold " << bitscore_cutoff[ taxid ] << ", stopping..." << std::endl;
				break;
			}
			node = working_list.top();
			working_list.pop();
		}
		return node;
}


void MyTaxonPredictionModel2::loadModel( const std::string& filename ) {
	std::string line;
	std::list< std::string > fields;
	std::list< std::string >::iterator field_it;
	std::string rank;
	float bthreshold;

	//TODO: check file existance etc.
	std::ifstream file_handle( filename.c_str() );

	while( std::getline( file_handle, line ) ) {
		if( ! ignoreLine( line ) ) {
			tokenizeSingleCharDelim( line, fields, FSEP, 2 );
			field_it = fields.begin();
			rank = boost::lexical_cast< std::string >( *field_it++ );
			bthreshold = boost::lexical_cast< float >( *field_it );
			ranks_cutoff[ &taxonomy->getRankInternal( rank ) ] = bthreshold;
			fields.clear();
		}
	}

	file_handle.close();
}



TaxonNode* MyTaxonPredictionModel2::predict( AlignmentRecord* record ) {
	TaxonNode* node = taxinter.getNode( record->reference_taxid );
	const TaxonNode* root_node = taxinter.getRoot();
	float bitscore = record->bitscore;

	// traverse decision tree
	while( node != root_node ) {
		if( node->data->mark_special ) {
			std::map< const std::string*, float >::iterator find_it = ranks_cutoff.find( &node->data->annotation->rank );
			if( find_it != ranks_cutoff.end() && bitscore >= find_it->second ) {
				break;
			}
		}
		node = node->parent;
	}
	return node;
}



void MyTaxonPredictionModel3::loadModel( const std::string& filename ) {
	std::string line;
	std::list< std::string > fields;
	std::list< std::string >::iterator field_it;
	float bthreshold;

	//TODO: check file existance etc.
	std::ifstream file_handle( filename.c_str() );

	while( std::getline( file_handle, line ) ) {
		if( ! ignoreLine( line ) ) {
			tokenizeSingleCharDelim( line, fields, FSEP, 2 );
			field_it = fields.begin();
			std::string rank = *field_it++;
			bthreshold = boost::lexical_cast< float >( *field_it );
			std::cerr << "storing ranks threshold for " << rank << ": " << bthreshold << std::endl;
			ranks_cutoff[ &taxonomy->getRankInternal( rank ) ] = bthreshold;
			fields.clear();
		}
	}

	file_handle.close();
}



TaxonNode* MyTaxonPredictionModel3::predict( AlignmentRecord* record ) {
	TaxonNode* node = taxinter.getNode( record->reference_taxid );
	const TaxonNode* root_node = taxinter.getRoot();
	float bitscore = record->bitscore;

	// traverse decision tree
	while( node != root_node ) {
		if( node->data->mark_special ) {
			std::map< const std::string*, float >::iterator find_it = ranks_cutoff.find( &node->data->annotation->rank );
			if( find_it != ranks_cutoff.end() && bitscore/max_bitscore >= find_it->second ) {
				break;
			}
		}
		node = node->parent;
	}
// 	std::cerr << "predicting single node with rel. bitscore " << bitscore/max_bitscore << " at level " << node->data->annotation->rank << " (" << record->reference_taxid << ")" << std::endl;
	return node;
}



void MyTaxonPredictionModel3::setMaxBitscore( const float bitscore ) {
	max_bitscore = bitscore;
}
