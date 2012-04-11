#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <assert.h>
#include <vector>
#include <string>
#include <set>
#include "utils.hh"
#include "ncbidata.hh"
#include "constants.hh"

Taxonomy* parseNCBIFlatFiles( const std::string& nodes_filename, const std::string& names_filename, const std::vector< std::string >* ranks_to_mark ) {

	std::cerr << "constructing taxonomy tree from NCBI dump files...";
	Taxonomy* tax = new Taxonomy;
	std::set< const std::string* > specialranks;

	if( ranks_to_mark ) {
		for( std::vector< std::string >::const_iterator rank_it = ranks_to_mark->begin(); rank_it != ranks_to_mark->end(); ++rank_it ) {
// 			std::cerr << "inserting special rank: " << *rank_it << std::endl;
			specialranks.insert( &(tax->insertRankInternal( *rank_it )) );
		}
	}

	std::multimap< TaxonID, TaxonID > children;
	std::map< TaxonID, TaxonAnnotation* > annotation;

	{
	 	std::string line, rank;
		TaxonID taxid, parent_taxid;

	   	{
		std::list<std::string> fields;
		std::list<std::string>::iterator field_it;

		// process nodes.dmp
		std::ifstream nodesfile( nodes_filename.c_str() );
		while( std::getline( nodesfile, line ) ) {
			tokenizeMultiCharDelim( line, fields, "\t|\t", 3 );
			field_it = fields.begin();
			taxid = boost::lexical_cast< TaxonID >( *field_it++ );
			parent_taxid = boost::lexical_cast< TaxonID >( *field_it++ );
			rank = *field_it;
			children.insert( std::make_pair( parent_taxid, taxid ) );
			annotation[ taxid ] = new TaxonAnnotation( tax->insertRankInternal( rank ) );
			fields.clear();
		}
		nodesfile.close(); //close nodes.dmp
		}

		// process names.dmp
		{
		std::vector<std::string> fields;

		std::ifstream namesfile( names_filename.c_str() );
		while( std::getline( namesfile, line ) ) {
			tokenizeMultiCharDelim( line, fields, "\t|\t", 4 );
			if( fields[3] == "scientific name\t|" ) { //stupid NCBI row separator
				taxid = boost::lexical_cast< TaxonID >( fields[0] );
				annotation[ taxid ]->name = fields[1];
			}
			fields.clear();
		}
		namesfile.close();
		}
	}

	// build tree from root (taxid 1) to leaves and set nested set and rootpath values for all nodes and max_depth_ for taxonomy
	std::multimap< TaxonID, TaxonID >::reverse_iterator child_rev_it;
	std::pair< std::multimap< TaxonID, TaxonID >::reverse_iterator, std::multimap< TaxonID, TaxonID >::reverse_iterator > children_range;

	small_unsigned_int depth_counter = 0;
	large_unsigned_int lrvalue_counter = 0;
	TaxonID node_taxid = 1; //NCBI root taxid
	small_unsigned_int max_depth = 0;

	// remove root self-link from children
	std::multimap< TaxonID, TaxonID >::iterator tmp_it = children.find( node_taxid );
	if( tmp_it != children.end() ) {
		children.erase( tmp_it ); //because keys are sorted this is the pair (1,1)
	} else {
		std::cerr << " (could not find any nodes linking to the root, this will be quite a small taxonomy!)";
	}

	Taxon* tmptaxon = new Taxon( annotation[ node_taxid ] );
	tmptaxon->taxid = node_taxid;
	tmptaxon->root_pathlength = depth_counter++;
	tmptaxon->leftvalue = lrvalue_counter;

	// check wether to mark the node as special
	if( specialranks.find( &(tmptaxon->annotation->rank) ) != specialranks.end() ) {
		tmptaxon->mark_special = true;
	}

	Taxonomy::iterator child_node_it, node_it = tax->set_head( tmptaxon );
	tax->addToIndex( node_taxid, node_it.node );
	Taxonomy::sibling_iterator sibling_node_it;
	const Taxonomy::iterator root_it = node_it;

// 	std::cout << std::setfill('-'); //DEBUG

	// depth-first construction
	do {

		// set node's leftvalue
		(*node_it)->leftvalue = ++lrvalue_counter;
// 		std::cout << "*" << std::setw( depth_counter ) << "*" << (*node_it)->annotation->name << std::endl;

		// check whether to set as node to unclassified

		// get children
		children_range = children.equal_range( node_taxid );

		if( children_range.first != children_range.second ) { //means that node has at least one child
			for( child_rev_it = children_range.second; child_rev_it != children_range.first; ++child_rev_it ) { //from last to first

				// make new node and append to parent node
				node_taxid = child_rev_it->second;
				tmptaxon = new Taxon( annotation[ node_taxid ] );
				tmptaxon->taxid = node_taxid;
				tmptaxon->root_pathlength = depth_counter;

				// check whether to mark the node as unclassified
				if( (*node_it)->is_unclassified ) {
					tmptaxon->is_unclassified = true;
				} else {
					if( tmptaxon->annotation->name.find( "unclassified" ) != std::string::npos ) {
						tmptaxon->is_unclassified = true;
					}
				}

				// check wether to mark the node as special
				if( specialranks.find( &(tmptaxon->annotation->rank) ) != specialranks.end() ) {
					tmptaxon->mark_special = true; //other children shall be marked, too
				}

				child_node_it = tax->prepend_child( node_it, tmptaxon );
				tax->addToIndex( node_taxid, child_node_it.node );

// 				std::cout << "new child: " << tmptaxon->annotation->name << "; taxid: " << node_taxid << std::endl;
			}
			//children.erase( children_range.first, children_range.second ); //balanced tree is better to clear at destruction time

			assert( node_it != child_node_it ); //DEBUG

			node_it = child_node_it; //go to first child
			++depth_counter;

		} else { //this is a leaf node, go to next (or any ancestor's) sibling
			max_depth = std::max( max_depth, (*node_it)->root_pathlength );
			do {
				(*node_it)->rightvalue = ++lrvalue_counter;
				if( node_it != root_it ) {
					sibling_node_it = node_it.node->next_sibling;
					if( sibling_node_it != tax->end( sibling_node_it ) ) { //this should really work on the root
// 						std::cout << "walking to next sibling: " << std::endl;
						node_it = sibling_node_it;
						node_taxid = (*node_it)->taxid;
						break;
					} else {
						node_it = node_it.node->parent;
// 						std::cout << "walking up to parent: " << (*node_it)->annotation->name << std::endl;
						--depth_counter;
					}
				} else {
					// construction end reached
					tax->setMaxDepth( max_depth );
// 					tax->setMaxDepth();
					std::cerr << " done" << std::endl;
					return tax; //bad style but efficient
				}
			} while( true );
		}
	} while( true ); //single exit condition is return

	std::cerr << " done" << std::endl;
	tax->setMaxDepth( max_depth );
// 	tax->setMaxDepth();
	return tax;
}



Taxonomy* loadTaxonomyFromEnvironment( const std::vector< std::string >* ranks_to_mark ) {
	char* env = getenv( ENVVAR_TAXONOMY_ROOT.c_str() );
	if( env == NULL ) {
		std::cerr << "Specify the folder containing the NCBI taxonomy dump files as " << ENVVAR_TAXONOMY_ROOT << " environment variable" << std::endl;
		return NULL;
	}

	const std::string ncbi_root_folder = env;
	const std::string nodes_filename = ncbi_root_folder + "/nodes.dmp";
	const std::string names_filename = ncbi_root_folder + "/names.dmp";

	if ( boost::filesystem::exists( nodes_filename ) ) {
		if ( boost::filesystem::exists( names_filename ) ) {
			return parseNCBIFlatFiles( nodes_filename, names_filename, ranks_to_mark );
		} else {
			std::cerr << " " << names_filename << "not found" << std::endl;
		}
	} else {
		std::cerr << " " << nodes_filename << "not found" << std::endl;
	}
	return NULL;
}



const std::string extractFastaCommentField( const std::string& comment, const std::string& key ) {
	const std::size_t key_length = key.size();
	std::list< std::string > fields;
	tokenizeSingleCharDelim( comment, fields, "|" ); //NCBI scheme

	std::list< std::string >::iterator field_it = fields.begin();
	while ( field_it != fields.end() ) {
		if ( *field_it == key || ( field_it->size() > key_length && field_it->substr( field_it->size() - key_length, key_length ) == key ) ) {
			++field_it;
			break;
		}
		++field_it;
	}
	if ( field_it != fields.end() ) return *field_it;
	return fields.back(); //default behavior if not found
}
