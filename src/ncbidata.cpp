#include <fstream>
// #include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <assert.h>
#include <vector>
#include <string>
#include <set>
#include "utils.hh"
#include "ncbidata.hh"

Taxonomy* parseNCBIFlatFiles( const std::string& nodes_filename, const std::string& names_filename, const std::vector< std::string >* ranks_to_mark ) {

	Taxonomy* tax = new Taxonomy;
	std::set< const TTPString* > specialranks;

	if( ranks_to_mark ) {
		for( std::vector< std::string >::const_iterator rank_it = ranks_to_mark->begin(); rank_it != ranks_to_mark->end(); ++rank_it ) {
// 			std::cerr << "inserting special rank: " << *rank_it << std::endl;
			specialranks.insert( &(tax->getRankInternal( *rank_it, true )) );
		}
	}

	std::multimap< unsigned int, unsigned int > children;
	std::map< unsigned int, TaxonAnnotation* > annotation;

	{
	 	std::string line, rank;
		unsigned int taxid, parent_taxid;

	   	{
		std::list<std::string> fields;
		std::list<std::string>::iterator field_it;

		// process nodes.dmp
		std::ifstream nodesfile( nodes_filename.c_str() );
		while( std::getline( nodesfile, line ) ) {
			tokenizeMultiCharDelim( line, fields, "\t|\t", 3 );
			field_it = fields.begin();
			taxid = boost::lexical_cast< unsigned int >( *field_it++ );
			parent_taxid = boost::lexical_cast< unsigned int >( *field_it++ );
			rank = *field_it;
			children.insert( std::make_pair( parent_taxid, taxid ) );
			annotation[ taxid ] = new TaxonAnnotation( tax->getRankInternal( rank, true ) );
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
				taxid = boost::lexical_cast< unsigned int >( fields[0] );
				annotation[ taxid ]->name = fields[1];
			}
			fields.clear();
		}
		namesfile.close();
		}
	}

	// build tree from root (taxid 1) to leaves and set nested set and rootpath values for all nodes
	std::multimap< unsigned int, unsigned int >::reverse_iterator child_rev_it;
	std::pair< std::multimap< unsigned int, unsigned int >::reverse_iterator, std::multimap< unsigned int, unsigned int >::reverse_iterator > children_range;

	int depth_counter = 0, lrvalue_counter = 0;
	unsigned int node_taxid = 1; //NCBI root taxid

	// remove root self-link from children
	std::multimap< unsigned int, unsigned int >::iterator tmp_it = children.find( node_taxid );
	if( tmp_it != children.end() ) {
		children.erase( tmp_it ); //because keys are sorted this is the pair (1,1)
	} else {
		std::cerr << "Could not find any nodes linking to the root, this will be quite a small taxonomy!";
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

	std::cout << std::setfill('-'); //DEBUG

	// depth-first construction
	do {

		// set node's leftvalue
		(*node_it)->leftvalue = ++lrvalue_counter;
// 		std::cout << "*" << std::setw( depth_counter ) << "*" << (*node_it)->annotation->name << std::endl;

		// check whether to set as node to anclassified

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

					return tax; //bad style but efficient
				}
			} while( true );
		}
	} while( true ); //single exit condition is return

	return tax;
}



Taxonomy* loadTaxonomyFromEnvironment( const std::vector< std::string >* ranks_to_mark ) {
	char* env = getenv("TTP_NCBI_ROOT");
	if( env == NULL ) {
		std::cerr << "Specify the folder containing the NCBI taxonomy dump files as TTP_NCBI_ROOT environment variable" << std::endl;
		return NULL;
	}

	const std::string ncbi_root_folder = env;

	//TODO: check if files exist

	return parseNCBIFlatFiles( ncbi_root_folder + "/nodes.dmp", ncbi_root_folder + "/names.dmp", ranks_to_mark );
}



const std::string extractFastaCommentField( const std::string& comment, const std::string& key ) {
	std::list< std::string > fields;
	tokenizeSingleCharDelim( comment, fields, "|" ); //NCBI scheme
	bool return_field = false;
	for( std::list< std::string >::iterator field_it = fields.begin(); field_it != fields.end(); ++field_it ) {
		if( return_field ) {
			return *field_it;
		}

		if( *field_it == key ) {
			return_field = true;
		}
	}
	return fields.back(); //default behavior if not found
}

