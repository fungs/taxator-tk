#include "taxontree.hh"
#include <algorithm>
#include <vector>



TaxonTree::~TaxonTree() {
	//delete all taxons and annotations
	for( iterator node_it = this->begin(); node_it != this->end(); ++node_it ) {
		delete *node_it;;
	}
}



void TaxonTree::addToIndex( unsigned int taxid , TaxonTree::Node* node ) {
	taxid2node[ taxid ] = node;
}



void TaxonTree::recreateNodeIndex() {
	taxid2node.clear();
	for( iterator node_it = this->begin(); node_it != this->end(); ++node_it ) {
		taxid2node.insert( std::make_pair( (*node_it)->taxid, node_it.node ) );
	}
}



// constant in time as apposed to size(), I think
int TaxonTree::indexSize() { //returns only real nodes (no dummies)
	return taxid2node.size();
}



const std::string& TaxonTree::getRankInternal ( const std::string& rankname, const bool insert ) {
	if( insert ) {
		return *ranks.insert( rankname ).first;
	}

	std::set< std::string >::iterator rank_it = ranks.find( rankname );
	if( rank_it == ranks.end() ) {
		return *ranks.insert( "" ).first;
	}

	return *rank_it;
}



void TaxonTree::deleteUnmarkedNodes() {
	iterator node_it = ++( this->begin() ); //root node
	while( node_it != this->end() ) {
		Taxon* deltaxon = *node_it;
		if( ! deltaxon->mark_special ) {
			taxid2node.erase( (*node_it)->taxid ); //delete from index
			reparent( iterator( node_it.node->parent ), node_it ); //move children
			node_it = erase( node_it ); //erase node
			delete deltaxon; //clear heap
		} else {
			++node_it;
		}
	}
}



void TaxonTree::setRankDistances( const std::vector< std::string >& ranks ) {
	//will only work if all possible ranks are contained in the given vector and are in the right order

	//build index structures
	int ranks_num = ranks.size();
	std::vector< const TTPString* > ranknames( ranks_num );
	std::map< const std::string*, int > rank2pos;
// 	std::cerr << "normalizing distances for ranks: ";
	for( int i = 0; i < ranks_num; ++i ) { //fill array with internal names in correct order
		const std::string& rankname = getRankInternal( ranks[ i ] );
// 		std::cerr << rankname << ", ";
		ranknames[ i ] = &rankname;
		rank2pos[ &rankname ] = i;
	}
// 	std::cerr << std::endl;

	// do it
	int running_index = 0;
	bool reset_running_index = false;
	iterator node_it = this->begin(); //root node
	TaxonNode* root_node = node_it.node;
	while ( ++node_it != this->end() ) { //skip root node
		if( (*node_it)->annotation->rank == node_it.node->parent->data->annotation->rank ) { //buggy NCBI taxonomy can have equal ranks as parents
			(*node_it)->root_pathlength = node_it.node->parent->data->root_pathlength;
			if( ! node_it.number_of_children() ) {
				reset_running_index = true;
			}
			continue;
		}
		(*node_it)->root_pathlength = node_it.node->parent->data->root_pathlength + 1; //keep distance info consistant
// 		std::cerr << "parent is: " << node_it.node->parent->data->annotation->name << "\t(" << node_it.node->parent->data->annotation->rank << ")\t" << node_it.node->parent->data->root_pathlength << std::endl;

		if( reset_running_index ) {
// 			std::cerr << "reset running index from: " << running_index;
			if( node_it.node->parent == root_node ) {
				running_index = 0;
			} else {
				running_index = rank2pos[ &node_it.node->parent->data->annotation->rank ] + 1;
			}
// 			std::cerr << " to " << running_index << std::endl;
			reset_running_index = false;
		}

		while( true ) { //in order to avoid unneccesary tests
			if( running_index == ranks_num ) {
// 				std::cerr << "current node: " << (*node_it)->annotation->name << "\t(" << (*node_it)->annotation->rank << ")\t" << (*node_it)->root_pathlength << std::endl;
// 				std::cerr << "fixing at this distance (a)!" << std::endl;
				reset_running_index = true;
				break;
			} else { //compare ranks
				if( &(*node_it)->annotation->rank == ranknames[ running_index ] ) {
// 					std::cerr << "current node: " << (*node_it)->annotation->name << "\t(" << (*node_it)->annotation->rank << ")\t" << (*node_it)->root_pathlength << std::endl;
// 					std::cerr << "fixing at this distance (b)!" << std::endl;
					// 						std::cerr << "this node has valid rank..." << std::endl;
					if( node_it.number_of_children() ) {
// 						std::cerr << "this node has still more children, not resetting index" << std::endl;
						++running_index;
					} else {
						reset_running_index = true;
					}
					break;
				} else {
					//insert virtual node by increasing the depth variable by 1
					(*node_it)->root_pathlength += 1;
					++running_index;
// 					std::cerr << "increasing distance to " << (*node_it)->root_pathlength << std::endl;
				}
			}
		}
	}
}



void TaxonTree::recalcNestedSetInfo() {
  //TODO: write
}
