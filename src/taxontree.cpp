#include "taxontree.hh"
#include <algorithm>
#include <vector>



TaxonTree::~TaxonTree() {
	//delete all taxons and annotations
	for( iterator node_it = this->begin(); node_it != this->end(); ++node_it ) {
		delete *node_it;
	}
}



void TaxonTree::addToIndex( unsigned int taxid , TaxonTree::Node* node ) {
	taxid2node_[ taxid ] = node;
}



void TaxonTree::recreateNodeIndex() {
	taxid2node_.clear();
	for( iterator node_it = this->begin(); node_it != this->end(); ++node_it ) {
		taxid2node_.insert( std::make_pair( (*node_it)->taxid, node_it.node ) );
	}
}



// constant in time as apposed to size(), I think
int TaxonTree::indexSize() const { //returns only real nodes (no dummies)
	return taxid2node_.size();
}



const std::string& TaxonTree::getRankInternal ( const std::string& rankname ) const {

	std::set< std::string >::const_iterator rank_it = ranks_.find( rankname );
	if( rank_it == ranks_.end() ) {
		return rank_not_found_;
	}
	return *rank_it;
}



const std::string& TaxonTree::insertRankInternal ( const std::string& rankname ) {
	return *ranks_.insert( rankname ).first;
}



void TaxonTree::deleteUnmarkedNodes() { //TODO: correct path_to_root: adjust at the end
	std::cerr << "deleting nodes that are not marked...";
	iterator node_it = ++( this->begin() ); //root node
	
	while( node_it != this->end() ) {
		Taxon* deltaxon = *node_it;
		if( ! deltaxon->mark_special ) {
			taxid2node_.erase( (*node_it)->taxid ); //delete from index
			reparent( iterator( node_it.node->parent ), node_it ); //move children
			node_it = erase( node_it ); //erase node
			delete deltaxon; //clear heap
		} else {
			++node_it;
		}
	}
	recalcDistToRoot( this->begin() ); //distances shrink
	std::cerr << " done" << std::endl;
}



void TaxonTree::setRankDistances( const std::vector< std::string >& ranklist ) {
	//will only work if all possible ranks are contained in the given vector and are in the right order

	//build index structures
	int ranks_num = ranklist.size();
	std::vector< const std::string* > ranknames( ranks_num );
	std::map< const std::string*, int > rank2pos;
// 	std::cerr << "normalizing distances for ranks: ";
	for( int i = 0; i < ranks_num; ++i ) { //fill array with internal names in correct order
		const std::string& rankname = getRankInternal( ranklist[ i ] );
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



void TaxonTree::setMaxDepth() {
	small_unsigned_int tmp = 0;
	for ( leaf_iterator it = begin_leaf(); it != end_leaf(); ++it ) {
		tmp = std::max( tmp, (*it)->root_pathlength );
	}
	max_depth_ = tmp;
}



void TaxonTree::recalcNestedSetInfo() {
  //TODO: write
}



void TaxonTree::recalcDistToRoot() {
	recalcDistToRoot( this->begin() ); //TODO: test if this works with nodes other than root
}


void TaxonTree::recalcDistToRoot( const iterator start ) {
	pre_order_iterator node_it( start );
	const small_unsigned_int start_depth = (*node_it)->root_pathlength;
	if( start_depth == 0 ) ++node_it;
	
// 	while( (*node_it)->root_pathlength < start_depth ) {
	while( node_it != this->end() ) {
		(*node_it)->root_pathlength = node_it.node->parent->data->root_pathlength + 1;
		++node_it;
	}
}