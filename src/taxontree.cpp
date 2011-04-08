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




TaxonNode* const TaxonomyInterface::getNode ( const unsigned int taxid ) {
	std::map< unsigned int, TaxonNode* >::iterator node_it = tax->taxid2node.find( taxid );
	if( node_it == tax->taxid2node.end() ) {
// 		std::cerr << "could not find node with TaxID " << taxid << " in taxonmy tree" << std::endl;
		return NULL;
	}
	return node_it->second;
}



TaxonNode* const TaxonomyInterface::getRoot() {
	return tax->begin().node;
}



const TTPString TaxonomyInterface::getRank ( const TaxonNode* node ) {
	return node->data->annotation->rank;
}



const TTPString TaxonomyInterface::getRank ( const unsigned int taxid ) {
	return tax->taxid2node[ taxid ]->data->annotation->rank;
}



const TTPString TaxonomyInterface::getName ( const TaxonNode* node ) {
	return node->data->annotation->name;
}



const TTPString TaxonomyInterface::getName ( const unsigned int taxid ) {
	return tax->taxid2node[ taxid ]->data->annotation->name;
}



bool TaxonomyInterface::isParentOf ( const TaxonNode* A, const TaxonNode* B ) {
	//if A's leftvalue is smaller than B's leftvalue and there is no range overlap this is just one comparison
	return A->data->rightvalue > B->data->leftvalue && A->data->leftvalue < B->data->leftvalue;
}



bool TaxonomyInterface::isParentOf ( const unsigned int A_taxid, const unsigned int B_taxid ) {
	TaxonNode* A = tax->taxid2node[ A_taxid ];
	TaxonNode* B = tax->taxid2node[ B_taxid ];
	return isParentOf( A, B );
}



TaxonNode* const TaxonomyInterface::getLCA ( const TaxonNode* A, const TaxonNode* B ) {
	unsigned int left_min = std::min( A->data->leftvalue, B->data->rightvalue );
	unsigned int right_max = std::max( A->data->rightvalue, B->data->rightvalue );

	TaxonNode* lca = (TaxonNode*) A;
	while(  lca->data->leftvalue > left_min || lca->data->rightvalue < right_max ) {
		lca = lca->parent;
	}

	return lca;
}



TaxonNode* const TaxonomyInterface::getLCA ( const unsigned int A_taxid, const unsigned int B_taxid ) { //TODO: handle taxid not found
	TaxonNode* A = tax->taxid2node[ A_taxid ];
	TaxonNode* B = tax->taxid2node[ B_taxid ];
	return getLCA( A, B );
}



TaxonNode* const TaxonomyInterface::getLCC ( TaxonNode* A, TaxonNode* B ) {
	if( isParentOf( B, A ) ) {
		return A;
	}
	return getLCA( A, B );
}



TaxonNode* const TaxonomyInterface::getLCC ( const unsigned int A_taxid, const unsigned int B_taxid ) { //TODO: handle taxid not found
	TaxonNode* A = tax->taxid2node[ A_taxid ];
	TaxonNode* B = tax->taxid2node[ B_taxid ];
	return getLCC( A, B );
}



std::pair< int, int > TaxonomyInterface::getPathLength ( const TaxonNode* A, const TaxonNode* B ) {
	if( A == B ) {
		return std::make_pair( 0, 0 );
	}
	if( isParentOf( A, B ) ) {
		return std::make_pair( getPathLengthToParent( B, A ), 0 );
	}
	if( isParentOf( B, A ) ) {
		return std::make_pair( 0, getPathLengthToParent( A, B ) );
	}

	//get distances to LCA
	TaxonNode* lca = getLCA( A, B );
	return std::make_pair( getPathLengthToParent( B, lca ), getPathLengthToParent( A, lca ) );
}



std::pair< int, int > TaxonomyInterface::getPathLength ( const unsigned int A_taxid, const unsigned int B_taxid ) {
	const TaxonNode* A = tax->taxid2node[ A_taxid ];
	const TaxonNode* B = tax->taxid2node[ B_taxid ];
	return getPathLength( A, B );
}



boost::tuple< int, int, int > TaxonomyInterface::getInterDistances( const TaxonNode* A, const TaxonNode* B ) {
	if( A == B ) {
		return boost::make_tuple( 0, A->data->root_pathlength, 0 );
	}
	if( isParentOf( A, B ) ) {
		return boost::make_tuple( 0, A->data->root_pathlength ,getPathLengthToParent( B, A ) );
	}
	if( isParentOf( B, A ) ) {
		return boost::make_tuple( getPathLengthToParent( A, B ), B->data->root_pathlength, 0 );
	}

	//get distances to LCA
	TaxonNode* lca = getLCA( A, B );
	return boost::make_tuple( getPathLengthToParent( B, lca ), lca->data->root_pathlength, getPathLengthToParent( A, lca ) );
}



boost::tuple< int, int, int > TaxonomyInterface::getInterDistances( const unsigned int A_taxid, const unsigned int B_taxid ) {
	const TaxonNode* A = tax->taxid2node[ A_taxid ];
	const TaxonNode* B = tax->taxid2node[ B_taxid ];
	return getInterDistances( A, B );
}



int TaxonomyInterface::getPathLengthToParent( const TaxonNode* A, const TaxonNode* B ) {
	return A->data->root_pathlength - B->data->root_pathlength;
}



int TaxonomyInterface::getPathLengthToParent( unsigned int A_taxid, unsigned int B_taxid ) {
	TaxonNode* A = tax->taxid2node[ A_taxid ];
	TaxonNode* B = tax->taxid2node[ B_taxid ];
	return getPathLengthToParent( A, B );
}


const std::string& TaxonomyInterface::getNameAtRank( TaxonNode* node, const std::string& rank ) {
	return getNameAtRank( node, &tax->getRankInternal( rank ) );
};



const std::string& TaxonomyInterface::getNameAtRank( TaxonNode* node, const std::string* internal_rank ) {
	TaxonNode* rn = getRoot();
	while( node != rn ) {
		if( node->data->annotation && ( &node->data->annotation->rank == internal_rank ) ) {
			return node->data->annotation->name;
		}
		node = node->parent;
	}
	return node->data->annotation->name;
};



TaxonNode* TaxonomyInterface::mapUnclassified( TaxonNode* node ) {
	const TaxonNode* root = getRoot();
	for(; node->data->is_unclassified && node != root; node=node->parent ) {};
	return node;
};



TaxonNode* TaxonomyInterface::mapUnclassified( unsigned int taxid ) {
	TaxonNode* tmp = getNode( taxid );
	if( tmp ) {
		return mapUnclassified( tmp );
	}
	return NULL;
};



void printPath( Taxonomy* tax, TaxonNode* node, TaxonNode* ancestor ) {
	TaxonomyInterface taxinter( tax );
	printPath( taxinter, node, ancestor );
}



void printPath( TaxonomyInterface& taxinter, TaxonNode* node, TaxonNode* ancestor ) {
	TaxonNode* tn = node;
	while( tn != ancestor ) {
		if( tn->data->mark_special && tn->data->annotation ) {
			std::cout << "name: " << tn->data->annotation->name <<
				", rank: " << tn->data->annotation->rank <<
				", dist_to_root: " << tn->data->root_pathlength <<
				std::endl;
			tn = tn->parent;
		}
	}
}
