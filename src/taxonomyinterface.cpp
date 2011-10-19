#include "taxonomyinterface.hh"
#include <algorithm>
#include <vector>



const TaxonNode* TaxonomyInterface::getNode ( const TaxonID taxid ) const {
	std::map< TaxonID, TaxonNode* >::const_iterator node_it = tax->taxid2node.find( taxid );
	if( node_it == tax->taxid2node.end() ) {
// 		std::cerr << "could not find node with TaxID " << taxid << " in taxonmy tree" << std::endl;
		return NULL;
	}
	return node_it->second;
}



const TaxonNode* TaxonomyInterface::getRoot() const {
	return tax->begin().node;
}



const TTPString TaxonomyInterface::getRank ( const TaxonNode* node ) const {
	return node->data->annotation->rank;
}



const TTPString TaxonomyInterface::getRank ( const TaxonID taxid ) const {
	return getNode( taxid )->data->annotation->rank;
}



const TTPString TaxonomyInterface::getName ( const TaxonNode* node ) const {
	return node->data->annotation->name;
}



const TTPString TaxonomyInterface::getName ( const TaxonID taxid ) const {
	return getNode( taxid )->data->annotation->name;
}



bool TaxonomyInterface::isParentOf ( const TaxonNode* A, const TaxonNode* B ) const {
	//if A's leftvalue is smaller than B's leftvalue and there is no range overlap this is just one comparison
	return A->data->rightvalue > B->data->leftvalue && A->data->leftvalue < B->data->leftvalue;
}



bool TaxonomyInterface::isParentOf ( const TaxonID A_taxid, const TaxonID B_taxid ) const {
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return isParentOf( A, B );
}



const TaxonNode* TaxonomyInterface::getLCA ( const TaxonNode* A, const TaxonNode* B ) const {
	large_unsigned_int left_min = std::min( A->data->leftvalue, B->data->rightvalue );
	large_unsigned_int right_max = std::max( A->data->rightvalue, B->data->rightvalue );

	const TaxonNode* lca = A;
	while(  lca->data->leftvalue > left_min || lca->data->rightvalue < right_max ) {
		lca = lca->parent;
	}

	return lca;
}



const TaxonNode* TaxonomyInterface::getLCA ( const TaxonID A_taxid, const TaxonID B_taxid ) const { //TODO: handle taxid not found
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return getLCA( A, B );
}



const TaxonNode* TaxonomyInterface::getLCC ( const TaxonNode* A, const TaxonNode* B ) const {
	if( isParentOf( B, A ) ) {
		return A;
	}
	return getLCA( A, B );
}



const TaxonNode* TaxonomyInterface::getLCC ( const TaxonID A_taxid, const TaxonID B_taxid ) const { //TODO: handle taxid not found
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return getLCC( A, B );
}



std::pair< small_unsigned_int, small_unsigned_int > TaxonomyInterface::getPathLength ( const TaxonNode* A, const TaxonNode* B ) const {
	if( A == B ) {
		return std::make_pair< small_unsigned_int, small_unsigned_int >( 0, 0 );
	}
	if( isParentOf( A, B ) ) {
		return std::make_pair< small_unsigned_int, small_unsigned_int >( getPathLengthToParent( B, A ), 0 );
	}
	if( isParentOf( B, A ) ) {
		return std::make_pair< small_unsigned_int, small_unsigned_int >( 0, getPathLengthToParent( A, B ) );
	}

	//get distances to LCA
	const TaxonNode* lca = getLCA( A, B );
	return std::make_pair( getPathLengthToParent( B, lca ), getPathLengthToParent( A, lca ) );
}



std::pair< small_unsigned_int, small_unsigned_int > TaxonomyInterface::getPathLength ( const TaxonID A_taxid, const TaxonID B_taxid ) const {
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return getPathLength( A, B );
}



boost::tuple< small_unsigned_int, small_unsigned_int, small_unsigned_int > TaxonomyInterface::getInterDistances( const TaxonNode* A, const TaxonNode* B ) const {
	if( A == B ) {
		return boost::make_tuple( 0, A->data->root_pathlength, 0 );
	}
	if( isParentOf( B, A ) ) {
		return boost::make_tuple( getPathLengthToParent( A, B ), B->data->root_pathlength, 0 );
	}
	if( isParentOf( A, B ) ) {
		return boost::make_tuple( 0, A->data->root_pathlength ,getPathLengthToParent( B, A ) );
	}

	//get distances to LCA
	const TaxonNode* lca = getLCA( A, B );
	return boost::make_tuple( getPathLengthToParent( A, lca ), lca->data->root_pathlength, getPathLengthToParent( B, lca ) );
}



boost::tuple< small_unsigned_int, small_unsigned_int, small_unsigned_int > TaxonomyInterface::getInterDistances( const TaxonID A_taxid, const TaxonID B_taxid ) const {
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return getInterDistances( A, B );
}



small_unsigned_int TaxonomyInterface::getPathLengthToParent( const TaxonNode* A, const TaxonNode* B ) const {
	return A->data->root_pathlength - B->data->root_pathlength;
}



small_unsigned_int TaxonomyInterface::getPathLengthToParent( TaxonID A_taxid, TaxonID B_taxid ) const {
	const TaxonNode* A = getNode(  A_taxid );
	const TaxonNode* B = getNode(  B_taxid );
	return getPathLengthToParent( A, B );
}


const std::string& TaxonomyInterface::getNameAtRank( const TaxonNode* node, const std::string& rank ) const {
	return getNameAtRank( node, &tax->getRankInternal( rank ) );
}



const std::string& TaxonomyInterface::getNameAtRank( const TaxonNode* node, const std::string* internal_rank ) const {
	const TaxonNode* rn = getRoot();
	while( node != rn ) {
		if( node->data->annotation && ( &node->data->annotation->rank == internal_rank ) ) {
			return node->data->annotation->name;
		}
		node = node->parent;
	}
	return node->data->annotation->name;
}



const TaxonNode* TaxonomyInterface::mapUnclassified( const TaxonNode* node ) const {
	const TaxonNode* root = getRoot();
	for(; node->data->is_unclassified && node != root; node=node->parent ) {};
	return node;
}



const TaxonNode* TaxonomyInterface::mapUnclassified( TaxonID taxid ) const {
	const TaxonNode* tmp = getNode( taxid );
	if( tmp ) {
		return mapUnclassified( tmp );
	}
	return NULL;
}



TaxonTree::PathUpIterator TaxonomyInterface::traverseUp( const TaxonNode* node ) const {
	return TaxonTree::PathUpIterator( node );
}



bool TaxonomyInterface::isLeaf(const TaxonNode* node) const {
	return ! node->first_child;
}



// void printPath( Taxonomy* tax, TaxonNode* node, TaxonNode* ancestor ) {
// 	TaxonomyInterface taxinter( tax );
// 	printPath( taxinter, node, ancestor );
// }
// 
// 
// 
// void printPath( TaxonomyInterface& taxinter, TaxonNode* node, TaxonNode* ancestor ) {
// 	TaxonNode* tn = node;
// 	while( tn != ancestor ) {
// 		if( tn->data->mark_special && tn->data->annotation ) {
// 			std::cout << "name: " << tn->data->annotation->name <<
// 				", rank: " << tn->data->annotation->rank <<
// 				", dist_to_root: " << tn->data->root_pathlength <<
// 				std::endl;
// 			tn = tn->parent;
// 		}
// 	}
// }
