#ifndef INCLUDED_pubtcrs_dbscan_HH
#define INCLUDED_pubtcrs_dbscan_HH

#include "types.hh"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
dbscan_cluster(
	Size const min_core_nbrcount,
	Size const min_cluster_size,
	vector< Sizes > const & all_nbrs,
	Sizes & centers, // the cluster member with the largest number of nbrs
	vector< Sizes > & all_members
)
{

	Size const N( all_nbrs.size() );
	Size const UNK( N );
	Size const BORDER( N+1 );
	Sizes all_clusterno( N );

	// identify 'core' points, sort by number of nbrs
	SizePairs l;
	for ( Size i=0; i< N; ++i ) {
		if ( all_nbrs[i].size() >= min_core_nbrcount ) {
			all_clusterno[i] = UNK; // core point, initially assigned to cluster UNK
			l.push_back( make_pair( all_nbrs[i].size(), i ) );
		} else {
			all_clusterno[i] = BORDER; // non-core point, initially assigned to cluster BORDER
		}
	}
	sort( l.begin(), l.end() );
	reverse( l.begin(), l.end() ); // now in decreasing order
	Sizes core_pts;
	foreach_( SizePair p,l ) {core_pts.push_back(p.second);}

	int clusterno(-1);
	while ( true ) {
		bool all_done( true );
		foreach_( Size center, core_pts ) {
			if ( all_clusterno[center] == UNK ) {
				all_done = false;
				// start with center, extend to all connected core points
				++clusterno;
				Sizes cluster_core;
				Sizes newl, oldl;
				newl.push_back( center );
				all_clusterno[ center ] = clusterno;
				cluster_core.push_back( center );
				while ( newl.size() ){
					oldl.swap(newl);
					runtime_assert( newl.empty() );
					foreach_( Size c, oldl ) {
						foreach_( Size n, all_nbrs[c] ) {
							if ( all_clusterno[n] == UNK ) { // unclassified core pt
								newl.push_back( n );
								all_clusterno[n] = clusterno;
								cluster_core.push_back(n);
							}
						}
					}
					oldl.clear();
				}
				// now extend the cluster to include non-core border points
				Sizes members( cluster_core );
				foreach_( Size c, cluster_core ) {
					foreach_( Size n, all_nbrs[c] ) {
						if ( all_clusterno[n] == BORDER ) { // unclassified border point
							all_clusterno[n] = clusterno;
							members.push_back(n);
						}
					}
				}

				if ( members.size() >= min_cluster_size ) {
					centers.push_back( center );
					all_members.push_back( members );
				}
				break;
			}
		}
		if ( all_done ) break;
	}

}

#endif
