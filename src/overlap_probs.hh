#ifndef INCLUDED_tcrdist_cxx_overlap_probs_HH
#define INCLUDED_tcrdist_cxx_overlap_probs_HH

#include "misc.hh"
#include <boost/math/distributions/hypergeometric.hpp>

inline
Real
product_cdf( Real const x )
{
	return x - x*log(x);
}

inline
Real
compute_overlap_pvalue(
	Size const overlap,
	Size const count1,
	Size const count2,
	Size const total
)
{
	using namespace boost::math;
	// whats the smallest possible overlap? max( 0, count1+count2-total )
	if ( overlap==0 || count1==total || count2==total || count1+count2 >= total+overlap ) return 1.0;
	hypergeometric_distribution<> hgd( count1, count2, total );
	return cdf( complement( hgd, overlap-1 ) ); // need overlap-1 to be a valid overlap value!
}


inline // silly?
Real
compute_overlap_pvalue_with_bias(
	bools const & occs1,
	bools const & occs2,
	Reals const & subject_bias_factor
)
{
	Real pval(1.);

	Size table[2][2];

	table[0][0] = table[0][1] = table[1][0] = table[1][1] = 0;

	Size const total( occs1.size() );
	runtime_assert( occs2.size() == total );
	runtime_assert( subject_bias_factor.size() == total );

	for ( Size i=0; i< total; ++i ) { ++table[ occs1[i] ][ occs2[i] ]; }

	Size const overlap( table[1][1] ), total1( table[1][0] + table[1][1] ), total2( table[0][1] + table[1][1] );

	if ( overlap > 0 && total1 < total && total2 < total ) {
		Real bias(1.0);
		for ( Size i=0; i< total; ++i ) {
			if ( occs1[i] && occs2[i] ) bias *= subject_bias_factor[i];
		}
		pval = min( 1.0, bias * compute_overlap_pvalue( overlap, total1, total2, total ) );
	}

	return pval;
}

#endif
