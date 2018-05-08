#ifndef INCLUDED_pubtcrs_cluster_scores_HH
#define INCLUDED_pubtcrs_cluster_scores_HH

#include <boost/math/distributions/binomial.hpp>

#include "misc.hh"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Size
choose_random_Size(
	Reals const & probs
)
{
	//using namespace boost::math;
	Real const f( uniform());

	Real total(0);

	for ( Size i=0; i<probs.size(); ++i ) {
		total += probs[i];
		if ( f <= total ) return i;
	}
	cerr << "ERROR choose_random_Size f out of bounds: " << f << ' ' << total << ' ' << probs.size() << endl;
	return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
choose_random_subjects_with_bias(
	Size const num_to_choose,
	Reals const & sampling_bias,
	bools & chosen
)
{
	Size const total_subjects( sampling_bias.size() );
	runtime_assert( num_to_choose <= total_subjects );

	chosen.resize( total_subjects );
	std::fill( chosen.begin(), chosen.end(), false );

	{ // sanity check
		Real total_bias(0);
		foreach_( Real f, sampling_bias ) total_bias += f;
		runtime_assert( fabs( 1.0 - total_bias )<1e-2 );
	}

	Size num_chosen(0);
	while ( num_chosen < num_to_choose ) {
		Size const s( choose_random_Size( sampling_bias ) );
		if ( !chosen[s] ) {
			chosen[s] = true;
			++num_chosen;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// both lists should sum to the same thing: the sum of the subject counts for all the tcrs in the cluster
//
//
Real // return the best score
score_l1_versus_l2(
	Reals const & l1,  // should be sorted
	Reals const & l2,  // -- ditto --
	Size & min_switchpoint // for i <= min_switchpoint, score += l1[i] - l2[i];
) //                            i >  min_switchpoint, score += l2[i] - l2[i];
{
	Real const epsilon( 1e-4 );

	runtime_assert( l1.front() >= l1.back()-epsilon ); // silly check for sorted
	runtime_assert( l2.front() >= l2.back()-epsilon ); // ditto

	Size const nsubjects( l2.size() );

	runtime_assert( l1.size() == nsubjects );

	Reals total_l1_minus_l2;

	Real total(0);
	for ( Size i=0; i< nsubjects; ++i ) {
		total += l1[i] - l2[i];
		total_l1_minus_l2.push_back( total );
	}

	runtime_assert( fabs( total_l1_minus_l2[ nsubjects-1 ] - 0.0 )<epsilon );

	Real best_score( 0.0 ), score;
	min_switchpoint = nsubjects-1; // this should give a score of 0.0

	total=0.;

	for ( Size i=nsubjects-1; i > 0; --i ) {
		total += l2[i] - l1[i];
		// this is the score if min_switchpoint were i-1:
		score = total_l1_minus_l2[ i-1 ] + total;
		if ( score >= best_score - epsilon ) { // choose smaller switchpoints with equal scores
			best_score = score;
			min_switchpoint = i-1;
		}
	}

	return best_score;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real // return diff_Zscore
analyze_cluster_occurrences(
	Sizes const & tcrs,
	vector< bools > const & all_occs, // could be just intra-HLA occs
	Sizes const & subject_indices, // for output, eg might be: ( hla_positive_subjects.find(hla)->second ) );
	Reals const & subject_sampling_bias,
	Size const nrep,
	string const tag, // empty for classic output
	Size & upper_tcr_count,
	Size & lower_tcr_count,
	ostream & out
)
{
	Size const total_subjects( all_occs.front().size() );
	runtime_assert( subject_sampling_bias.size() == total_subjects );
	runtime_assert( subject_indices.size() == total_subjects );

	Reals avg_counts( total_subjects, 0. );
	Real diff_score_mean, diff_score_sdev;

	bools chosen( total_subjects ); // re-use this array
	Sizes counts( total_subjects );
	Reals realcounts( total_subjects );


	Reals rand_diff_scores;

	// go through this loop 2*nrep times, building up the average the first time, then comparing the counts
	// 	to the average the second time, so we can get a Z-score for the observed counts
	for ( Size rep=1; rep<= 2*nrep; ++rep ) {
		std::fill( counts.begin(), counts.end(), Size(0) );

		// randomly choose tcr occurrences
		foreach_( Size t, tcrs ) {
			Size num_subjects_this_tcr(0);
			foreach_( bool occ, all_occs[t] ) { if ( occ ) ++num_subjects_this_tcr; }
			if ( num_subjects_this_tcr ) {
				choose_random_subjects_with_bias( num_subjects_this_tcr, subject_sampling_bias, chosen );
				for ( Size i=0; i< total_subjects; ++i ) {
					if ( chosen[i] ) ++counts[i];
				}
			}
		}

		sort( counts.begin(), counts.end() ); // increasing order
		std::reverse( counts.begin(), counts.end() ); // now decreasing order

		if ( rep <= nrep ) {
			for ( Size i=0; i< total_subjects; ++i ) avg_counts[i] += counts[i];

			if ( rep == nrep ) { // normalize avg scores, write them out
				for ( Size i=0; i< total_subjects; ++i ) avg_counts[i] /= nrep;
				out << "avg_subject_tcr_counts" << tag << ": total_subjects= " << total_subjects;
				for ( Size i=0; i< total_subjects; ++i ) out << ' ' << F(9,3,avg_counts[i]);
				out << '\n';
			}
		} else {
			/// score our counts against avg_counts
			for ( Size i=0; i< total_subjects; ++i ) realcounts[i] = Real(counts[i]);

			Size tmp_switchpoint;
			Real pos_score( score_l1_versus_l2( realcounts, avg_counts, tmp_switchpoint ) ),
				neg_score( score_l1_versus_l2( avg_counts, realcounts, tmp_switchpoint ) );

			if ( pos_score > neg_score ) {
				rand_diff_scores.push_back( pos_score );
			} else {
				rand_diff_scores.push_back( -1 * neg_score );
			}
			runtime_assert( rand_diff_scores.size() == rep-nrep );
		}
	} // rep=1, 2*nrep
	runtime_assert( rand_diff_scores.size() == nrep );
	get_mean_sdev( rand_diff_scores, diff_score_mean, diff_score_sdev );


	Sizes subject_tcr_counts( total_subjects, 0 );
	for ( Size i=0; i< tcrs.size(); ++i ) {
		bools const & i_occs( all_occs[ tcrs[i] ] );
		for ( Size k=0; k< total_subjects; ++k ) {
			if ( i_occs[k] ) ++subject_tcr_counts[k];
		}
	}
	SizePairs countslist;
	Size num_nonzero(0);
	for ( Size k=0; k< total_subjects; ++k ) {
		if ( subject_tcr_counts[k]>0 ) ++num_nonzero;
		countslist.push_back( make_pair( subject_tcr_counts[k], k ) );
	}
	sort( countslist.begin(), countslist.end () );
	reverse( countslist.begin(), countslist.end () );

	//Reals realcounts( total_subjects, 0.0 );
	for ( Size i=0; i< total_subjects; ++i ) {
		realcounts[i] = countslist[i].first;
	}

	Size pos_switchpoint, neg_switchpoint;
	Real pos_score( score_l1_versus_l2( realcounts, avg_counts, pos_switchpoint ) ),
		neg_score( score_l1_versus_l2( avg_counts, realcounts, neg_switchpoint ) ),
		diff_score( pos_score>neg_score ? pos_score : -1*neg_score ),
		diff_Zscore( ( diff_score - diff_score_mean )/ diff_score_sdev );

	upper_tcr_count = countslist[ pos_switchpoint ].first;
	lower_tcr_count = ( pos_switchpoint < total_subjects ? countslist[ pos_switchpoint+1 ].first : 0 );


	out << "subject_tcr_counts" << tag << ": total_subjects= " << total_subjects <<
		" num_nonzero_tcr_counts: " << num_nonzero <<
		" D_CO: " << F(9,3,diff_score) <<
		" Z_CO: " << F(9,3,diff_Zscore) <<
		" D_CO_rand_mean: " << F(9,3,diff_score_mean) <<
		" D_CO_rand_sdev: " << F(9,3,diff_score_sdev) <<
		" pos_switchpoint: " << pos_switchpoint <<
		" upper_tcr_count: " << upper_tcr_count <<
		" lower_tcr_count: " << lower_tcr_count;
	// Sizes const & hla_subs( hla_positive_subjects.find(hla)->second );
	// runtime_assert( hla_subs.size() == total_subjects );
	for ( Size k=0; k< total_subjects; ++k ) {
		if ( countslist[k].first == 0 ) break; // don't show the zeros
		Size const subject( subject_indices[ countslist[k].second ] );
		out << ' ' << subject << ':' << countslist[k].first;
			//<< ':' << subject_sampling_bias[ countslist[k].second ]; // now show the bias...
	}
	out << '\n';

	return diff_Zscore;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// returns 1 if we don't have info for this (similarity_mode,nbr_pval_threshold)
//
//
// probability of seeing an equal or greater number of neighbors for a random TCR, given num_clustered_tcrs
//
Real
compute_cluster_size_score(
	Size const num_center_nbrs, // max nbr number for cluster member
	Size const num_clustered_tcrs,
	Size const similarity_mode,
	Real const nbr_pval_threshold
)
{
	using namespace boost::math;
	static bool init( false );
	static map< pair<Size,Real>, Reals > precomputed_nbr_rates;

	if ( !init ) {
		init = true;
		string const filename( misc::dbdir + "nbr_rates_for_P_size_calc.txt" );
		strings lines;
		read_lines_from_file( filename, lines );
		foreach_( string line, lines ) {
			istringstream l(line);
			string tmp;
			Size l_sim_mode;
			Real l_nbr_pval_threshold;
			l >> tmp >> l_sim_mode >> tmp >> l_nbr_pval_threshold >> tmp;
			Reals nbr_rates;
			while ( !l.fail() ) {
				Real rate;
				l >> rate;
				if ( !l.fail() ) nbr_rates.push_back( rate );
			}
			cout << "Read " << nbr_rates.size() << " precomputed_nbr_rates for similarity_mode: " << l_sim_mode <<
				" and nbr_pval_threshold: " << l_nbr_pval_threshold << endl;
			precomputed_nbr_rates[ make_pair( l_sim_mode, l_nbr_pval_threshold ) ] = nbr_rates;
		}
	}

	Reals nbr_rates;

	if ( similarity_mode == 1 ) {
		nbr_rates.push_back( nbr_pval_threshold );
	} else { // look for info in the file
		for ( map< pair<Size,Real>, Reals >::const_iterator it= precomputed_nbr_rates.begin();
					it!= precomputed_nbr_rates.end(); ++it ) {
			if ( it->first.first == similarity_mode && fabs( log( it->first.second ) - log( nbr_pval_threshold ) )<.1 ) {
				// match
				nbr_rates = it->second;
			}
		}
	}

	if ( nbr_rates.empty() ) {
		cerr << "Failed to find precomputed information for computing cluster size score:: " <<
			" similarity_mode: " << similarity_mode << " nbr_pval_threshold: " << nbr_pval_threshold << endl;
		return 1.0;
	}

	//
	Real P_size(0); // odds of seeing an equal or greater neighbor number, given the total number of tcrs clustered

	Real const wt( 1.0/ nbr_rates.size() );
	foreach_( Real rate, nbr_rates ) {
		P_size += wt * cdf( complement( binomial( num_clustered_tcrs-1, rate ), num_center_nbrs-1 ) );
	}

	return P_size;
}

#endif
