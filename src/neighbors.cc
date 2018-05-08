#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"
#include "overlap_probs.hh"

#include <tclap/CmdLine.h>

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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Size
choose_random_Size(
	Reals const & probs
)
{
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

	try {

		TCLAP::CmdLine cmd("Find neighbors (e.g., for clustering) by comparing TCR occurrences "
			"in matrix1 to those in matrix2. If matrix2 is not provided, nbrs will be found within matrix1. "
			"If --feature_mask argument is provided, the subjects used for clustering will be the intersection "
			"of the positive subject sets for all the features in the provided file. Otherwise all the subjects will "
			"be used.", ' ', "0.1");

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"./db/","path",cmd);

 		TCLAP::ValueArg<int> min_cluster_size_arg("z","min_cluster_size","Minimum cluster size "
			" for DBSCAN clustering.",
			false, -1, "unsigned integer", cmd );
 		TCLAP::ValueArg<string> min_core_nbrcounts_file_arg("n","min_core_nbrcounts_file","File mapping from number "
			"of clustered TCRs to min_core_nbrcount values.",
			false, "", "filename", cmd );
 		TCLAP::ValueArg<int> min_core_nbrcount_arg("m","min_core_nbrcount","Minimum neighbor-number for a TCR "
			"to be considered a 'core' point for DBSCAN clustering (aka N_core).",
			false, -1, "unsigned integer", cmd );
		TCLAP::SwitchArg cluster_arg("c","cluster","Use nbrs to do DBSCAN clustering", cmd, false);
 		TCLAP::ValueArg<string> feature_mask_arg("f","feature_mask","For subsetting the subjects to be used",
			false,"","filename",cmd);
 		TCLAP::ValueArg<string> subject_bias_arg("b","subject_bias","For heuristic p-value adjustment",
			true,"","filename",cmd);
 		TCLAP::ValueArg<Real> nbr_pval_threshold_arg("p","nbr_pval_threshold",
			"P-value threshold for defining neighbors (aka T_sim)",
			true, Real(1e-4), "double", cmd );
 		TCLAP::ValueArg<string> matrix2_arg("j","matrix2","File containing matrix of TCR occurrences, "
			"to be compared with matrix1",false,"","filename",cmd);
 		TCLAP::ValueArg<string> matrix1_arg("i","matrix1","File containing matrix of TCR occurrences, "
			"to be compared with matrix2 (or self if matrix2 is not provided)",true,"","filename",cmd);
		Sizes allowed_similarity_modes( make_vector(Size(1),Size(2)) );
		TCLAP::ValuesConstraint<Size> similarity_mode_arg_cst(allowed_similarity_modes);
 		TCLAP::ValueArg<Size> similarity_mode_arg("s","similarity_mode","Mode switch for similarity calculation.\n"
			" 1 == co-occurrence only\n 2 == co-occurrence and TCRdist",
			false, 1, &similarity_mode_arg_cst, cmd );

		cmd.parse( argc, argv );

		set_dbdir( database_arg.getValue() );

		TCRdistCalculator tcrdist; // after set_dbdir

		Size const similarity_mode( similarity_mode_arg.getValue() );
		string const matrix_file1( matrix1_arg.getValue() );
		string const matrix_file2( matrix2_arg.getValue() );
		string const subject_bias_file( subject_bias_arg.getValue() );
		string const feature_file( feature_mask_arg.getValue() );
		Real const nbr_pval_threshold( nbr_pval_threshold_arg.getValue() );

		Real const min_tcrdist_pval( 1.2e-4 );
		bool const do_clustering( cluster_arg.getValue() );
		bool const write_nbrs_to_stdout( !do_clustering );

		if ( do_clustering && ( min_core_nbrcount_arg.getValue() < 0 && min_core_nbrcounts_file_arg.getValue().empty() ) ) {
			cout << "Please pass --min_core_nbrcount or --min_core_nbrcounts_file for clustering" << endl;
			exit(1);
		}
		if ( do_clustering && min_cluster_size_arg.getValue() < 0 ) {
			cout << "Please pass --min_cluster_size for clustering" << endl;
			exit(1);
		}
		if ( do_clustering && !matrix_file2.empty() ) {
			cout << "No need for the --matrix2 arg when clustering (it only uses the matrix1 TCRs)" << endl;
			exit(1);
		}

		// setup
		Sizes subjects_subset; // empty is signal to use all subjects
		if ( !feature_file.empty() ) {
			Features fl;
			read_features( feature_file, fl );
			subjects_subset = fl.front().poslist;
			if ( fl.size()>1 ) {
				for ( Size j=1; j< fl.size(); ++j ) {
					Sizes const other_poslist( fl[j].poslist );
					Sizes delme;
					foreach_( Size s, subjects_subset ) {
						if ( !has_element( s, other_poslist ) ) delme.push_back(s);
					}
					foreach_( Size s, delme ) {
						subjects_subset.erase( find( subjects_subset.begin(), subjects_subset.end(), s ) );
					}
				}
			}
			cout << "final_subjects_subset: " << subjects_subset.size() << endl;
			runtime_assert( !subjects_subset.empty() );
		}


		bool const self_nbrs( matrix_file2.empty() );

		vector< bools > all_occs1, all_occs2;
		strings all_tcrs1, all_tcrs2;
		read_matrix_file( matrix_file1, all_occs1, all_tcrs1 );

		if ( !self_nbrs ) {
			read_matrix_file( matrix_file2, all_occs2, all_tcrs2 );
		} else {
			all_tcrs2 = all_tcrs1; // but don't make a copy of all_occs1
		}

		Reals subject_bias_factors;
		read_floats_from_file( subject_bias_file, subject_bias_factors ); // whitespace separated
		runtime_assert( subject_bias_factors.size() == all_occs1.front().size() );

		/// may have to subset the subjects
		if ( !subjects_subset.empty() ) {
			sort( subjects_subset.begin(), subjects_subset.end() );
			foreach_( bools & occs, all_occs1 ) {
				bools new_occs;
				foreach_( Size s, subjects_subset ) new_occs.push_back( occs[s] );
				occs.swap( new_occs );
			}

			if ( !self_nbrs ) {
				foreach_( bools & occs, all_occs2 ) {
					bools new_occs;
					foreach_( Size s, subjects_subset ) new_occs.push_back( occs[s] );
					occs.swap( new_occs );
				}
			}


			Reals new_bias_factors;
			foreach_( Size s, subjects_subset ) new_bias_factors.push_back( subject_bias_factors[s] );
			subject_bias_factors.swap( new_bias_factors );
		} else {
			// fake subset info for clustering
			for ( Size i=0; i< subject_bias_factors.size(); ++i ) {
				subjects_subset.push_back( i );
			}
		}

		Size const total_subjects( all_occs1.front().size() );
		runtime_assert( subject_bias_factors.size() == total_subjects );
		runtime_assert( total_subjects == subjects_subset.size() );

		/// normalize subject_bias_factors
		{
			Real avg(0);
			foreach_( Real f, subject_bias_factors ) {
				runtime_assert( f>=0 );
				avg += f;
			}
			avg /= total_subjects;
			foreach_( Real &f, subject_bias_factors ) {
				f /= avg;
			}
		}

		vector< DistanceTCR_f > all_dtcrs1, all_dtcrs2;

		foreach_( string tcr, all_tcrs1 ) {
			all_dtcrs1.push_back( tcrdist.create_distance_tcr_f( tcr ) );
		}
		foreach_( string tcr, all_tcrs2 ) {
			all_dtcrs2.push_back( tcrdist.create_distance_tcr_f( tcr ) );
		}

		vector<Sizes> all_nbrs1;
		if ( do_clustering ) all_nbrs1.resize( all_tcrs1.size() );

		for ( Size ii=0; ii< all_tcrs1.size(); ++ii ) {
			DistanceTCR_f const & ii_dtcr( all_dtcrs1[ii] );
			bools const & ii_occs( all_occs1[ii] );
			for ( Size jj=0; jj< all_tcrs2.size(); ++jj ) {
				if ( self_nbrs && jj<=ii ) continue;

				bools const & jj_occs( self_nbrs ? all_occs1[jj] : all_occs2[jj] );

				Real const co_pval( compute_overlap_pvalue_with_bias( ii_occs, jj_occs, subject_bias_factors ) );
				Real combined_pval( co_pval ), tcrd(0), tcrdist_pval(1);
				if ( similarity_mode == 2 ) { // include TCRdist similarity
					tcrd = tcrdist( ii_dtcr, all_dtcrs2[jj] );
					tcrdist_pval = tcrdist.prob_equal_or_smaller_tcrdist( tcrd );
					combined_pval = product_cdf( max( min_tcrdist_pval, tcrdist_pval ) * co_pval );
				}

				if ( combined_pval <= nbr_pval_threshold ) {
					if ( write_nbrs_to_stdout ) {
						cout << "nbrs: " << ii << ' ' << jj <<
							" tcr1: " << all_tcrs1[ii] <<
							" tcr2: " << all_tcrs2[jj] <<
							" co_pval: " << co_pval;
						if ( similarity_mode == 1 ) cout << '\n';
						else cout << " tcrdist: " << tcrd << " combined_pval: " << combined_pval << '\n';
					}
					if ( do_clustering ) {
						all_nbrs1[ii].push_back( jj );
						all_nbrs1[jj].push_back( ii );
					}
				}
			}
		}


		if ( do_clustering ) {

			Size const min_cluster_size( min_cluster_size_arg.getValue() );

			/// There are two ways to set min_core_nbrcount:
			///
			///  --min_core_nbrcount <value>
			///        or
			///  --min_core_nbrcounts_file <filename>
			///
			///  In the latter case the file has a mapping from number of clustered TCRs to min_core_nbrcount values
			///

			Size min_core_nbrcount(0);
			if ( min_core_nbrcount_arg.getValue() >= 0 ) { // passed on cmdline
				min_core_nbrcount = min_core_nbrcount_arg.getValue();
			} else { // provided by a file, works for different numbers of clustered tcrs...
				string const filename( min_core_nbrcounts_file_arg.getValue() );
				Size const num_clustered_tcrs( all_tcrs1.size() );
				strings lines;
				read_lines_from_file( filename, lines );
				foreach_( string & line, lines ) {
					istringstream l(line);
					Size l_num_tcrs, l_min_core_nbrcount;
					l >> l_num_tcrs >> l_min_core_nbrcount;
					if ( l.fail() ) {
						cerr << "parse error " << filename << ' ' << line << endl;
						exit(1);
					}
					if ( l_num_tcrs == num_clustered_tcrs ) {
						min_core_nbrcount = l_min_core_nbrcount;
						break;
					}
				}
				if ( !min_core_nbrcount ) {
					cerr << "Failed to find num_clustered_tcrs (" << num_clustered_tcrs << ") in file " << filename << endl;
					exit(1);
				}
			}
			runtime_assert( min_core_nbrcount );


			/// get rid of the '1' suffix since there is only one set of tcrs/occs anyhow:
			vector< Sizes > const & all_nbrs( all_nbrs1 );
			vector< bools > const & all_occs( all_occs1 );
			strings const & all_tcrs( all_tcrs1 );

			cout << "Clustering " << all_tcrs.size() << " TCRs " <<
				" min_core_nbrcount= " << min_core_nbrcount <<
				" nbr_pval_threshold= " << nbr_pval_threshold << endl;

			Sizes centers;
			vector< Sizes > all_members;
			dbscan_cluster( min_core_nbrcount, min_cluster_size, all_nbrs, centers, all_members );

			cout << "Found " << centers.size() << " DBSCAN clusters." << endl;

			/// fill probability array for choosing each subject during resampling (for Z_CO score calc)
			Reals subject_sampling_bias;
			Real tot(0);
			foreach_( Real f, subject_bias_factors ) {
				subject_sampling_bias.push_back( f/total_subjects );
				tot += f;
			}
			runtime_assert( fabs( tot - total_subjects )<1e-2 ); // sanity check on normalization

			for ( Size ic=0; ic< centers.size(); ++ic ) {
				Sizes const & members( all_members[ic] );
				Size const center( centers[ic] );
				// for output, cluster numbers are 1-indexed:
				cout << "cluster: " << ic+1 << " size: " << members.size() << " center: " << all_tcrs[ center ] <<
					" num_center_nbrs: " << all_nbrs[center].size() << " members:";
				foreach_( Size m, members ) {
					cout << ' ' << all_tcrs[ m ];
				}
				cout << endl;

				// compute Z_CO score
				Size const Z_CO_random_repeats(1000);
				Size upper_tcr_count, lower_tcr_count;
				analyze_cluster_occurrences( members, all_occs, subjects_subset, subject_sampling_bias, Z_CO_random_repeats, "",
					upper_tcr_count, lower_tcr_count, cout );
			}
		}


	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
