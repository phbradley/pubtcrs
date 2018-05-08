#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"
#include "dbscan.hh"
#include "cluster_scores.hh"
#include "overlap_probs.hh"

#include <tclap/CmdLine.h>


int main(int argc, char** argv)
{

	try {

		TCLAP::CmdLine cmd("Find neighbors (e.g., for clustering) by comparing TCR occurrences "
			"in matrix1 to those in matrix2. If matrix2 is not provided, nbrs will be found within matrix1. "
			"If --feature_mask argument is provided, the subjects used for neighbors/clustering will be the intersection "
			"of the positive subject sets for all the features in the provided file. Otherwise all the subjects will "
			"be used. If the --cluster flag is used, DBSCAN clustering will be performed using the identified "
			"neighbor relationships.", ' ', "0.1");

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"./db/","path",cmd);

 		TCLAP::ValueArg<Size> min_cluster_size_arg("z","min_cluster_size","Minimum cluster size "
			" for DBSCAN clustering.",
			false, 0, "unsigned integer", cmd );

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
			false, 0.0, "double", cmd );

 		TCLAP::ValueArg<string> matrix2_arg("j","matrix2","File containing matrix of TCR occurrences, "
			"to be compared with matrix1",false,"","filename",cmd);

 		TCLAP::ValueArg<string> matrix1_arg("i","matrix1","File containing matrix of TCR occurrences, "
			"to be compared with matrix2 (or self if matrix2 is not provided)",true,"","filename",cmd);

		Sizes allowed_similarity_modes( make_vector(Size(1),Size(2),Size(3)) );
		TCLAP::ValuesConstraint<Size> similarity_mode_arg_cst(allowed_similarity_modes);
 		TCLAP::ValueArg<Size> similarity_mode_arg("s","similarity_mode","Mode switch for similarity calculation.\n"
			" 1 == co-occurrence only\n 2 == co-occurrence and TCRdist\n 3 == TCRdist only",
			true, 0, &similarity_mode_arg_cst, cmd );

		cmd.parse( argc, argv );

		set_dbdir( database_arg.getValue() );

		TCRdistCalculator tcrdist; // after set_dbdir

		Size const similarity_mode( similarity_mode_arg.getValue() );
		string const matrix_file1( matrix1_arg.getValue() );
		string const matrix_file2( matrix2_arg.getValue() );
		string const subject_bias_file( subject_bias_arg.getValue() );
		string const feature_file( feature_mask_arg.getValue() );
		Real nbr_pval_threshold( nbr_pval_threshold_arg.getValue() );
		if ( nbr_pval_threshold == 0 ) {
			/// try to set a sensible value
			if ( similarity_mode == 1 || similarity_mode == 2 ) {
				nbr_pval_threshold = 1e-4;
			} else {
				runtime_assert( similarity_mode == 3 );
				// This is a value of nbr_pval_threshold for which we have precomputed Ncore values
				// It also seems to give reasonable-looking clusters.
				// It corresponds to a raw TCRdist threshold of 38 with the current CDF.
				//
				nbr_pval_threshold = 0.0009;
			}
			cout <<"[ WARNING ] Setting a default value for nbr_pval_threshold: " << nbr_pval_threshold <<
				" similarity_mode: " << similarity_mode << endl;
		}

		Real const min_tcrdist_pval( 1.2e-4 );
		bool const do_clustering( cluster_arg.getValue() );
		bool const write_nbrs_to_stdout( !do_clustering );
		Size const min_cluster_size( min_cluster_size_arg.getValue() );

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

				Real co_pval(1), tcrd(0), tcrdist_pval(1), combined_pval(1);
				if ( similarity_mode ==1 ) {
					combined_pval = co_pval = compute_overlap_pvalue_with_bias( ii_occs, jj_occs, subject_bias_factors );
				} else if ( similarity_mode == 2 ) {
					co_pval = compute_overlap_pvalue_with_bias( ii_occs, jj_occs, subject_bias_factors );
					tcrd = tcrdist( ii_dtcr, all_dtcrs2[jj] );
					tcrdist_pval = tcrdist.prob_equal_or_smaller_tcrdist( tcrd );
					combined_pval = product_cdf( max( min_tcrdist_pval, tcrdist_pval ) * co_pval );
				} else if ( similarity_mode == 3 ) {
					tcrd = tcrdist( ii_dtcr, all_dtcrs2[jj] );
					combined_pval = tcrdist_pval = tcrdist.prob_equal_or_smaller_tcrdist( tcrd );
				}

				if ( combined_pval <= nbr_pval_threshold ) {
					if ( write_nbrs_to_stdout ) {
						cout << "nbrs: " << ii << ' ' << jj <<
							" tcr1: " << all_tcrs1[ii] <<
							" tcr2: " << all_tcrs2[jj] <<
							" pval: " << combined_pval;
						if ( similarity_mode == 1 ) {
							cout << " co_pval: " << co_pval << '\n';
						} else if ( similarity_mode == 2 ) {
							cout << " co_pval: " << co_pval << " tcrdist: " << tcrd << " tcrdist_pval: " << tcrdist_pval << '\n';
						} else if ( similarity_mode == 3 ) {
							cout << " tcrdist: " << tcrd << " tcrdist_pval: " << tcrdist_pval << '\n';
						}
					}
					if ( do_clustering ) {
						all_nbrs1[ii].push_back( jj );
						all_nbrs1[jj].push_back( ii );
					}
				}
			}
		}


		if ( do_clustering ) {

			Size const num_clustered_tcrs( all_tcrs1.size() );

			if ( num_clustered_tcrs < min_cluster_size ) {
				cerr << "Too few TCRs for clustering: " << num_clustered_tcrs << " < min_cluster_size= " << min_cluster_size <<
					endl;
				exit(0);
			}

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
			} else { // try to find a match in the database file
				string const filename( misc::dbdir + "DBSCAN_params_versus_num_clustered_TCRs_max_2000.txt" );
				strings lines;
				read_lines_from_file( filename, lines );
				foreach_( string & line, lines ) {
					istringstream l(line);
					string tmp;
					Size l_sim_mode, l_num_tcrs, l_min_core_nbrcount;
					Real l_nbr_pval_threshold;
					l >> tmp >> l_sim_mode >> tmp >> l_nbr_pval_threshold >> tmp >> l_num_tcrs >> tmp >> l_min_core_nbrcount;
					runtime_assert( !l.fail() );
					if ( similarity_mode == l_sim_mode && l_num_tcrs == num_clustered_tcrs &&
						fabs( log( l_nbr_pval_threshold ) - log( nbr_pval_threshold ) )<0.1 ) {
						min_core_nbrcount = l_min_core_nbrcount;
						break;
					}
				}
				if ( !min_core_nbrcount ) {
					cerr << "[ERROR] Trying to set min_core_nbrcount automagically but failed to find a match to\n" <<
						"your clustering params in the db file " << filename <<
						"\nsimilarity_mode: " << similarity_mode <<
						" nbr_pval_threshold: " << nbr_pval_threshold <<
						" num_clustered_tcrs: " << num_clustered_tcrs << '\n' <<
						"Try setting --min_core_nbrcount directly on the command line" << endl;
					exit(1);
				}
			}
			runtime_assert( min_core_nbrcount );


			if ( num_clustered_tcrs < min_core_nbrcount+1 ) {
				cerr << "Too few TCRs for clustering: " << num_clustered_tcrs <<
					" < min_core_nbrcount+1= " << min_core_nbrcount+1 << endl;
				exit(0);
			}


			/// get rid of the '1' suffix since there is only one set of tcrs/occs anyhow:
			vector< Sizes > const & all_nbrs( all_nbrs1 );
			vector< bools > const & all_occs( all_occs1 );
			strings const & all_tcrs( all_tcrs1 );

			cout << "Clustering " << num_clustered_tcrs << " TCRs " <<
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
