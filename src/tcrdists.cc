#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Real
compute_NNdistance_nosorting(
	Size const num,
	Reals const & dists // should already be sorted in increasing order
)
{
	runtime_assert( dists.size() >= num );

	Real nbrdist(0), totalwt(0);
	for ( Size i=0; i<num; ++i ) {
		Real const wt( 1.0 - Real(i)/num );
		totalwt += wt;
		nbrdist += wt * dists[i];
	}
	return nbrdist / totalwt;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Real
compute_NNdistance_nosorting(
	Reals const & dists // should already be sorted in increasing order
)
{
	return compute_NNdistance_nosorting( max( Size(1), dists.size()/10 ), dists );
}



// this is templated because there are different flavors of TCR depending on whether we know the V gene at the family
// or allele level
template< typename T >
void
compute_distances(
	bool const self_distances,
	strings const & f1_tcrs,
	vector< T > const & f1_dtcrs,
	strings const & f2_tcrs,
	vector< T > const & f2_dtcrs,
	TCRdistCalculator const & tcrdist
)
{


	for ( Size ii=0; ii< f1_tcrs.size(); ++ii ) {
		cout << "dists: " << ii << ' '<< f1_tcrs[ii];
		Reals dists;
		for ( Size jj=0; jj< f2_tcrs.size(); ++jj ) {
			Real const dist( tcrdist(f1_dtcrs[ii], f2_dtcrs[jj]) );
			cout << ' ' << dist;
			if ( !self_distances || ii!=jj ) dists.push_back( dist );
		}
		cout << '\n';
		if ( !self_distances ) runtime_assert( dists.size() == f2_tcrs.size() );
		sort( dists.begin(), dists.end() ); // now in increasing order

		cout << "NNDIST tcr1_index: " << ii << " tcr1: " << f1_tcrs[ii] <<
			" nndistance_10P_wrt_file2_tcrs: " << F(9,3,compute_NNdistance_nosorting( dists )) << endl;
	}

}



int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "Calculate TCRdist distances between TCRs in one file or between two files "
			"(for example, annotation TCRs versus public TCRs).", ' ', "0.1" );

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"./db/","path",cmd);

 		TCLAP::ValueArg<string> tcrs_file2_arg("j","tcrs_file2","File containing the second set of TCRs. If "
			"not provided will compute file1-vs-file1 matrix of distances.",false,"","string",cmd);

 		TCLAP::ValueArg<string> tcrs_file1_arg("i","tcrs_file1","File containing TCRs for distance computation. "
			"Will compute matrix of distances between TCRs in this file and the TCRs in tcrs_file2 if that "
			"argument is provided, otherwise just the matrix of self-distances.",true,"unk","string",cmd);

		cmd.parse( argc, argv );

		set_dbdir( database_arg.getValue() );


		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );

		bool const self_distances( tcrs_file2.empty() );


		strings f1_tcrs, f2_tcrs;

		read_tcrs_from_file( tcrs_file1, f1_tcrs );
		if ( self_distances ) {
			f2_tcrs = f1_tcrs;
		} else {
			read_tcrs_from_file( tcrs_file2, f2_tcrs );
		}

		// which kind of input do we have? V-family or exact V-gene?
		bool by_family;
		char tcr_chain( 'B' ); // the default

		{ // hacky stuff here
			string const first_tcr( f1_tcrs.front() );
			if ( first_tcr[0] == 'V' ) {
				by_family = true;
				// sanity checking:
				foreach_( string tcr, f1_tcrs ) {
					runtime_assert( tcr[0] == 'V' && tcr[3] == ',' && split_to_vector(tcr,",").size() == 2 );
				}
				foreach_( string tcr, f2_tcrs ) {
					runtime_assert( tcr[0] == 'V' && tcr[3] == ',' && split_to_vector(tcr,",").size() == 2 );
				}
			} else {
				by_family = false; // full allele-level gene information
				runtime_assert( first_tcr.substr(0,2) == "TR" );
				tcr_chain = first_tcr[2];
				// sanity checking:
				runtime_assert( string("AB").find(tcr_chain) != string::npos );
				foreach_( string tcr, f1_tcrs ) {
					runtime_assert( tcr.substr(0,2) == "TR" && tcr[2] == tcr_chain && split_to_vector(tcr,",").size() == 2 );
				}
				foreach_( string tcr, f2_tcrs ) {
					runtime_assert( tcr.substr(0,2) == "TR" && tcr[2] == tcr_chain && split_to_vector(tcr,",").size() == 2 );
				}
			}
		} // scope for io checking

		TCRdistCalculator tcrdist( tcr_chain ); // after calling set_dbdir, and tcr_chain determined

		if ( by_family ) {
			vector< DistanceTCR_f > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			compute_distances( self_distances, f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist );
		} else {

			vector< DistanceTCR_g > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			compute_distances( self_distances, f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist );
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
