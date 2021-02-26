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
	bool const only_nndists, // if TRUE, don't write out the distances
	strings const & f1_tcrs,
	vector< T > const & f1_dtcrs,
	strings const & f2_tcrs,
	vector< T > const & f2_dtcrs,
	TCRdistCalculator const & tcrdist
)
{


	for ( Size ii=0; ii< f1_tcrs.size(); ++ii ) {
		if ( !( only_nndists || terse() ) ) cout << "dists: " << ii << ' '<< f1_tcrs[ii];
		Reals dists;
		for ( Size jj=0; jj< f2_tcrs.size(); ++jj ) {
			Real const dist( tcrdist(f1_dtcrs[ii], f2_dtcrs[jj]) );
			if ( !only_nndists ) cout << ' ' << dist;
			if ( ii != jj || (!self_distances) ) {
				dists.push_back( dist ); // could include 0 for ii==jj if self_distances
			}
		}
		if ( !only_nndists ) cout << '\n';

		sort( dists.begin(), dists.end() ); // now in increasing order

		if ( !terse() ) {
			cout << "NNDIST tcr1_index: " << ii << " tcr1: " << f1_tcrs[ii] <<
				" mindist_to_file2_tcrs: " << F(9,3,dists.front()) <<
				" nndistance_10P_wrt_file2_tcrs: " << F(9,3,compute_NNdistance_nosorting( dists )) << endl;
		}
	}
}



int main(int argc, char** argv)
{
	try { // to catch tclap exceptions

		TCLAP::CmdLine cmd( "Calculate TCRdist distances between TCRs in one file or between two files "
			"(for example, annotation TCRs versus public TCRs). The TCRs can be defined at the V-beta family "
			"level (e.g. V19,CASSIRSSYEQYF) or the V-allele level (e.g. TRBV19*01,CASSIRSSYEQYF or "
			"TRAV3*01,CAVPPDSWGKLQF). If two files are provided, NN-distance scores for the TCRs in the "
			"first file will be computed with respect to the repertoire of TCR chains in the second file.",
			' ', "0.1" );

 		TCLAP::ValueArg<string> chain_arg("c","chain","TCR chain (A or B) in case this can't be determined "
			"from the tcr strings themselves", false, "X", "string", cmd);

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"./db/","path",cmd);

		TCLAP::SwitchArg terse_arg("t","terse", "terse output", cmd, false);

		TCLAP::SwitchArg only_nndists_arg("n","only_nndists","Don't write out the distances, "
			"just the NNDIST scores", cmd, false);

 		TCLAP::ValueArg<string> tcrs_file2_arg("j","tcrs_file2","File containing the second set of TCRs. If "
			"not provided will compute file1-vs-file1 matrix of distances.",false,"","string",cmd);

 		TCLAP::ValueArg<string> tcrs_file1_arg("i","tcrs_file1","File containing TCRs for distance computation. "
			"Will compute matrix of distances between TCRs in this file and the TCRs in tcrs_file2 if that "
			"argument is provided, otherwise just the matrix of self-distances.",true,"unk","string",cmd);

		cmd.parse( argc, argv );

		set_dbdir( database_arg.getValue() );

		set_terse( terse_arg.getValue() );

		string const tcrs_file1( tcrs_file1_arg.getValue() );
		string const tcrs_file2( tcrs_file2_arg.getValue() );
		bool const only_nndists( only_nndists_arg.getValue() );

		bool const self_distances( tcrs_file2.empty() || tcrs_file1 == tcrs_file2 );

		strings f1_tcrs, f2_tcrs;

		read_tcrs_from_file( tcrs_file1, f1_tcrs );
		if ( self_distances ) {
			f2_tcrs = f1_tcrs;
		} else {
			read_tcrs_from_file( tcrs_file2, f2_tcrs );
		}

		// which kind of input do we have? V-family or exact V-gene?
		bool by_family;
		char tcr_chain( chain_arg.getValue()[0] ); // will be 'X' if not specified

		{ // hacky stuff here
			string const random_tcr( f1_tcrs[ f1_tcrs.size()/2 ] ); // the first could be a csv header, for example
			if ( random_tcr[0] == 'V' ) {
				by_family = true;
				if ( tcr_chain=='X') tcr_chain='B';
			} else {
				by_family = false; // full allele-level gene information
				runtime_assert( random_tcr.substr(0,2) == "TR" );
				if (tcr_chain=='X') tcr_chain = random_tcr[2];
			}
		} // scope for io checking

		TCRdistCalculator tcrdist( tcr_chain ); // after calling set_dbdir, and tcr_chain determined

		for ( int i=f1_tcrs.size()-1; i>=0; --i ) {
			if ( !tcrdist.check_tcr_string_ok( f1_tcrs[i] ) ) {
				cout << "[WARNING] bad file1 tcr: " << i << ' ' << f1_tcrs[i] << endl;
				f1_tcrs.erase( f1_tcrs.begin()+i );
			}
		}
		for ( int i=f2_tcrs.size()-1; i>=0; --i ) {
			if ( !tcrdist.check_tcr_string_ok( f2_tcrs[i] ) ) {
				cout << "[WARNING] bad file2 tcr: " << i << ' ' << f2_tcrs[i] << endl;
				f2_tcrs.erase( f2_tcrs.begin()+i );
			}
		}


		if ( by_family ) {
			vector< DistanceTCR_f > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
			}
			compute_distances( self_distances, only_nndists, f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist );
		} else {

			vector< DistanceTCR_g > f1_dtcrs, f2_dtcrs;
			foreach_( string tcr, f1_tcrs ) {
				f1_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			foreach_( string tcr, f2_tcrs ) {
				f2_dtcrs.push_back( tcrdist.create_distance_tcr_g( tcr ) );
			}
			compute_distances( self_distances, only_nndists, f1_tcrs, f1_dtcrs, f2_tcrs, f2_dtcrs, tcrdist );
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
