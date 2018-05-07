#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "tcrdist.hh"


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

		TCRdistCalculator tcrdist; // after set_dbdir!

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

		vector< DistanceTCR_f > f1_dtcrs, f2_dtcrs;

		foreach_( string tcr, f1_tcrs ) {
			f1_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
		}
		foreach_( string tcr, f2_tcrs ) {
			f2_dtcrs.push_back( tcrdist.create_distance_tcr_f( tcr ) );
		}

		for ( Size ii=0; ii< f1_tcrs.size(); ++ii ) {
			cout << "dists: " << ii << ' '<< f1_tcrs[ii];
			DistanceTCR_f const & dtcr1( f1_dtcrs[ii] );
			for ( Size jj=0; jj< f2_tcrs.size(); ++jj ) {
				cout << ' ' << tcrdist(dtcr1, f2_dtcrs[jj]);
			}
			cout << '\n';
		}
	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
