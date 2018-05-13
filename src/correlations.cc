#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "overlap_probs.hh"

#include <tclap/CmdLine.h>

int main(int argc, char** argv)
{

	try {

		TCLAP::CmdLine cmd("Find correlations between TCR occurrences provided in the matrix and feature "
			"occurrences provided in the features file", ' ', "0.1");

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"./db/","path",cmd);

 		TCLAP::ValueArg<Real> min_pval_threshold_arg("q","min_pval_threshold","P-value threshold for writing out the "
			"*best* TCR-feature correlation for each TCR",false,0,"double",cmd);

 		TCLAP::ValueArg<Real> pval_threshold_arg("p","pval_threshold","P-value threshold for writing out a "
			"TCR-feature correlation",true,0,"double",cmd);

 		TCLAP::ValueArg<string> features_arg("f","features","File containing the set of features for comparing "
			"to the TCR occurrences",true,"unk","string",cmd);

 		TCLAP::ValueArg<string> matrix_arg("m","matrix","File containing matrix of TCR occurrences, "
			"to be compared with the feature-positive subsets provided in the features file",true,"unk","string",cmd);

		cmd.parse( argc, argv );

		set_dbdir( database_arg.getValue() );

		Real const pval_threshold( pval_threshold_arg.getValue() );
		Real const min_pval_threshold( min_pval_threshold_arg.getValue() );

		string const matrix_filename( matrix_arg.getValue() );
		string const features_filename( features_arg.getValue() );

		// read a list of features
		Features features;
		read_features( features_filename, features );

		// read a matrix file
		ifstream data( matrix_filename.c_str() );
		check_file( data, matrix_filename );

		string line, tcr;
		bools tcr_occs;

		while ( getline( data, line ) ) {
			if ( parse_matrix_line( line, tcr, tcr_occs ) ) {
				string best_line;
				Real best_pval( 10 );
				foreach_( Feature const & f, features ) {

					Size fpos_with(0), fpos_wout(0), fneg_with(0), fneg_wout(0); // _with or _wout the TCR

					foreach_( Size i, f.poslist ) {
						if ( tcr_occs[i] ) ++fpos_with;
						else               ++fpos_wout;
					}

					foreach_( Size i, f.neglist ) {
						if ( tcr_occs[i] ) ++fneg_with;
						else               ++fneg_wout;
					}

					Size const fpos_total( fpos_with + fpos_wout );
					Size const total_with( fpos_with + fneg_with );
					Size const total( fpos_with + fpos_wout + fneg_with + fneg_wout );

					// what are the odds of seeing that large of an intersection, or larger??
					Real const pval
						( compute_overlap_pvalue( fpos_with, fpos_total, total_with, total ) );

					if ( pval < best_pval || pval <= pval_threshold ) {
						ostringstream out;
						out << "pval: " << pval <<
							" tcr: " << tcr <<
							" feature: " << f.name <<
							" overlap: " << fpos_with <<
							" feature_total: " << fpos_total <<
							" tcr_total: " << total_with <<
							" total: " << total << '\n';

						if ( pval <= pval_threshold ) cout << out.str();

						if ( pval < best_pval ) {
							best_pval = pval;
							best_line = "min_" + out.str();
						}

					}
				} // f in features

				if ( best_pval <= min_pval_threshold ) {
					cout << best_line;
				}
			} // parse success
		} //getline
		data.close();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
