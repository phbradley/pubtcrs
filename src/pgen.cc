#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "vdj_probs.hh"

#include <tclap/CmdLine.h>



int main(int argc, char** argv)
{

	try {

		TCLAP::CmdLine cmd("Compute VDJ rearrangement probabilities for TCR beta chain AA seqs "
			"using the V-family,CDR3seq format, like V19,CASSIRSSYEQYF", ' ', "0.1");

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",false,
			"/home/pbradley/tcr_scripts/cxx/db/","string",cmd);

		// option to specify an alternative model params file
 		TCLAP::ValueArg<std::string> model_file_arg("m","model_file", "File containing model params",
			false, "", "string", cmd);

		// the input file containing the
 		TCLAP::ValueArg<std::string> input_file_arg("i","input_file","File containing TCRs",true,"unk","string",cmd);

		cmd.parse( argc, argv );

		// Get the value parsed by each arg.
		set_dbdir( database_arg.getValue() );

		string const input_file( input_file_arg.getValue() );

		strings lines;
		read_lines_from_file( input_file, lines ); // the tcrs

		JunctionCounts sampling_jc;
		if ( model_file_arg.getValue().empty() ) {
			//this is an awesome name
			string const basename("tmp.junction_stats_full_new.log.probs_adjusted_adjusted_adjusted_adjusted_adjusted");
			sampling_jc.read_probs_and_counts_from_file( misc::dbdir + basename );
		} else {
			sampling_jc.read_probs_and_counts_from_file( model_file_arg.getValue() );
		}

		map<string,strings> const v_family2v_genes( setup_v_family2v_genes( sampling_jc.v_genes() ) );

		foreach_( string tcr, lines ) {
			strings const l( split_to_vector( tcr, "," ) );
			runtime_assert( l.size() == 2 );
			runtime_assert( l[0].size() == 3 && l[0][0] == 'V' );

			string const vfam( l[0] ), cdr3_protseq( l[1] );
			if ( v_family2v_genes.count(vfam) == 0 ) {
				cout << "bad v_family " << tcr << endl;
				continue;
			}

			strings const & v_ties( v_family2v_genes.find(vfam)->second );
			strings const & j_ties( sampling_jc.j_genes() );

			Real const protprob( sampling_jc.calc_cdr3_protseq_probability( cdr3_protseq, v_ties, j_ties ) );
			cout << "pgen: " << protprob << " tcr: " << tcr << endl;
		}
	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
