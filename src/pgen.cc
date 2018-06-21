#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "vdj_probs.hh"

#include <tclap/CmdLine.h>

// allow either ',' or ';' as separator
strings
get_alleles_from_string( string const & tag )
{
	if ( tag.find(",") != string::npos ) return split_to_vector( tag, "," );
	if ( tag.find(";") != string::npos ) return split_to_vector( tag, ";" );
	return make_vector( tag );
}



int main(int argc, char** argv)
{

	try {

		TCLAP::CmdLine cmd("Compute VDJ rearrangement probabilities for TCR beta chains. Three options "
			"are available for the input TCRs:\n"
			" (1) the V-family,CDR3aa format, like V19,CASSIRSSYEQYF (plain text input file)\n"
			" (2) V-alleles and J-alleles (1 or more of each) and the CDR3 protein sequence (tsv formatted input file "
			" containing the fields: v_alleles j_alleles cdr3_protseq)\n"
			" (3) V-alleles and J-alleles (1 or more of each ) and the CDR3 nucleotide sequence (tsv formatted input file "
			" containing the fields: v_alleles j_alleles cdr3_nucseq)\n"
			"See examples in pubtcrs/test/pgen/run.bash",
			' ', "0.1");

		// path to database files
 		TCLAP::ValueArg<std::string> database_arg("d","database","Path to database directory",
			false, "./db/", "path", cmd );

		// option to specify an alternative model params file
 		TCLAP::ValueArg<std::string> model_file_arg("m","model_file", "File containing model params",
			false, "", "string", cmd );

		// the input file containing the tcrs
 		TCLAP::ValueArg<std::string> input_file_arg("i","input_file","File containing TCRs",
			true, "unk", "string", cmd );

		// an optional output file for writing the pgen output-- make tsv output cleaner
 		TCLAP::ValueArg<std::string> output_file_arg("o","output_file","File for output of probabilites",
			false, "", "string",cmd );

		cmd.parse( argc, argv );

		// Get the value parsed by each arg.
		set_dbdir( database_arg.getValue() );
		set_verbose( false );

		string const input_file( input_file_arg.getValue() ), output_file( output_file_arg.getValue() );

		strings lines;
		read_lines_from_file( input_file, lines ); // the tcrs

		ofstream file_out;
		if ( output_file.size() ) file_out.open( output_file.c_str() );
		ostream & out( output_file.empty() ? cout : file_out );

		JunctionCounts sampling_jc;
		if ( model_file_arg.getValue().empty() ) {
			//this is an awesome name
			string const basename("tmp.junction_stats_full_new.log.probs_adjusted_adjusted_adjusted_adjusted_adjusted");
			sampling_jc.read_probs_and_counts_from_file( misc::dbdir + basename );
		} else {
			sampling_jc.read_probs_and_counts_from_file( model_file_arg.getValue() );
		}

		// what kind of input did we get?? silly check: no tabs in first line
		bool const vfamily_io( split_to_vector( lines.front(), "\t" ).size() == 1 );

		if ( vfamily_io ) {
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
				out << "pgen: " << protprob << " tcr: " << tcr << endl;
			}
		} else { // .tsv file input
			map<string,Size> tag2index;
			{ // parse the header
				strings const l( split_to_vector( lines.front(), "\t" ) );
				for ( Size i=0; i<l.size(); ++i ) tag2index[ l[i] ] = i;
			}
			foreach_( string tag, make_vector( string("v_alleles"), string( "j_alleles" ) ) ) { // required tags in header
				if ( tag2index.count(tag) == 0 ) {
					cerr << "[ERROR] Missing required tag: " << tag << " in input .tsv file" << endl;
					exit(1);
				}
			}
			// output a new header for the file
			out << lines.front() << "\tPgen" << endl;
			bool const nucseq_io( tag2index.count("cdr3_nucseq") );
			if ( !nucseq_io && tag2index.count("cdr3_protseq") == 0 ) {
				cerr << "[ERROR] Need to provide either cdr3_nucseq or cdr3_protseq for Pgen calculation" << endl;
				exit(1);
			}
			for ( Size li=1; li<lines.size(); ++li ) {
				strings const l( split_to_vector( lines[li], "\t" ) );
				// allows ";" or "," as separator between multiple alleles
				strings const v_alleles( get_alleles_from_string( l[ tag2index["v_alleles"] ] ) );
				strings const j_alleles( get_alleles_from_string( l[ tag2index["j_alleles"] ] ) );
				if ( nucseq_io ) {
					string const cdr3_nucseq( l[ tag2index["cdr3_nucseq"] ] );
					Real const nucprob( sampling_jc.calc_cdr3_nucseq_probability( cdr3_nucseq, v_alleles, j_alleles ) );
					out << lines[li] << '\t' << nucprob << endl;
				} else {
					string const cdr3_protseq( l[ tag2index["cdr3_protseq"] ] );
					Real const protprob( sampling_jc.calc_cdr3_protseq_probability( cdr3_protseq, v_alleles, j_alleles ) );
					out << lines[li] << '\t' << protprob << endl;
				}
			}
		}
		if ( output_file.size() ) file_out.close();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
