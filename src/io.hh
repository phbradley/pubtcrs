#ifndef INCLUDED_tcrdist_cxx_io_HH
#define INCLUDED_tcrdist_cxx_io_HH

#include "misc.hh"


void
check_file(
	ifstream const & data,
	string const & filename
)
{
	if ( !data.good() ) {
		cerr << "\n\n[ERROR]  Unable to open " << filename << "\n\n" << endl;
		exit(1);
	}
}


void
read_fasta(
	string const & filename,
	map< string, string > & all_fasta
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );

	string line, id;
	while ( getline( data, line ) ) {
		strings const l( split_to_vector(line) );
		if ( line[0] == '>' ) {
			id = l[0].substr(1);
			cout << "read_fasta: " << id << ' ' << filename <<endl;
			all_fasta[id];
		} else {
			all_fasta[id] += l[0];
		}
	}
	data.close();
}

void
read_fasta_and_info(
	string const & filename,
	map< string, string > & all_fasta, // id to sequence
	map< string, string > & all_info // id to the string that follows the id on the ">" line
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );

	string line, id;
	while ( getline( data, line ) ) {
		strings const l( split_to_vector(line) );
		if ( line[0] == '>' ) {
			id = l[0].substr(1);
			all_fasta[id];
			all_info[id] = l[1];
			cout << "read_fasta: " << id << " info: " << l[1] << ' ' << filename <<endl;
		} else {
			all_fasta[id] += l[0];
		}
	}
	data.close();
}


// human_D_nucleotide_sequences.fasta
// human_J_nucleotide_sequences_with_info.fasta
// human_V_nucleotide_sequences_with_info.fasta


void
read_vj_fastas_and_info(
	map< string, string > & v_fasta,
	map< string, string > & j_fasta,
	map< string, int > & all_num_cdr3_nucleotides
)
{
	map<string,string> v_info, j_info;

	read_fasta_and_info( misc::dbdir+"human_V_nucleotide_sequences_with_info.fasta", v_fasta, v_info );
	read_fasta_and_info( misc::dbdir+"human_J_nucleotide_sequences_with_info.fasta", j_fasta, j_info );

	// parse the info
	foreach_( string v, get_keys( v_fasta ) ) {
		strings const l( split_to_vector( v_info[ v ], ":" ) );
		runtime_assert( l[0] == "num_cdr3_nucleotides" );
		runtime_assert( is_int( l[1] ) );
		runtime_assert( int_of( l[1] ) >= 0 );
		all_num_cdr3_nucleotides[v] = int_of(l[1]);
	}
	foreach_( string j, get_keys( j_fasta ) ) {
		strings const l( split_to_vector( j_info[ j ], ":" ) );
		runtime_assert( l[0] == "num_cdr3_nucleotides" );
		runtime_assert( is_int( l[1] ) );
		runtime_assert( int_of( l[1] ) >= 0 );
		all_num_cdr3_nucleotides[j] = int_of(l[1] );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool // true on success
parse_matrix_line(
	string const & line,
	string & tcr,
	bools & occs
)
{
	if ( line.substr(0,4)!="tcr:" ) return false;

	istringstream l(line);
	//tcr: V02,CAAAAGANVLTF num_subjects: 3 subjects:  53 349 456
	Size num_subjects, isub, total_subjects;
	string tag1,tag2,tag3,tag4;
	l >> tag1 >> tcr >> tag2 >> num_subjects >> tag3 >> total_subjects >> tag4;
	if ( !( tag1 == "tcr:" && tag2 == "num_subjects:" && tag3 == "of:" && tag4 == "subjects:" ) ) return false;

	occs.clear();
	occs.resize( total_subjects, false );

	for ( Size i=0; i< num_subjects; ++i ) {
		l >> isub;
		runtime_assert( isub < total_subjects ); // isub is 0-indexed !!!
		occs[isub] = true;
	}
	if ( l.fail() ) return false;
	return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
read_matrix_file(
	string const & filename,
	vector< bools > & all_occs,
	strings & all_tcrs
)
{
	all_occs.reserve( 100000 );

	cout << "start reading " << filename << endl;
	ifstream data( filename.c_str() );
	check_file( data, filename );
	string line;
	Size line_number(0);
	while ( getline(data,line) ) {
		++line_number;
		if ( line_number%1000000 == 0 ) cout << "read file " << filename << ' ' << line_number << endl;
		bools occs;
		string tcr;
		if ( parse_matrix_line( line, tcr, occs ) ) { // success
			all_tcrs.push_back( tcr );
			all_occs.resize( all_tcrs.size() );
			all_occs.back().swap( occs );
		}
	}
	data.close();

	runtime_assert( all_tcrs.size() == all_occs.size() );

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
read_features(
	string const feature_filename,
	Features & features
)
{
	ifstream data( feature_filename.c_str() );
	check_file( data, feature_filename );

	string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,8) == "feature:" ) {
			istringstream l(line);
			Feature f;
			string tag1,tag2,tag3;
			Size num_positives,num_negatives,isub;
			l >> tag1 >> f.name >> tag2 >> num_positives >> tag3;
			runtime_assert( !l.fail() && tag1 == "feature:" && tag2 == "num_positives:" && tag3 == "positives:" );
			for ( Size i=0; i<num_positives; ++i ) {
				l >> isub;
				f.poslist.push_back(isub);
			}

			l >> tag2 >> num_negatives >> tag3;
			runtime_assert( tag2 == "num_negatives:" && tag3 == "negatives:" );
			for ( Size i=0; i<num_negatives; ++i ) {
				l >> isub;
				f.neglist.push_back(isub);
			}
			features.push_back(f);
			cout << "Read new feature: " << f.name << " num_positives: " << f.poslist.size() <<
				" num_negatives: " << f.neglist.size() << endl;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// just looks for "tcr: " in line and takes the next thing...
// or if not "tcr:" in the line it takes the whole line
//
void
read_tcrs_from_file(
	string const filename,
	strings & tcrs
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );

	string line;
	while (getline( data, line ) ) {
		Size pos( line.find("tcr: ") );
		if ( pos != string::npos ) {
			istringstream l( line.substr(pos+5) );
			string tcr;
			l >> tcr;
			tcrs.push_back( tcr );
		} else {
			tcrs.push_back( line );
		}
	}
	data.close();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
void
read_floats_from_file(
	string const filename,
	Reals & vals
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );

	while ( !data.fail() ) {
		Real f;
		data >> f;
		if ( !data.fail() ) vals.push_back(f);
	}
	data.close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
read_lines_from_file(
	string const filename,
	strings & lines
)
{
	ifstream data( filename.c_str() );
	check_file( data, filename );
	string line;
	while ( getline( data, line ) ) lines.push_back( line );
	data.close();
}

#endif
