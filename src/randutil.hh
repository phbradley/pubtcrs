// Some utility functions involving random selections
//

#ifndef INCLUDED_randutil_HH
#define INCLUDED_randutil_HH

#include "types.hh"

#include <random>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
Real
uniform()
{
	static std::random_device rd;  //Will be used to obtain a seed for the random number engine
	static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	static std::uniform_real_distribution<> dis(0.0,1.0);
	return dis(gen);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
char
choose_random_char(
	vector< Real > const & probs,
	string const & alphabet
)
{
	runtime_assert( alphabet.size() == probs.size() );

	Real const f( uniform());

	Real total(0);

	for ( Size i=0; i<probs.size(); ++i ) {
		total += probs[i];
		if ( f <= total ) return alphabet[i];
	}
	cout << "ERROR choose_random_char f out of bounds: " << f << ' ' << total << endl;
	return alphabet[0];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
string
choose_random_string(
	map< string,Real > const & probs
)
{
	Real const f( uniform());

	Real total(0);

	for ( map<string,Real>::const_iterator it=probs.begin(), ite=probs.end(); it!=ite; ++it ) {
		total += it->second;
		if ( f <= total ) return it->first;
	}
	cout << "ERROR choose_random_string f out of bounds: " << f << ' ' << total << ' ' << probs.begin()->first << endl;
	return probs.begin()->first;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
int
choose_random_int(
	vector< Real > const & probs
)
{
	Real const f( uniform());

	Real total(0);

	for ( Size i=0; i<probs.size(); ++i ) {
		total += probs[i];
		if ( f <= total ) return i;
	}
	cout << "ERROR choose_random_int f out of bounds: " << f << ' ' << total << ' ' << probs.size() << endl;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline
void
choose_random_int_pair(
	vector< vector< Real > > const & probs,
	int & a,
	int & b
)
{
	Real const f( uniform());

	Real total(0);

	for ( Size i=0; i<probs.size(); ++i ) {
		for ( Size j=0, j_end=probs[i].size(); j<j_end; ++j ) {
			total += probs[i][j];
			if ( f<total ) {
				a = i;
				b = j;
				return;
			}
		}
	}

	cout << "ERROR choose_random_int_pair f out of bounds: " << f << ' ' << total << ' ' << probs.size() << endl;
	a=b=0;
	return;

}


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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
