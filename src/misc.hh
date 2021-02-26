// Miscellaneous generic helper functions
//

#ifndef INCLUDED_misc_HH
#define INCLUDED_misc_HH

#include "types.hh"

#include <random>
#include <iomanip>

#include <tclap/CmdLine.h>

// mutable global state is bad

namespace misc {
string dbdir;
bool verbose_(false);
bool terse_(false);
}

inline
bool
verbose()
{
	return misc::verbose_;
}

void
set_verbose(
	bool const setting
)
{
	misc::verbose_ = setting;
}

inline
bool
terse()
{
	return misc::terse_;
}

void
set_terse(
	bool const setting
)
{
	misc::terse_ = setting;
}

void
set_dbdir(
	string const & setting
)
{
	misc::dbdir = setting;
	if ( misc::dbdir.back() != '/' ) misc::dbdir += "/";
}



void
my_exit(
	string const & file,
	int const line,
	string const & message
)
{

	ostringstream oss;
	if ( ! message.empty() ) oss << "\n" << "ERROR: " << message << "\n";
	oss << "ERROR:: Exit from: " << file << " line: " << line << "\n";
	string failure_message = oss.str();
	cerr << failure_message << flush;
	assert( false ); // for gdb?
	exit(1);
}

#define runtime_assert(_Expression) if ( !(_Expression) ) my_exit(__FILE__, __LINE__, #_Expression)

/// silly helper
void
get_mean_sdev( Reals const & vals, Real & mean, Real & sdev )
{
	mean = sdev = 0.0;
	if ( vals.empty() ) return;

	foreach_( Real val, vals ) mean += val;
	mean /= vals.size();
	foreach_( Real val, vals ) sdev += ( val - mean ) * ( val - mean );
	sdev = std::sqrt( sdev / vals.size() );
}

/// @details  Split a string to a vector
vector< string >
split_to_vector( string const & s )
{
	vector< string > v;

	istringstream l( s );
	string tag;
	l >> tag;
	while ( !l.fail() ) {
		v.push_back( tag );
		l >> tag;
	}
	return v;
}

/// @details  Split a string to a vector using sep as a separator
vector< string >
split_to_vector(
	string s,
	string const & sep
)
{
	vector< string > v;

	Size pos( s.find( sep ) );
	while ( pos != string::npos ) {
		v.push_back( s.substr( 0, pos ) );
		s.erase( 0, pos + sep.size() );
		pos = s.find( sep );
	}
	assert( s.find( sep ) == string::npos );
	v.push_back( s );
	return v;
}

// more hacky stuff
inline
bool
is_whitespace( string const & s )
{
	if ( s.empty() ) {
		return true;
	} else {
		return ( s.find_last_not_of( " \t\000" ) == string::npos );
	}
}

template< typename T >
inline
bool
is_type( string const & s )
{
	if ( is_whitespace( s ) ) {
		return false;
	} else {
		istringstream t_stream(s);
		T t;
		t_stream >> t;
		return ( ( t_stream ) && ( t_stream.eof() ) );
	}
}

inline
bool
is_int( string const & s )
{
	return is_type< int >( s );
}

inline
int
int_of( string const & s ) {
	istringstream l(s);
	int val;
	l >> val;
	runtime_assert( !l.fail() );
	return val;
}

inline
Real // actually a double
float_of( string const & s ) {
	istringstream l(s);
	Real val;
	l >> val;
	runtime_assert( !l.fail() );
	return val;
}

// template< typename T >
// inline
// string
// lead_zero_string_of(
// 	T const & t,
// 	int const w // Minimum width
// )
// {
// 	ostringstream t_stream;
// 	t_stream << internal << uppercase
// 	 << setw( w ) << setfill( '0' ) << setprecision( TypeTraits< T >::precision() ) << t;
// 	return t_stream.str();
// }

// annoying but don't want to rely on c++11 to_string function...
// template< typename T >
// inline
// string
// my_string_of(
// 	T const & t
// )
// {
// 	ostringstream os;
// 	os << t;
// 	return os.str();
// }

//
// eg, "V04"
//
string
get_v_family_from_v_gene( string const & g )
{
	runtime_assert( g.substr(0,2) == "TR" && g[3] == 'V' );
	// runtime_assert( g.substr(0,4) == "TRBV" );
	Size numlen(0);
	while ( is_int( g.substr(4,numlen+1) ) )++numlen;
	if (!numlen) { // for example, TRGVA (which is a pseudogene)
		return "V00";
	}
	runtime_assert( numlen );
	Size const vno( int_of( g.substr(4,numlen) ) );
	string const zeropad( vno < 10 ? "0" : "" );
	string const v_family( "V" + zeropad + to_string( vno ) );
	return v_family;
}


map<string,strings>
setup_v_family2v_genes( strings const & v_genes )
{
	map<string,strings> v_family2v_genes;
	foreach_( string g, v_genes ) {
		v_family2v_genes[ get_v_family_from_v_gene( g ) ].push_back( g );
	}
	return v_family2v_genes;
}

template < typename T1, typename T2 >
vector< T1 >
get_keys( map< T1, T2 > const & m )
{
	vector< T1 > ks;
	for ( typename map< T1,T2 >::const_iterator it= m.begin(); it != m.end(); ++it ) ks.push_back( it->first );
	return ks;
}


template < class T >
bool
has_element( vector< T > const & v, T const & t )
{
	return ( find( v.begin(), v.end(), t ) != v.end() );
}

/// never can remember the order...
template < class T >
bool
has_element( T const & t, vector< T > const & v )
{
	return ( find( v.begin(), v.end(), t ) != v.end() );
}


template < class T >
Size
vector_index( T const & t, vector< T > const & v )
{
	return find( v.begin(), v.end(), t ) - v.end();
}


template < class T >
Size
vector_index( vector< T > const & v, T const & t )
{
	return find( v.begin(), v.end(), t ) - v.end();
}


template < class T >
vector< T >
make_vector( T const & t )
{
	vector<T> v;
	v.push_back(t);
	return v;
}

template < class T >
vector< T >
make_vector( T const & t1,  T const & t2 )
{
	vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	return v;
}

template < class T >
vector< T >
make_vector( T const & t1,  T const & t2 ,  T const & t3 )
{
	vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	return v;
}



/// these are inspired by the output functions from Rosetta
inline
string
F( int const width, int const decimals, float const & t )
{
	stringstream out;
	out << fixed << showpoint << setprecision( decimals ) << setw( width ) << t;
	return out.str();
}

inline
string
F( int const width, int const decimals, double const & t )
{
	stringstream out;
	out << fixed << showpoint << setprecision( decimals ) << setw( width ) << t;
	return out.str();
}

template< typename T >
inline
string
I( int const width, T const & value )
{
	ostringstream out;
	out << right << setw( width ) << value;
	return out.str();
}


template < class T >
T
max( vector<T> const & v )
{
	return *max_element( v.begin(), v.end() );
}

template < class T >
T
min( vector<T> const & v )
{
	return *min_element( v.begin(), v.end() );
}

#endif
