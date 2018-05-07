#ifndef INCLUDED_tcrdist_cxx_tcrdist_HH
#define INCLUDED_tcrdist_cxx_tcrdist_HH

#include "misc.hh"

struct DistanceTCR_g {
	Size v_num;
	string cdr3; // from C to 'F'
};


struct DistanceTCR_f {
	Size vfam_num;
	string cdr3; // from C to 'F'
};

struct DistanceTCR_gs {
	Sizes v_nums;
	string cdr3; // from C to 'F'
};


class TCRdistCalculator {
public:

	// ctor
	TCRdistCalculator( char const ab = 'B' );

	// do error checking inside here:
	DistanceTCR_g
	create_distance_tcr_g( string const & v_gene, string const & cdr3 ) const;

	DistanceTCR_gs
	create_distance_tcr_gs( strings const & v_genes, string const & cdr3 ) const;

	DistanceTCR_f
	create_distance_tcr_f( string const & tcr ) const; // tcr looks like "V07,CASSIRSSAYEQFF"

	DistanceTCR_f
	create_distance_tcr_f( string const & v_family, string const & cdr3 ) const;

	Real
	operator()( DistanceTCR_f const & t1, DistanceTCR_f const & t2 ) const;

	Real
	operator()( DistanceTCR_g const & t1, DistanceTCR_g const & t2 ) const;

	Real
	operator()( DistanceTCR_gs const & t1, DistanceTCR_gs const & t2 ) const;

	inline
	Real
	aa_distance( char const a, char const b ) const;

	inline
	Real
	cdr3_distance( string const & a, string const & b ) const;

	string const &
	v_gene( Size const v_num ) const {
		runtime_assert( v_num < v_genes_.size() );
		return v_genes_[ v_num ];
	}

	string const &
	v_family( Size const vfam_num ) const {
		runtime_assert( vfam_num < v_families_.size() );
		return v_families_[ vfam_num ];
	}

	inline
	Real
	prob_equal_or_smaller_tcrdist( Real const dist ) const;


private:
	char const ab_;
	vector< string > v_genes_;
	vector< string > v_families_;

	map<string, Size> v_gene2v_num_;
	map<string, Size> v_family2vfam_num_;

	vector< vector< Real > > V_dist_matrix_;

	vector< vector< Real > > Vfam_dist_matrix_;

	vector< vector< Real > > AA_dist_matrix_;

	vector< Real > cdf_;

	// the V region weights and distances are rolled into the v distances
	Real weight_cdr3_region_;
	Real gap_penalty_cdr3_region_;

	// helpful
	string amino_acids_;

};


TCRdistCalculator::TCRdistCalculator(
	char const ab // = 'B'
):
	ab_( ab )
{
	runtime_assert( ab=='A' || ab=='B');
	weight_cdr3_region_ = 3;
	gap_penalty_cdr3_region_ = 12; // make optional

	Real const big_dist( 1e6 );

	amino_acids_ = "ACDEFGHIKLMNPQRSTVWY";
	runtime_assert( amino_acids_.size() == 20 );

	// load the data
	string const filename( "/home/pbradley/tcr_scripts/tcrdist_info_both_chains.txt" );
	ifstream data( filename.c_str() );

	// read the aa distance matrix
	AA_dist_matrix_.resize(26);
	// 20 distance lines
	string line;
	for ( Size i=0; i<20; ++i ) {
		getline( data, line );
		istringstream l( line );
		string tag, aastring;
		l >> tag >> aastring;
		runtime_assert( tag == "AAdist");
		runtime_assert( aastring.size() == 1 );
		runtime_assert( amino_acids_.find( aastring ) != string::npos );
		char const aa( aastring[0] );
		runtime_assert( 'A' <= aa && aa <= 'Z' );
		int const aa_ind( aa-'A' );
		runtime_assert( aa_ind >= 0 && aa_ind<26 );
		AA_dist_matrix_[ aa_ind ].resize( 26, big_dist );
		for ( Size j=0; j<20; ++j ) {
			char const bb( amino_acids_[j] );
			runtime_assert( 'A' <= bb && bb <= 'Z' );
			int const bb_ind( bb-'A' );
			l >> AA_dist_matrix_[ aa_ind ][ bb_ind ];
			// cout << "AAdist: " << char('A'+aa_ind) << ' ' << char('A'+bb_ind) << ' ' <<
			// 	AA_dist_matrix_[ aa_ind ][ bb_ind ] << endl;
		}
		runtime_assert( !l.fail() );
	}

	string const abstring( 1,ab);
	string const numtag( "num_V" + abstring + "_genes" ), disttag( "V"+abstring+"dist" );
	getline(data,line);
	while ( line.substr(0,numtag.size()) != numtag ) getline(data,line);

	Size num_v_genes(0);
	{ // parse num_v_genes
		istringstream l(line);
		string tag;
		l>> tag >> num_v_genes;
		runtime_assert( tag == numtag && !l.fail() );
		// cout << line << endl;
	}

	V_dist_matrix_.resize( num_v_genes );

	Size vnum(0);
	v_genes_.clear(); // a vector<string>
	while ( getline( data, line ) ) {
		istringstream l( line );
		string tag,g;
		l >> tag >> g;
		if ( tag != disttag ) break;
		v_gene2v_num_[ g ] = vnum;
		v_genes_.push_back( g );
		runtime_assert( vnum < num_v_genes );
		V_dist_matrix_[ vnum ].resize( num_v_genes, big_dist );
		for ( Size i=0; i<num_v_genes; ++i ) {
			l >> V_dist_matrix_[ vnum ][ i ];
		}
		runtime_assert( !l.fail() );
		++vnum;
	}
	data.close();

	// could now get a v-family distance matrix?
	map<string,strings> const v_family2v_genes( setup_v_family2v_genes( get_keys( v_gene2v_num_ ) ) );

	v_families_ = get_keys( v_family2v_genes );
	std::sort( v_families_.begin(), v_families_.end() );

	Size const num_v_families( v_families_.size() );

	// cout << "num_v_families: " << num_v_families << endl;

	Vfam_dist_matrix_.resize( num_v_families );
	for ( Size i=0; i< num_v_families; ++i ) {
		string const ifam( v_families_[i] ); // note i+1
		v_family2vfam_num_[ ifam ] = i;
		Vfam_dist_matrix_[i].resize( num_v_families, big_dist );
		for ( Size j=0; j< num_v_families; ++j ) {
			string const jfam( v_families_[j] );
			Real mindis( big_dist ), maxdis(0);
			foreach_( string ig, v_family2v_genes.find( ifam )->second ) {
				foreach_( string jg, v_family2v_genes.find( jfam )->second ) {
					mindis = min( mindis, V_dist_matrix_[ v_gene2v_num_[ ig ] ][ v_gene2v_num_[ jg ] ] );
					maxdis = max( maxdis, V_dist_matrix_[ v_gene2v_num_[ ig ] ][ v_gene2v_num_[ jg ] ] );
				}
			}
			// cout << "Vfam_dist_matrix: " << ifam << ' ' << jfam << ' ' << mindis << " maxdis: " << maxdis <<
			// 	endl;
			Vfam_dist_matrix_[ i ][ j ] = mindis;
		}
	}


	{ // this part only works for beta right now!
		string filename;
		if ( fabs( gap_penalty_cdr3_region_ - 8. ) < 1e-3 ) {
			filename = "/home/pbradley/tcr_scripts/tcrdist_cdf_randpubtcrs_v1.txt" ;
		} else {
			//runtime_assert( fabs( gap_penalty_cdr3_region_-12 ) < 1e-3 );
			filename = "/home/pbradley/tcr_scripts/tcrdist_cdf_randpubtcrs_v1_cdr3_gap_penalty_12_ge_5_subjects.txt";
			//filename = "/home/pbradley/tcr_scripts/tcrdist_cdf_randpubtcrs_v1_cdr3_gap_penalty_12.txt";
		}
		ifstream data( filename.c_str() );
		string line, tmp;
		vector<Size> totals;
		Size newtotal(0), idist, count, total;
		while ( getline( data, line ) ) {
			istringstream l(line);
			l >> tmp >> idist >> tmp >> count >> tmp >> total;
			newtotal += count;
			runtime_assert( newtotal == total );
			runtime_assert( idist == totals.size() );
			totals.push_back( total );
		}
		data.close();
		Size maxdist( totals.size()*2 );
		cdf_.clear();
		cdf_.resize( maxdist+1, 1.0 );
		for ( Size i=0; i<totals.size(); ++i ) {
			cdf_[i] = Real(totals[i]) / newtotal;
			// cout << "tcrdist_cdf: " << i << ' ' << cdf_[i] << endl;
		}
	}


}

inline
Real
TCRdistCalculator::prob_equal_or_smaller_tcrdist( Real const dist ) const
{ // dist should basically be an integer
	runtime_assert( ab_ == 'B' );
	return cdf_[ int( dist+0.5 ) ];
}

inline
Real
TCRdistCalculator::aa_distance( char const a, char const b ) const
{
	// could do some sanity checking here...
	return AA_dist_matrix_[a-'A'][b-'A'];
}

inline
Real
TCRdistCalculator::cdr3_distance( string const & a, string const & b ) const
{
	static int const ntrim(3), ctrim(2); // params

	int const alen( a.size() ), blen( b.size() );

	string const & shortseq( alen <= blen ? a : b );
	string const &  longseq( alen <= blen ? b : a );

	int const lenshort( min(alen,blen) ), lenlong( max(alen,blen) ), lendiff( lenlong-lenshort ),
		gappos( min( 6, 3 + (lenshort-5)/2 ) ), remainder( lenshort-gappos );

	runtime_assert( lenshort > 5 );

	Real dist( 0 );

	for ( int i=ntrim; i<gappos; ++i ) {
		dist += aa_distance( shortseq[i], longseq[i] );
	}
	for ( int i=ctrim; i<remainder; ++i ) {
		dist += aa_distance( shortseq[ lenshort-1-i ], longseq[ lenlong-1-i ] );
	}

	return weight_cdr3_region_ * dist + lendiff * gap_penalty_cdr3_region_; // gap penalty is not also weighted
}

Real
TCRdistCalculator::operator()(
	DistanceTCR_g const & t1,
	DistanceTCR_g const & t2
) const
{
	return V_dist_matrix_[ t1.v_num ][ t2.v_num ] + cdr3_distance( t1.cdr3, t2.cdr3 );
}

Real
TCRdistCalculator::operator()(
	DistanceTCR_gs const & t1,
	DistanceTCR_gs const & t2
) const
{
	Real min_vdist(1e6);
	foreach_( Size v1, t1.v_nums ) {
		foreach_( Size v2, t2.v_nums ) {
			min_vdist = min( min_vdist, V_dist_matrix_[ v1 ][ v2 ] );
		}
	}
	return min_vdist + cdr3_distance( t1.cdr3, t2.cdr3 );
}

Real
TCRdistCalculator::operator()(
	DistanceTCR_f const & t1,
	DistanceTCR_f const & t2
) const
{
	return Vfam_dist_matrix_[ t1.vfam_num ][ t2.vfam_num ] + cdr3_distance( t1.cdr3, t2.cdr3 );
}

//
DistanceTCR_f
TCRdistCalculator::create_distance_tcr_f(
	string const & vfam,
	string const & cdr3
) const
{
	// sanity check
	foreach_( char const aa, cdr3 ) {
		runtime_assert( amino_acids_.find(aa) != string::npos );
	}
	runtime_assert( cdr3.size() > 5 );
	runtime_assert( v_family2vfam_num_.count( vfam ) );

	DistanceTCR_f dtcr;
	dtcr.vfam_num = v_family2vfam_num_.find( vfam )->second;
	dtcr.cdr3 = cdr3;

	return dtcr;
}

DistanceTCR_f
TCRdistCalculator::create_distance_tcr_f(
	string const & tcr
) const
{
	strings const l( split_to_vector( tcr, "," ) );
	runtime_assert( l.size() == 2 );
	runtime_assert( l[0].size() == 3 && l[0][0] == 'V' );

	return create_distance_tcr_f( l[0], l[1] );
}

DistanceTCR_g
TCRdistCalculator::create_distance_tcr_g(
	string const & vgene,
	string const & cdr3
) const
{
	// sanity check
	foreach_( char const aa, cdr3 ) {
		runtime_assert( amino_acids_.find(aa) != string::npos );
	}
	runtime_assert( cdr3.size() > 5 );
	runtime_assert( v_gene2v_num_.count( vgene ) );

	DistanceTCR_g dtcr;
	dtcr.v_num = v_gene2v_num_.find( vgene )->second;
	dtcr.cdr3 = cdr3;

	return dtcr;
}

DistanceTCR_gs
TCRdistCalculator::create_distance_tcr_gs(
	strings const & vgenes,
	string const & cdr3
) const
{
	// sanity check
	foreach_( char const aa, cdr3 ) {
		runtime_assert( amino_acids_.find(aa) != string::npos );
	}
	runtime_assert( cdr3.size() > 5 );

	DistanceTCR_gs dtcr;
	foreach_( string const & vgene, vgenes ) {
		runtime_assert( v_gene2v_num_.count( vgene ) );
		dtcr.v_nums.push_back( v_gene2v_num_.find( vgene )->second );
	}
	dtcr.cdr3 = cdr3;

	return dtcr;
}



#endif
