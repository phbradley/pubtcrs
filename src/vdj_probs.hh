// Routines for computing rearrangement probabilities (P_gen)
//


#ifndef INCLUDED_vdj_probs_HH
#define INCLUDED_vdj_probs_HH


#include "types.hh"
#include "misc.hh"
#include "io.hh"
#include "sequtil.hh"
#include "randutil.hh"

// Namespace for hard-coded parameters. Should make these more configurable at some point.
//
// Some of these are a little restrictive: they are probably fine for looking at public TCRs, which are generally
// pretty close to germline. But private TCRs may have more insertions or trims than we allow here.
// It wouldn't be hard to recollect new stats... email pbradley@fredhutch.org if interested.
//
namespace vdj_probs {
Size const max_trim( 15 ); // does not include pnucs
Size const max_vj_insert( 25 );
Size const max_vdj_insert( 15 ); // max insertion on either side of d gene
Size const min_d_trimmed_length( 3 ); // we do also allow for not finding a d gene at all, since they can be trimmed away
Size const max_p_nucleotides( 2 );
Size const max_extra_vj_trim_for_probs(3);
Size const max_extra_d_trim_for_probs(3);
Size const max_d_gap( 2 );
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple struct to hold information on a rearrangement parse. Used when collecting stats.
//
struct JunctionParse {
	 // d_gene is empty string if trimmed too far back to find
	string v_gene, d_gene, j_gene,
		cdr3_nucseq, cdr3_newseq,
		vd_insertseq, dj_insertseq, vj_insertseq,
		v_ties, d_ties, j_ties,
		cdr3_nucseq_from_v, cdr3_nucseq_from_d, cdr3_nucseq_from_j;

	int v_trim, d0_trim, d1_trim, j_trim, vd_insert, dj_insert, vj_insert;

	void
	clear(); // just clear the strings, don't worry about the ints

	void
	clear_dashes();

};

inline
string
dash_if_empty( string const & s ) {
	if ( s.empty() ) return "-";
	else return s;
}

ostream &
operator<< ( ostream & out, JunctionParse const & jp )
{
	runtime_assert( jp.cdr3_nucseq.size() == jp.cdr3_newseq.size() );
	string cdr3_mmseq, a( boost::to_upper_copy( jp.cdr3_nucseq ) ), b( boost::to_upper_copy( jp.cdr3_newseq ) );
	for ( Size i=0; i<a.size(); ++i ) {
		cdr3_mmseq.push_back( ( a[i]==b[i] ? '0' : '1' ) );
	}
	out <<
		" cdr3_nucseq: " << dash_if_empty( jp.cdr3_nucseq ) <<
		" cdr3_newseq: " << dash_if_empty( jp.cdr3_newseq ) <<
		" cdr3_mmseq: " << dash_if_empty( cdr3_mmseq ) <<
		" cdr3_mismatches: " << std::count( cdr3_mmseq.begin(), cdr3_mmseq.end(), '1' ) <<
		" in_frame: " << ( jp.cdr3_nucseq.size()%3==0 ) <<
		" v_gene: " << dash_if_empty( jp.v_gene ) <<
		" d_gene: " << dash_if_empty( jp.d_gene ) <<
		" j_gene: " << dash_if_empty( jp.j_gene ) <<
		" v_trim: "  <<  jp.v_trim <<
		" d0_trim: " << jp.d0_trim <<
		" d1_trim: " << jp.d1_trim <<
		" j_trim: "  <<  jp.j_trim <<
		" vd_insert: " << jp.vd_insert <<
		" dj_insert: " << jp.dj_insert <<
		" vj_insert: " << jp.vj_insert <<
		" vd_insertseq: " << dash_if_empty( jp.vd_insertseq ) <<
		" dj_insertseq: " << dash_if_empty( jp.dj_insertseq ) <<
		" vj_insertseq: " << dash_if_empty( jp.vj_insertseq ) <<
		" v_ties: " << dash_if_empty( jp.v_ties ) <<
		" d_ties: " << dash_if_empty( jp.d_ties ) <<
		" j_ties: " << dash_if_empty( jp.j_ties ) <<
		" cdr3_nucseq_from_v: " << dash_if_empty( jp.cdr3_nucseq_from_v ) <<
		" cdr3_nucseq_from_d: " << dash_if_empty( jp.cdr3_nucseq_from_d ) <<
		" cdr3_nucseq_from_j: " << dash_if_empty( jp.cdr3_nucseq_from_j );

	return out;
}

istream &
operator>> ( istream & data, JunctionParse & jp )
{
	string tmp;
	data >> tmp;
	if ( tmp != "cdr3_nucseq:" ) {
		data.setstate( std::ios_base::failbit );
		return data;
	}
	string cdr3_mmseq;
	Size cdr3_mismatches;
	bool in_frame;
	data >> jp.cdr3_nucseq >>
		tmp >> jp.cdr3_newseq >>
		tmp >> cdr3_mmseq >>
		tmp >> cdr3_mismatches >>
		tmp >> in_frame >>
		tmp >> jp.v_gene >>
		tmp >> jp.d_gene >>
		tmp >> jp.j_gene >>
		tmp >> jp.v_trim >>
		tmp >> jp.d0_trim >>
		tmp >> jp.d1_trim >>
		tmp >> jp.j_trim >>
		tmp >> jp.vd_insert >>
		tmp >> jp.dj_insert >>
		tmp >> jp.vj_insert >>
		tmp >> jp.vd_insertseq >>
		tmp >> jp.dj_insertseq >>
		tmp >> jp.vj_insertseq  >>
		tmp >> jp.v_ties >>
		tmp >> jp.d_ties >>
		tmp >> jp.j_ties >>
		tmp >> jp.cdr3_nucseq_from_v >>
		tmp >> jp.cdr3_nucseq_from_d >>
		tmp >> jp.cdr3_nucseq_from_j;

	if ( data.fail() ) jp.clear();
	else jp.clear_dashes(); // so dashes only live for file I/O

	return data;
}



void
JunctionParse::clear()
{
	v_gene.clear();
	d_gene.clear();
	j_gene.clear();
	cdr3_nucseq.clear();
	cdr3_newseq.clear();
	vd_insertseq.clear();
	dj_insertseq.clear();
	vj_insertseq.clear();
	cdr3_nucseq_from_v.clear();
	cdr3_nucseq_from_d.clear();
	cdr3_nucseq_from_j.clear();
	v_ties.clear();
	d_ties.clear();
	j_ties.clear();
	v_trim = d0_trim = d1_trim = j_trim = vd_insert = dj_insert = vj_insert = 0;
}

inline
void
clear_if_dash( string & s )
{
	if ( s== "-" ) s.clear();
}

void
JunctionParse::clear_dashes()
{
	clear_if_dash( v_gene );
	clear_if_dash( j_gene );
	clear_if_dash( d_gene );
	clear_if_dash( cdr3_nucseq );
	clear_if_dash( cdr3_newseq );
	clear_if_dash( vd_insertseq );
	clear_if_dash( dj_insertseq );
	clear_if_dash( vj_insertseq );
	clear_if_dash( cdr3_nucseq_from_v );
	clear_if_dash( cdr3_nucseq_from_d );
	clear_if_dash( cdr3_nucseq_from_j );
	clear_if_dash( v_ties );
	clear_if_dash( d_ties );
	clear_if_dash( j_ties );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is an all-purpose class used for collecting VDJ rearrangement statistics and also
// computing rearrangement probabilities.
//
class JunctionCounts {
public:

	JunctionCounts();

	bool // true if it was added to the counts
	store_junction(
		JunctionParse const & jp,
		string const & nucseq
	);

	void
	compute_probs();

	void
	show_probs( string const & tag, ostream & out ) const;

	void
	read_probs_and_counts_from_file(
		string const & filename,
		bool const save_old_counts
	);

	void
	resample_junction_and_read(
		Size const readlen,
		JunctionParse & jp,
		string & readseq
	) const;

	void
	initialize_counts_and_probs_arrays();

	void
	adjust_probs_to_move_sampled_toward_target(
		JunctionCounts const & target,
		JunctionCounts const & sampled
	);

	void
	to_probs_map( map< string, Reals > & pmap ) const;

	void
	from_probs_map( map< string, Reals > const & pmap );


	Real
	calc_ndn_nucseq_probability_given_j_gene(
		string const & ndn_nucseq,
		string const & j_gene,
		vector< Real > const & v_probs_by_prefix,
		vector< Real > const & j_probs_by_prefix // j_prefix is reverse-complemented wrt ndn_nucseq
	) const;

	Real
	calc_cdr3_nucseq_probability(
		string const & cdr3_nucseq,
		strings const & v_ties,
		strings const & j_ties
	) const;

	Real
	calc_cdr3_protseq_probability(
		string const & cdr3_nucseq,
		strings const & v_ties,
		strings const & j_ties
	) const;

	strings const &
	j_genes() const { return j_genes_; }

	strings const &
	v_genes() const { return v_genes_; }

	string
	resample_nucseq(
		Size const len,
		char const prefix
	) const;

	inline
	void
	increment_vj_trim_count(
		string const & g, // V or J determined from gene name
		int const ptrim
	);


	inline
	void
	increment_d_trim_count(
		string const & g,
		int const ptrim0,
		int const ptrim1
	);

	void
	increment_nuc_counts(
		char const prefix,
		string const & nucseq
	);

	Real
	nucseq_probability(
		string const & nucseq,
		char const prefix
	) const;

	inline
	Real
	gene_probability( string const & g ) const;

private:
	int max_trim_, max_vj_insert_, max_vdj_insert_, max_pnucs_;
	string alphabet_;

	strings v_genes_, j_genes_, d_genes_;

	map<string,string> v_fasta_, d_fasta_, j_fasta_, v_fasta_cdr3_region_, j_fasta_cdr3_region_,
		v_fasta_cdr3_region_inc_pnucs_, j_fasta_cdr3_region_inc_pnucs_, d_fasta_inc_pnucs_;

	map<string,int> all_num_cdr3_nucleotides_;

	// counts
	map< string, Size > v_counts_, d_counts_, j_counts_;
	map< string, vector< Size > > v_trim_counts_, j_trim_counts_, j_trail_counts_;
	map< string, vector< vector< Size > > > d_trim_counts_;
	vector< Size > vd_insert_counts_, dj_insert_counts_, vj_insert_counts_;  // need to start at 0
	vector< vector< Size > > nuc_counts_;

	// probabilities
	map< string, Real > v_probs_, d_probs_, j_probs_;
	map< string, vector< Real > > v_trim_probs_, j_trim_probs_, j_trail_probs_;
	map< string, vector< vector< Real > > > d_trim_probs_;
	vector< Real > vd_insert_probs_, dj_insert_probs_, vj_insert_probs_;  // need to start at 0
	vector< vector< Real > > nuc_probs_;

	Real d_success_rate_;

};


/// ctor
JunctionCounts::JunctionCounts():
	max_trim_( vdj_probs::max_trim ), // does not include pnucs
	max_vj_insert_( vdj_probs::max_vj_insert ),
	max_vdj_insert_( vdj_probs::max_vdj_insert ), // max insertion on either side of d gene
	max_pnucs_( vdj_probs::max_p_nucleotides ), // seems to capture most of them
	alphabet_( "acgt" )
{

	read_vj_fastas_and_info( v_fasta_, j_fasta_, all_num_cdr3_nucleotides_ );
	read_fasta( misc::dbdir+"human_D_nucleotide_sequences.fasta", d_fasta_ );

	v_genes_ = get_keys( v_fasta_ );
	d_genes_ = get_keys( d_fasta_ );
	j_genes_ = get_keys( j_fasta_ );
	std::sort( v_genes_.begin(), v_genes_.end() );
	std::sort( d_genes_.begin(), d_genes_.end() );
	std::sort( j_genes_.begin(), j_genes_.end() );

	foreach_( string g, v_genes_ ) {
		string const & fullseq( v_fasta_.find(g)->second );
		Size const sublen( all_num_cdr3_nucleotides_[g] );
		string const cdr3seq( fullseq.substr(fullseq.size() - sublen ) );
		v_fasta_cdr3_region_[ g ] = cdr3seq;
		runtime_assert( (int) fullseq.size() >= max_pnucs_ );
		string const pnucseq( reverse_complement( fullseq.substr( fullseq.size()-max_pnucs_ ) ) );
		v_fasta_cdr3_region_inc_pnucs_[g] = cdr3seq + pnucseq;
	}

	foreach_( string g, j_genes_ ) {
		string const & fullseq( j_fasta_.find(g)->second );
		Size const sublen( all_num_cdr3_nucleotides_[g] );
		string const cdr3seq( fullseq.substr( 0, sublen ) );
		runtime_assert( (int) fullseq.size() >= max_pnucs_ );
		j_fasta_cdr3_region_[ g ] = cdr3seq;
		j_fasta_cdr3_region_inc_pnucs_[ g ] = reverse_complement( fullseq.substr(0,max_pnucs_) ) + cdr3seq;
	}

	foreach_( string g, d_genes_ ) {
		string const fullseq( d_fasta_.find( g )->second );
		runtime_assert( (int) fullseq.size() >= max_pnucs_ );
		d_fasta_inc_pnucs_[g] = reverse_complement( fullseq.substr(0,max_pnucs_) ) + fullseq +
			reverse_complement( fullseq.substr( fullseq.size() - max_pnucs_ ) );
	}

	initialize_counts_and_probs_arrays();
}

void
JunctionCounts::initialize_counts_and_probs_arrays()
{
	Size const trim_array_size( max_trim_ + max_pnucs_ + 1 );

	/// size and initialize the counts arrays
	vd_insert_counts_.clear();
	dj_insert_counts_.clear();
	vj_insert_counts_.clear();
	vd_insert_counts_.resize( max_vdj_insert_+1, 0 );
	dj_insert_counts_.resize( max_vdj_insert_+1, 0 );
	vj_insert_counts_.resize( max_vj_insert_+1, 0 );

	foreach_ ( string const g, v_genes_ ) {
		v_counts_[g] = 0;
		v_trim_counts_[g];
		v_trim_counts_.find(g)->second.clear();
		v_trim_counts_.find(g)->second.resize( trim_array_size, 0 );
	}

	foreach_ ( string const g, j_genes_ ) {
		j_counts_[g] = 0;
		j_trim_counts_[g];
		j_trail_counts_[g];
		j_trim_counts_.find(g)->second.clear();
		j_trail_counts_.find(g)->second.clear();
		j_trim_counts_.find(g)->second.resize( trim_array_size, 0 );
		j_trail_counts_.find(g)->second.resize( j_fasta_.find(g)->second.size(), 0 );
	}

	foreach_ ( string const g, d_genes_ ) {
		d_counts_[g] = 0;
		Size const seqlen( d_fasta_.find(g)->second.size() );
		d_trim_counts_[g];
		d_trim_counts_.find(g)->second.clear();
		d_trim_counts_.find(g)->second.resize( seqlen+1+max_pnucs_ );
		for ( Size i=0; i< seqlen+1+max_pnucs_; ++i ) {
			d_trim_counts_.find(g)->second[ i ].clear();
			d_trim_counts_.find(g)->second[ i ].resize( seqlen+1+max_pnucs_, 0 );
		}
	}

	nuc_counts_.clear();
	nuc_counts_.resize( alphabet_.size() );
	for ( Size i=0; i<alphabet_.size(); ++i ) {
		nuc_counts_[i].resize( alphabet_.size(), 0 );
	}

	/// size and initialize the probs arrays (don't really need to initialize)
	//
	// don't be so paranoid about zeroing out here...
	//
	vd_insert_probs_.resize( max_vdj_insert_+1, 0 );
	dj_insert_probs_.resize( max_vdj_insert_+1, 0 );
	vj_insert_probs_.resize( max_vj_insert_+1, 0 );

	foreach_ ( string const g, v_genes_ ) {
		v_probs_[g] = 0;
		v_trim_probs_[g];
		v_trim_probs_.find(g)->second.resize( trim_array_size, 0 );
	}

	foreach_ ( string const g, j_genes_ ) {
		j_probs_[g] = 0;
		j_trim_probs_[g];
		j_trail_probs_[g];
		j_trim_probs_.find(g)->second.resize( trim_array_size, 0 );
		j_trail_probs_.find(g)->second.resize( j_fasta_.find(g)->second.size(), 0 );
	}

	foreach_ ( string const g, d_genes_ ) {
		d_probs_[g] = 0;
		Size const seqlen( d_fasta_.find(g)->second.size() );
		d_trim_probs_[g];
		d_trim_probs_.find(g)->second.resize( seqlen+1+max_pnucs_ );
		for ( Size i=0; i< seqlen+1+max_pnucs_; ++i ) {
			d_trim_probs_.find(g)->second[ i ].resize( seqlen+1+max_pnucs_, 0 );
		}
	}

	nuc_probs_.clear();
	nuc_probs_.resize( alphabet_.size() );
	for ( Size i=0; i<alphabet_.size(); ++i ) {
		nuc_probs_[i].resize( alphabet_.size(), 0 );
	}
}

inline
Real
JunctionCounts::gene_probability( string const & g ) const
{
	if ( g[3] == 'V' ) {
		runtime_assert( v_probs_.count(g) == 1 );
		return v_probs_.find(g)->second;
	} else if ( g[3] == 'J' ) {
		runtime_assert( j_probs_.count(g) == 1 );
		return j_probs_.find(g)->second;
	} else if ( g[3] == 'D' ) {
		runtime_assert( d_probs_.count(g) == 1 );
		return d_probs_.find(g)->second;
	}

	runtime_assert( false ); // bad gene
	return 0.0;
}

inline
void
JunctionCounts::increment_vj_trim_count(
	string const & g,
	int const ptrim // could be negative, to indicate p-nucs
)
{
	char const vj( g[3] );
	runtime_assert( vj == 'V' || vj == 'J' );
	map< string, vector< Size > > & counts( vj == 'V' ? v_trim_counts_ : j_trim_counts_ );
	runtime_assert( counts.find(g) != counts.end() );
	runtime_assert( ptrim + max_pnucs_ >= 0 );
	Size const trim( ptrim + max_pnucs_ );
	vector< Size > & gcounts( counts.find(g)->second );
	runtime_assert( trim < gcounts.size() );
	gcounts[ trim ] = gcounts[ trim ] + 1;
}


inline
void
JunctionCounts::increment_d_trim_count(
	string const & g,
	int const ptrim0,
	int const ptrim1
)
{
	runtime_assert( d_trim_counts_.find(g) != d_trim_counts_.end() );
	runtime_assert( ptrim0 + max_pnucs_ >= 0  && ptrim1 + max_pnucs_ >= 0 );
	Size const trim0( ptrim0 + max_pnucs_ );
	Size const trim1( ptrim1 + max_pnucs_ );
	vector< vector< Size > > & gcounts( d_trim_counts_.find(g)->second );
	runtime_assert( trim0 < gcounts.size() );
	runtime_assert( trim1 < gcounts[trim0].size() );
	gcounts[ trim0 ][ trim1 ] = gcounts[ trim0 ][ trim1 ] + 1;
}


void
JunctionCounts::increment_nuc_counts(
	char const prefix,
	string const & nucseq
)
{
	Size oldnuc( alphabet_.find(prefix) );
	runtime_assert( oldnuc != string::npos );
	foreach_( char const nuc, nucseq ) {
		Size const newnuc( alphabet_.find( nuc ) );
		runtime_assert( newnuc != string::npos );
		++nuc_counts_[ oldnuc ][ newnuc ];
		oldnuc = newnuc;
	}

}

/// assumes we are oriented with gene_nucseq coming first, then insertseq, like P-nucs after a V gene
///
/// if there are palindromic nucleotides we shift them from the beginning of insertseq to the end of gene_nucseq
///
///
void
get_palindromic_trim_adjust_nucseqs(
	Size const trim,
	string & gene_nucseq,
	string & insertseq,
	int & ptrim
)
{
	ptrim = trim;
	if ( trim > 0 || insertseq.empty() || gene_nucseq.empty() ) return; // can't be any pnucs

	Size const gene_nucseq_oldsize( gene_nucseq.size() ), insertseq_oldsize( insertseq.size() );

	while ( (!insertseq.empty()) && abs(ptrim) < int( vdj_probs::max_p_nucleotides ) ) {
		// first pnuc: gene_nucseq[ size-1 ]
		// second pnuc: gene_nucseq[ size-3 ] since we add the first pnuc to the end
		if ( gene_nucseq[ int(gene_nucseq.size()) - 1 + 2*ptrim ] == reverse_complement_nucleotide( insertseq[ 0 ] ) ) {
			gene_nucseq.push_back( insertseq[0] );
			insertseq.erase( insertseq.begin() );
			--ptrim;
		} else {
			break;
		}
	}

	runtime_assert( (int) gene_nucseq.size() == (int ) gene_nucseq_oldsize - ptrim );
	runtime_assert( (int) insertseq.size() == (int ) insertseq_oldsize + ptrim );

}

bool // true if it was added to the counts
JunctionCounts::store_junction(
	JunctionParse const & jp,
	string const & nucseq
)
{
	if ( jp.v_gene.empty() || jp.j_gene.empty() ) return false;

	if ( v_counts_.find( jp.v_gene ) == v_counts_.end() ||
		j_counts_.find( jp.j_gene ) == j_counts_.end() ) {
		cout << "bad genes? " << jp << endl;
		return false;
	}

	if ( jp.v_trim > max_trim_ || jp.j_trim > max_trim_ ) return false;

	if ( jp.d_gene.empty() ) {
		// failed to find d gene
		if ( jp.vj_insert > max_vj_insert_ ) return false;
	} else {
		if ( d_counts_.find( jp.d_gene ) == d_counts_.end() ) {
			cout << "bad D gene? " << jp << endl;
			return false;
		}
		if ( jp.vd_insert > max_vdj_insert_ || jp.dj_insert > max_vdj_insert_ ) return false;
	}


	Size const cdr3_offset( nucseq.find( jp.cdr3_nucseq ) );
	if ( cdr3_offset == string::npos ) {
		cout << "nucseq does not contain cdr3_nucseq: " << nucseq << ' ' << jp.cdr3_nucseq << endl;
		return false;
	}

	Size const ntrail( nucseq.size() - jp.cdr3_nucseq.size() - cdr3_offset );
	runtime_assert( ntrail < nucseq.size() );
	if ( ntrail >= j_trail_counts_.find( jp.j_gene )->second.size() ) {
		cout << "too many trailing nucs: " << nucseq << ' ' << jp.cdr3_nucseq << endl;
		return false;
	}


	// ok, now we are going to add to our counts ///////////////////////////////////////////////////////
	++v_counts_[ jp.v_gene ];
	++j_counts_[ jp.j_gene ];


	//if ( TR.Trace.visible() ) TR.Trace << "JunctionCounts::store_junction: " << jp << endl;


	// split logic by whether or not we have d gene
	int v_ptrim, j_ptrim, d0_ptrim, d1_ptrim;
	if ( jp.d_gene.empty() ) { //////////////////////////////// no D gene /////////////////////////////////////////////////
		runtime_assert( jp.vd_insertseq.empty() && jp.dj_insertseq.empty() );  // not '-'

		string vj_insertseq( jp.vj_insertseq ), cdr3_nucseq_from_v( jp.cdr3_nucseq_from_v );
		get_palindromic_trim_adjust_nucseqs( jp.v_trim, cdr3_nucseq_from_v, vj_insertseq, v_ptrim );

		string jv_insertseq( reverse_complement( vj_insertseq ) ),
			cdr3_nucseq_from_j_rev( reverse_complement( jp.cdr3_nucseq_from_j ) );

		get_palindromic_trim_adjust_nucseqs( jp.j_trim, cdr3_nucseq_from_j_rev, jv_insertseq, j_ptrim );

		vj_insertseq = reverse_complement( jv_insertseq ); // update other direction in case trimmed

		increment_vj_trim_count( jp.v_gene, v_ptrim );
		increment_vj_trim_count( jp.j_gene, j_ptrim );

		++vj_insert_counts_[ vj_insertseq.size() ];

		// get aa counts for vj_insertseq fwd and reverse directions
		// just cut it in half, give front to v and back to j
		string const v_insertseq( vj_insertseq.substr( 0, vj_insertseq.size()/2 ) ),
			j_insertseq_rev( reverse_complement( vj_insertseq.substr( vj_insertseq.size()/2 ) ) );

		if ( !( cdr3_nucseq_from_v.empty() || v_insertseq.empty() ) ) {
			increment_nuc_counts( cdr3_nucseq_from_v.back(), v_insertseq );
		}
		if ( !( cdr3_nucseq_from_j_rev.empty() || j_insertseq_rev.empty() ) ) {
			increment_nuc_counts( cdr3_nucseq_from_j_rev.back(), j_insertseq_rev );
		}

	} else { ////////////////////////////// found a D gene ////////////////////////////////////////////////////////////////
		string vd_insertseq( jp.vd_insertseq ), jd_insertseq( reverse_complement( jp.dj_insertseq ) ),
			cdr3_nucseq_from_v( jp.cdr3_nucseq_from_v ),
			cdr3_nucseq_from_j_rev( reverse_complement( jp.cdr3_nucseq_from_j ) );

		get_palindromic_trim_adjust_nucseqs( jp.v_trim, cdr3_nucseq_from_v    , vd_insertseq, v_ptrim );
		get_palindromic_trim_adjust_nucseqs( jp.j_trim, cdr3_nucseq_from_j_rev, jd_insertseq, j_ptrim );

		increment_vj_trim_count( jp.v_gene, v_ptrim );
		increment_vj_trim_count( jp.j_gene, j_ptrim );

		// now that vd_insertseq and jd_insertseq have been updated we can get the d ptrims
		string
			dv_insertseq( reverse_complement( vd_insertseq ) ),
			dj_insertseq( reverse_complement( jd_insertseq ) ),
			cdr3_nucseq_from_d( jp.cdr3_nucseq_from_d ), // these two cdr3_nucseq_from_d and rev will go out of sync
			cdr3_nucseq_from_d_rev( reverse_complement( jp.cdr3_nucseq_from_d ) ); // but they are not used...

		get_palindromic_trim_adjust_nucseqs( jp.d0_trim, cdr3_nucseq_from_d_rev, dv_insertseq, d0_ptrim );
		get_palindromic_trim_adjust_nucseqs( jp.d1_trim, cdr3_nucseq_from_d    , dj_insertseq, d1_ptrim );

		vd_insertseq = reverse_complement( dv_insertseq ); // update these two, to...
		jd_insertseq = reverse_complement( dj_insertseq ); // reflect removal of pnucs

		++d_counts_[ jp.d_gene ];
		increment_d_trim_count( jp.d_gene, d0_ptrim, d1_ptrim ); //, d_trim_counts_ );

		/// be sure to use the most recent insertseqs here
		++vd_insert_counts_[ dv_insertseq.size() ];
		++dj_insert_counts_[ dj_insertseq.size() ];

		if ( !( cdr3_nucseq_from_v.empty() || vd_insertseq.empty() ) ) {
			increment_nuc_counts( cdr3_nucseq_from_v.back(), vd_insertseq );
		}
		if ( !( cdr3_nucseq_from_j_rev.empty() || jd_insertseq.empty() ) ) {
			increment_nuc_counts( cdr3_nucseq_from_j_rev.back(), jd_insertseq );
		}

	}


	// store info on length of trailing seqs
	//increment_vj_trim_count( jp.j_gene, ntrail, j_trail_counts_ );
	++j_trail_counts_.find( jp.j_gene )->second[ ntrail ];

	return true;
}


void
JunctionCounts::compute_probs()
{
	Size v_total(0), j_total(0), d_total(0);
	foreach_ ( string g, v_genes_ ) { v_total += v_counts_[g]; }
	foreach_ ( string g, d_genes_ ) { d_total += d_counts_[g]; }
	foreach_ ( string g, j_genes_ ) { j_total += j_counts_[g]; }

	runtime_assert( v_total == j_total );

	d_success_rate_ = Real(d_total) / v_total;

	foreach_( string g, v_genes_ ) {
		v_probs_[g] = Real( v_counts_[g] ) / v_total;
		Size total( v_counts_[g] ), retotal(0);
		vector<Real> & trim_probs( v_trim_probs_.find(g)->second );
		vector<Size> const & trim_counts( v_trim_counts_.find(g)->second );
		runtime_assert( trim_probs.size() == trim_counts.size() );
		if ( total ) {
			for ( Size i=0; i< trim_counts.size(); ++i ) {
				trim_probs[i] = Real( trim_counts[i] )/total;
				retotal += trim_counts[i];
			}
			runtime_assert( retotal == total );
		} else {
			std::fill( trim_probs.begin(), trim_probs.end(), 0.0 );
		}
	}

	foreach_( string g, j_genes_ ) {
		j_probs_[g] = Real( j_counts_[g] ) / j_total;
		Size total( j_counts_[g] ), retotal(0);
		vector<Real> &trim_probs( j_trim_probs_.find(g)->second ), &trail_probs( j_trail_probs_.find(g)->second );
		vector<Size> const & trim_counts( j_trim_counts_.find(g)->second ),
			&trail_counts( j_trail_counts_.find( g )->second );
		runtime_assert( trim_counts.size() == trim_probs.size() );
		if ( total ) {
			for ( Size i=0; i< trim_counts.size(); ++i ) {
				trim_probs[i] = Real( trim_counts[i] )/total;
				retotal += trim_counts[i];
			}
			runtime_assert( retotal == total );
			for ( Size i=0; i< trail_counts.size(); ++i ) {
				trail_probs[i] = Real( trail_counts[i] )/total;
			}
		} else {
			std::fill( trim_probs.begin(), trim_probs.end(), 0.0 );
			std::fill( trail_probs.begin(), trail_probs.end(), 0.0 );
		}
	}

	foreach_( string g, d_genes_ ) {
		d_probs_[g] = Real( d_counts_[g] ) / d_total;
		Size total( d_counts_[g] ), retotal(0);
		vector< vector<Real> > & trim_probs( d_trim_probs_.find(g)->second );
		vector< vector<Size> > const & trim_counts( d_trim_counts_.find(g)->second );
		Size const sz( trim_counts.size() );
		runtime_assert( trim_probs.size() == sz );
		for ( Size i=0; i< sz; ++i ) {
			runtime_assert( trim_counts[i].size() == sz );
			runtime_assert( trim_probs[i].size() == sz );
			for ( Size j=0; j< sz; ++j ) {
				if ( total ) {
					trim_probs[i][j] = Real( trim_counts[i][j] )/total;
					retotal += trim_counts[i][j];
				} else {
					trim_probs[i][j] = 0.0;
				}
			}
		}
		runtime_assert( retotal == total );
	}

	// now inserts
	{
		Size total( v_total - d_total ), retotal(0);

		for ( int i=0; i<= max_vj_insert_; ++i ) {
			vj_insert_probs_[i] = Real(vj_insert_counts_[i]) / total;
			retotal += vj_insert_counts_[i];
		}
		runtime_assert( total == retotal );
	}

	{
		Size total( d_total ), retotal(0);

		for ( int i=0; i<= max_vdj_insert_; ++i ) {
			vd_insert_probs_[i] = Real(vd_insert_counts_[i]) / total;
			retotal += vd_insert_counts_[i];
		}
		runtime_assert( total == retotal );

		retotal=0;
		for ( int i=0; i<= max_vdj_insert_; ++i ) {
			dj_insert_probs_[i] = Real(dj_insert_counts_[i]) / total;
			retotal += dj_insert_counts_[i];
		}
		runtime_assert( total == retotal );
	}

	// now nucleotide counts
	runtime_assert( nuc_counts_.size() == alphabet_.size() );
	runtime_assert( nuc_probs_.size() == alphabet_.size() );
	for ( Size i=0; i< alphabet_.size(); ++i ) { // compute probs conditioned on alphabet_[i] as 1st na
		Size total(0);
		for ( Size j=0; j< alphabet_.size(); ++j ) total += nuc_counts_[i][j];
		for ( Size j=0; j< alphabet_.size(); ++j ) {
			nuc_probs_[i][j] = Real( nuc_counts_[i][j] ) / total;
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// take ratio of target divided by sampled, multiply probs by this
void
adjust_sampling_probs(
	Reals const & target,
	Reals const & sampled,
	Reals & probs, //  used to make the distn seen in sampled
	string const tag = string()
)
{
	Real const sampling_adjustment_factor( 1.0 ); // -sampling_adjustment_factor in tmp.compare_probs.v[456].log
	Size const n( probs.size() );
	runtime_assert( target.size() == n );
	runtime_assert( sampled.size() == n );

	Real total(0), totdev(0), jsd(0.);
	for ( Size i=0; i< n; ++i ) {
		Real oldvalue( probs[i] ), newvalue( probs[i] ); // default is to reuse old value
		if ( sampled[i] < 1e-9 ) {
			if ( target[i] > 1e-9 ) {
				newvalue = max( oldvalue, target[i] );
			}
		} else {
			Real const ratio( target[i] / sampled[i] ),
				rescale( exp( sampling_adjustment_factor * log ( max( 0.5, min( 2.0, ratio ) ) ) ) );
			newvalue = rescale * oldvalue;
		}
		probs[i] = newvalue;
		total += newvalue;
		if ( !tag.empty() ) {
			cout << "adjust_sampling_probs: " << I(3,i) <<
				" target: " << F(9,6,target[i]) <<
				" sampled: " << F(9,6,sampled[i]) <<
				" old_probs: " << F(9,6,oldvalue) <<
				" new_probs: " << F(9,6,newvalue) << ' ' << tag << endl;
		}
		{ // diagnostics
			totdev += fabs( target[i] - sampled[i] );
			Real const avgprob( 0.5 *( target[i] + sampled[i] ) ), tinyprob(1e-9);
			if ( avgprob > tinyprob ) {
				if (  target[i] > tinyprob ) jsd += 0.5 * (  target[i] * log(  target[i]/avgprob ) / log(2.0) );
				if ( sampled[i] > tinyprob ) jsd += 0.5 * ( sampled[i] * log( sampled[i]/avgprob ) / log(2.0) );
			}
		}
	}
	if ( !tag.empty() ) {
		cout << "adjust_sampling_probs: totdev: " << F(12,6,totdev) << " jsd: " << F(12,6,jsd) << ' ' << tag << endl;
	}

	if ( total >1e-9 ) { // now normalize
		for ( Size i=0; i< n; ++i ) {
			probs[i] /= total;
		}
	}
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// after this, counts will be invalid
//
void
JunctionCounts::adjust_probs_to_move_sampled_toward_target(
	JunctionCounts const & target,
	JunctionCounts const & sampled
)
{
	map< string, Reals > target_pmap, sampled_pmap, pmap;
	target.to_probs_map( target_pmap );
	sampled.to_probs_map( sampled_pmap );
	this->to_probs_map( pmap );

	strings const tags( get_keys( pmap ) );
	foreach_( string tag, tags ) {
		Reals const & target_probs( target_pmap.find( tag )->second );
		Reals const & sampled_probs( sampled_pmap.find( tag )->second );
		Reals & probs( pmap.find( tag )->second );
		string cout_tag( tag );
		if ( tag.find("_trim_probs:")!=string::npos ) { // better status info
			ostringstream out;
			strings const l( split_to_vector(tag,":") );
			runtime_assert( l.size() == 2 );
			string const g( l.front() );
			out << " target_gprob: " << F(9,6,target.gene_probability(g) ) <<
				" observed_gprob: " << F(9,6,sampled.gene_probability(g) ) <<
				" sampling_gprob: " << F(9,6,this->gene_probability(g)) <<
				' ' << tag;
			cout_tag = out.str();
		}
		adjust_sampling_probs( target_probs, sampled_probs, probs, cout_tag );
	}

	this->from_probs_map( pmap );
}




void
JunctionCounts::show_probs(
	string const & tag,
	ostream & out
) const
{
	foreach_( string g, v_genes_ ) {
		out << "v_probs: " << F(12,9,v_probs_.find(g)->second) << ' ' << g << ' ' << v_counts_.find(g)->second << ' ' <<
			tag << '\n';
		vector<Real> const & trim_probs( v_trim_probs_.find(g)->second );
		vector<Size> const & trim_counts( v_trim_counts_.find(g)->second );
		for ( Size i=0; i< trim_counts.size(); ++i ) {
			out << "v_trim_probs: " << I(2,i) << F(12,9,trim_probs[i]) << ' ' << g << ' ' << trim_counts[i] << ' ' << tag <<
				'\n';
		}
	}

	foreach_( string g, j_genes_ ) {
		out << "j_probs: " << F(12,9,j_probs_.find(g)->second) << ' ' << g << ' ' << j_counts_.find(g)->second << ' ' <<
			tag << '\n';
		vector<Real> const & trim_probs( j_trim_probs_.find(g)->second );
		vector<Size> const & trim_counts( j_trim_counts_.find(g)->second );
		for ( Size i=0; i< trim_probs.size(); ++i ) {
			out << "j_trim_probs: " << I(2,i) << F(12,9,trim_probs[i]) << ' ' << g << ' ' << trim_counts[i] << ' ' << tag <<
				'\n';
		}
		vector<Real> const & trail_probs( j_trail_probs_.find(g)->second );
		vector<Size> const & trail_counts( j_trail_counts_.find(g)->second );
		for ( Size i=0; i< trail_counts.size(); ++i ) {
			out << "j_trail_probs: " << I(2,i) << F(12,9,trail_probs[i]) << ' ' << g << ' ' << trail_counts[i] << ' ' << tag <<
				'\n';
		}
	}

	foreach_( string g, d_genes_ ) {
		out << "d_probs: " << F(12,9,d_probs_.find(g)->second) << ' ' << g << ' ' << d_counts_.find(g)->second << ' ' <<
			tag << '\n';
		vector< vector<Real> > const & trim_probs( d_trim_probs_.find(g)->second );
		vector< vector<Size> > const & trim_counts( d_trim_counts_.find(g)->second );
		Size const seqlen( trim_counts.size()-1 );
		for ( Size i=0; i<= seqlen; ++i ) {
			for ( Size j=0; j<= seqlen; ++j ) {
				out << "d_trim_probs: " << I(3,seqlen) << I(3,i) << I(3,j) << F(12,9,trim_probs[i][j]) << ' ' << g << ' ' <<
					trim_counts[i][j] << ' ' << tag << '\n';
			}
		}
	}

	// now inserts
	for ( int i=0; i<= max_vj_insert_; ++i ) {
		out << "vj_insert_probs: " << I(2,i) << F(12,9,vj_insert_probs_[i] ) << ' ' << vj_insert_counts_[i] << ' ' << tag <<
			'\n';
	}
	for ( int i=0; i<= max_vdj_insert_; ++i ) {
		out << "vd_insert_probs: " << I(2,i) << F(12,9,vd_insert_probs_[i] ) << ' ' << vd_insert_counts_[i] << ' ' << tag <<
			'\n';
	}
	for ( int i=0; i<= max_vdj_insert_; ++i ) {
		out << "dj_insert_probs: " << I(2,i) << F(12,9,dj_insert_probs_[i] ) << ' ' << dj_insert_counts_[i] << ' ' << tag <<
			'\n';
	}

	// nuc_probs_
	for ( Size i=0; i< alphabet_.size(); ++i ) {
		for ( Size j=0; j< alphabet_.size(); ++j ) {
			out << "nuc_probs: " << i << ' ' << j << ' ' << F(12,9,nuc_probs_[i][j]) << ' ' << nuc_counts_[i][j] << ' ' <<
				tag << '\n';
		}
	}



	out << "d_success_rate: " << F(12,9,d_success_rate_) << ' ' << tag << '\n';

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// doesn't have counts information
//
void
JunctionCounts::to_probs_map(
	map< string, Reals > & pmap
) const
{

	pmap.clear();

	foreach_( string g, v_genes_ ) pmap["v_probs"].push_back( v_probs_.find(g)->second );
	foreach_( string g, d_genes_ ) pmap["d_probs"].push_back( d_probs_.find(g)->second );
	foreach_( string g, j_genes_ ) pmap["j_probs"].push_back( j_probs_.find(g)->second );

	// trims
	foreach_( string g, v_genes_ ) {
		string const tag( "v_trim_probs:"+g );
		vector<Real> const & trim_probs( v_trim_probs_.find(g)->second );
		for ( Size i=0; i< trim_probs.size(); ++i ) pmap[tag].push_back( trim_probs[i] );
	}

	foreach_( string g, j_genes_ ) {
		string tag( "j_trim_probs:"+g );
		vector<Real> const & trim_probs( j_trim_probs_.find(g)->second );
		for ( Size i=0; i< trim_probs.size(); ++i ) pmap[tag].push_back( trim_probs[i] );
		tag = "j_trail_probs:"+g;
		Size const jlen( j_fasta_.find(g)->second.size() );
		runtime_assert( j_trail_probs_.find(g)->second.size() == jlen );
		for ( Size i=0; i< jlen; ++i ) pmap[tag].push_back( j_trail_probs_.find(g)->second[i] );
	}


	foreach_( string g, d_genes_ ) {
		string const tag( "d_trim_probs:"+g );
		vector< vector<Real> > const & trim_probs( d_trim_probs_.find(g)->second );
		Size const sz( trim_probs.size() );
		for ( Size i=0; i< sz; ++i ) {
			for ( Size j=0; j< sz; ++j ) {
				pmap[tag].push_back(trim_probs[i][j]);
			}
		}
	}

	// now inserts
	for ( int i=0; i<= max_vj_insert_; ++i ) pmap[ "vj_insert_probs" ].push_back( vj_insert_probs_[i] );
	for ( int i=0; i<= max_vdj_insert_; ++i ) pmap[ "vd_insert_probs" ].push_back( vd_insert_probs_[i] );
	for ( int i=0; i<= max_vdj_insert_; ++i ) pmap[ "dj_insert_probs" ].push_back( dj_insert_probs_[i] );

	// nuc probs
	for ( Size i=0; i< alphabet_.size(); ++i ) {
		string const tag( "nuc_probs_"+alphabet_.substr(i,1));
		for ( Size j=0; j< alphabet_.size(); ++j ) {
			pmap[ tag ].push_back( nuc_probs_[i][j] );
		}
	}

	pmap["d_success"] = make_vector( 1.0 - d_success_rate_, d_success_rate_ );

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// in Rosetta version Reals was 1-indexed, need to check this carefully!
//
void
JunctionCounts::from_probs_map(
	map< string, Reals > const & pmap
)
{
	initialize_counts_and_probs_arrays(); // now the counts will be zeroed out-- signal that they are invalid

	for ( Size i=0; i< v_genes_.size(); ++i ) v_probs_[ v_genes_[i] ] = pmap.find("v_probs")->second[i];
	for ( Size i=0; i< d_genes_.size(); ++i ) d_probs_[ d_genes_[i] ] = pmap.find("d_probs")->second[i];
	for ( Size i=0; i< j_genes_.size(); ++i ) j_probs_[ j_genes_[i] ] = pmap.find("j_probs")->second[i];

	// trims
	foreach_( string g, v_genes_ ) {
		string const tag( "v_trim_probs:"+g );
		vector<Real> & trim_probs( v_trim_probs_.find(g)->second );
		for ( Size i=0; i< trim_probs.size(); ++i ) trim_probs[i] = pmap.find(tag)->second[i];
	}

	foreach_( string g, j_genes_ ) {
		string tag( "j_trim_probs:"+g );
		vector<Real> & trim_probs( j_trim_probs_.find(g)->second );
		for ( Size i=0; i< trim_probs.size(); ++i ) trim_probs[i] = pmap.find(tag)->second[i];
		tag = "j_trail_probs:"+g;
		Size const jlen( j_fasta_.find(g)->second.size() );
		runtime_assert( j_trail_probs_.find(g)->second.size() == jlen );
		for ( Size i=0; i< jlen; ++i ) j_trail_probs_.find(g)->second[i] = pmap.find(tag)->second[i];
	}


	foreach_( string g, d_genes_ ) {
		string const tag( "d_trim_probs:"+g );
		vector< vector<Real> > & trim_probs( d_trim_probs_.find(g)->second );
		Size const sz( trim_probs.size() );
		Size counter(0);
		for ( Size i=0; i< sz; ++i ) {
			for ( Size j=0; j< sz; ++j ) {
				trim_probs[i][j] = pmap.find(tag)->second[ counter ];
				++counter;
			}
		}
	}

	// now inserts
	for ( int i=0; i<= max_vj_insert_ ; ++i ) vj_insert_probs_[i] = pmap.find("vj_insert_probs")->second[i];
	for ( int i=0; i<= max_vdj_insert_; ++i ) vd_insert_probs_[i] = pmap.find("vd_insert_probs")->second[i];
	for ( int i=0; i<= max_vdj_insert_; ++i ) dj_insert_probs_[i] = pmap.find("dj_insert_probs")->second[i];

	// nuc probs
	{
		for ( Size i=0; i< alphabet_.size(); ++i ) {
			string const tag( "nuc_probs_"+alphabet_.substr(i,1));
			for ( Size j=0; j< alphabet_.size(); ++j ) {
				nuc_probs_[i][j] = pmap.find( tag )->second[ j ];
			}
		}
	}

	d_success_rate_ = pmap.find("d_success")->second.back();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Note that if save_old_counts is TRUE -- the counts and probs at the end may be out of sync, so if you are
// going to read and sum counts from multiple files, at the end you had better recompute the probs...
//
//
void
JunctionCounts::read_probs_and_counts_from_file(
	string const & filename,
	bool const save_old_counts = false
)
{
	if ( verbose() ) cout << "JunctionCounts::read_probs_and_counts_from_file: " << filename << endl;
	ifstream data( filename.c_str() );

	if ( !save_old_counts ) initialize_counts_and_probs_arrays();

	string line, tag;
	Size num_success_rate_lines( 0 );
	while ( getline( data, line ) ) {
		istringstream l( line );
		l >> tag;
		Real p;
		string g;
		Size ii, jj, c, seqlen;
		if ( tag == "v_probs:" ) {
			l >> p >> g >> c;
			v_counts_[g] += c;
			v_probs_[g] = p;
		} else if ( tag == "j_probs:" ) {
			l >> p >> g >> c;
			j_counts_[g] += c;
			j_probs_[g] = p;
		} else if ( tag == "d_probs:" ) {
			l >> p >> g >> c;
			d_counts_[g] += c;
			d_probs_[g] = p;
		} else if ( tag == "v_trim_probs:" ) {
			l >> ii >> p >> g >> c;
			v_trim_counts_.find(g)->second[ii] += c;
			v_trim_probs_.find(g)->second[ii] = p;
		} else if ( tag == "j_trim_probs:" ) {
			l >> ii >> p >> g >> c;
			j_trim_counts_.find(g)->second[ii] += c;
			j_trim_probs_.find(g)->second[ii] = p;
		} else if ( tag == "j_trail_probs:" ) {
			l >> ii >> p >> g >> c;
			j_trail_counts_.find(g)->second[ii] += c;
			j_trail_probs_.find(g)->second[ii] = p;
		} else if ( tag == "d_trim_probs:" ) {
			l >> seqlen >> ii >> jj >> p >> g >> c;
			d_trim_counts_.find(g)->second[ii][jj] += c;
			d_trim_probs_.find(g)->second[ii][jj] = p;
		} else if ( tag == "vj_insert_probs:" ) {
			l >> ii >> p >> c;
			vj_insert_counts_[ii] += c;
			vj_insert_probs_[ii] = p;
		} else if ( tag == "vd_insert_probs:" ) {
			l >> ii >> p >> c;
			vd_insert_counts_[ii] += c;
			vd_insert_probs_[ii] = p;
		} else if ( tag == "dj_insert_probs:" ) {
			l >> ii >> p >> c;
			dj_insert_counts_[ii] += c;
			dj_insert_probs_[ii] = p;
		} else if ( tag == "d_success_rate:" ) {
			l >> d_success_rate_;
			++num_success_rate_lines;
		} else if ( tag == "nuc_probs:" ) {
			l >> ii >> jj >> p >> c;
			nuc_probs_[ii][jj] = p;
			nuc_counts_[ii][jj] += c;
		}
		runtime_assert( !l.fail() );
	} // while getline

	data.close();

	if ( num_success_rate_lines != 1 ) {
		cout << "ERROR bad number of d_success_rate_ lines: " << num_success_rate_lines << ' ' << filename << endl;
		cerr << "ERROR bad number of d_success_rate_ lines: " << num_success_rate_lines << ' ' << filename << endl;
		runtime_assert( false );
	}

	//compute_probs(); // based on the counts...

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// string
// random_nucseq( Size const len )
// {
// 	string const bases( "acgt" );
// 	string seq;
// 	for ( Size i=0; i<len; ++i ) {
// 		seq += bases[ random_range(0,bases.size()-1) ];
// 	}
// 	return seq;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string
JunctionCounts::resample_nucseq(
	Size const len,
	char const prefix
) const
{
	string seq;
	Size oldnuc( alphabet_.find(prefix));
	runtime_assert( oldnuc != string::npos );

	for ( Size i=0; i<len; ++i ) {
		char const newna( choose_random_char( nuc_probs_[ oldnuc ], alphabet_ ) );
		seq.push_back( newna );
		oldnuc = alphabet_.find( newna );
		runtime_assert( oldnuc != string::npos );
	}
	return seq;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
JunctionCounts::resample_junction_and_read(
	Size const readlen,
	JunctionParse & jp,
	string & readseq
) const
{
	using boost::to_upper_copy;

	jp.v_gene = choose_random_string( v_probs_ );
	jp.j_gene = choose_random_string( j_probs_ );

	int const v_ptrim( choose_random_int( v_trim_probs_.find( jp.v_gene )->second ) - max_pnucs_ );
	int const j_ptrim( choose_random_int( j_trim_probs_.find( jp.j_gene )->second ) - max_pnucs_ );
	int const v_pnucs( v_ptrim < 0 ? abs( v_ptrim ) : 0 ), j_pnucs( j_ptrim < 0 ? abs( j_ptrim ) : 0 );

	jp.v_trim = ( v_ptrim <= 0 ? 0 : v_ptrim );
	jp.j_trim = ( j_ptrim <= 0 ? 0 : j_ptrim );

	int const num_cdr3_nucleotides_left_from_v( all_num_cdr3_nucleotides_.find(jp.v_gene)->second - jp.v_trim );
	int const num_cdr3_nucleotides_left_from_j( all_num_cdr3_nucleotides_.find(jp.j_gene)->second - jp.j_trim );
	runtime_assert( num_cdr3_nucleotides_left_from_v >= 0 );
	runtime_assert( num_cdr3_nucleotides_left_from_j >= 0 );

	string const & v_seq( v_fasta_.find( jp.v_gene )->second ), &j_seq( j_fasta_.find( jp.j_gene )->second );
	jp.cdr3_nucseq_from_v = v_seq.substr( v_seq.size() - num_cdr3_nucleotides_left_from_v - jp.v_trim,
		num_cdr3_nucleotides_left_from_v );
	jp.cdr3_nucseq_from_j = j_seq.substr( jp.j_trim, num_cdr3_nucleotides_left_from_j );

	string const v_pnucseq( v_pnucs > 0 ? reverse_complement( jp.cdr3_nucseq_from_v ).substr(0,v_pnucs ): string("") );
	string const j_pnucseq( j_pnucs > 0 ? reverse_complement( jp.cdr3_nucseq_from_j.substr(0,j_pnucs ) ): string("") );

	// char const v_prefix( v_pnucs>0 ? *v_pnucseq.rbegin() : *jp.cdr3_nucseq_from_v.rbegin() );
	char const v_prefix( v_pnucs > 0 ? v_pnucseq.back() :
		( num_cdr3_nucleotides_left_from_v > 0 ? jp.cdr3_nucseq_from_v.back() :
			v_seq[ v_seq.size() - all_num_cdr3_nucleotides_.find(jp.v_gene)->second - 1 ] ) );

	char const j_prefix
		( reverse_complement_nucleotide(
			( j_pnucs > 0 ? j_pnucseq[0] :
				( num_cdr3_nucleotides_left_from_j > 0 ? jp.cdr3_nucseq_from_j[0] :
					j_seq[ all_num_cdr3_nucleotides_.find(jp.j_gene)->second ] ) ) ) );


	// d gene-- did we find one?
	bool const d_success( uniform() < d_success_rate_ );
	string ndn_nucseq, ndn_newseq;
	if ( d_success ) {
		// remember the D/J interaction...
		while ( true ) {
			jp.d_gene = choose_random_string( d_probs_ );
			Size const jno( int_of( jp.j_gene.substr(4,1) ) );
			Size const dno( int_of( jp.d_gene.substr(4,1) ) );
			runtime_assert( ( jno==1 || jno==2 ) && ( dno==1 || dno==2 ) );
			if ( jno == 2 || dno == 1 ) break;
		}

		int d0_ptrim, d1_ptrim;
		choose_random_int_pair( d_trim_probs_.find( jp.d_gene )->second, d0_ptrim, d1_ptrim );
		d0_ptrim -= max_pnucs_;
		d1_ptrim -= max_pnucs_;
		int const d0_pnucs( d0_ptrim < 0 ? abs( d0_ptrim ) : 0 ), d1_pnucs( d1_ptrim < 0 ? abs( d1_ptrim ) : 0 );
		jp.d0_trim = ( d0_ptrim < 0 ? 0 : d0_ptrim );
		jp.d1_trim = ( d1_ptrim < 0 ? 0 : d1_ptrim );
		string const d_seq( d_fasta_.find( jp.d_gene )->second );
		Size const dlen( d_seq.size() - (jp.d0_trim+jp.d1_trim) );
		runtime_assert( dlen >= vdj_probs::min_d_trimmed_length );

		jp.cdr3_nucseq_from_d = d_seq.substr( jp.d0_trim, dlen );

		string const d0_pnucseq( d0_pnucs > 0 ? reverse_complement( jp.cdr3_nucseq_from_d.substr(0,d0_pnucs ) ): string(""));
		string const d1_pnucseq( d1_pnucs > 0 ? reverse_complement( jp.cdr3_nucseq_from_d ).substr(0,d1_pnucs ): string(""));


		// junction parse doesn't know about p-nucleotides, so pnucs count as inserted nucleotides for jp
		//
		int const vd_insert_wo_pnucs( choose_random_int( vd_insert_probs_ ) );
		int const dj_insert_wo_pnucs( choose_random_int( dj_insert_probs_ ) );

		jp.vd_insertseq =  v_pnucseq + resample_nucseq( vd_insert_wo_pnucs, v_prefix ) + d0_pnucseq;
		jp.dj_insertseq = d1_pnucseq + reverse_complement( resample_nucseq( dj_insert_wo_pnucs, j_prefix ) ) + j_pnucseq;

		ndn_nucseq =                jp.vd_insertseq   + jp.cdr3_nucseq_from_d +                jp.dj_insertseq;
		ndn_newseq = to_upper_copy( jp.vd_insertseq ) + jp.cdr3_nucseq_from_d + to_upper_copy( jp.dj_insertseq );


		jp.vd_insert = jp.vd_insertseq.size();
		jp.dj_insert = jp.dj_insertseq.size();



		// zero out others
		jp.vj_insert = 0;
		jp.vj_insertseq.clear();

	} else {
		jp.d_gene.clear();

		int const vj_insert_wo_pnucs( choose_random_int( vj_insert_probs_ ) ),
			v_insert( vj_insert_wo_pnucs/2 ), j_insert( vj_insert_wo_pnucs - v_insert );

		jp.vj_insertseq = v_pnucseq + resample_nucseq( v_insert, v_prefix ) +
			reverse_complement( resample_nucseq( j_insert, j_prefix ) ) + j_pnucseq;
		jp.vj_insert = jp.vj_insertseq.size();

		ndn_nucseq = jp.vj_insertseq;
		ndn_newseq = to_upper_copy( jp.vj_insertseq );

		// zero out the other stuff
		jp.d0_trim = jp.d1_trim = jp.vd_insert = jp.dj_insert = 0;
		jp.vd_insertseq.clear();
		jp.dj_insertseq.clear();
		jp.cdr3_nucseq_from_d.clear();
	}

	// finalize the cdr3_nucseq
	//
	jp.cdr3_nucseq = jp.cdr3_nucseq_from_v + ndn_nucseq + jp.cdr3_nucseq_from_j;
	jp.cdr3_newseq = jp.cdr3_nucseq_from_v + ndn_newseq + jp.cdr3_nucseq_from_j;

	// now finalize the read sequence
	Size const num_after( choose_random_int( j_trail_probs_.find( jp.j_gene )->second ) );
	string const fullseq(
		v_seq.substr(0,v_seq.size() - all_num_cdr3_nucleotides_.find(jp.v_gene)->second ) +
		jp.cdr3_nucseq +
		j_seq.substr( all_num_cdr3_nucleotides_.find( jp.j_gene )->second, num_after ) );
	runtime_assert( readlen <= fullseq.size() );

	readseq = fullseq.substr( fullseq.size() - readlen );
	runtime_assert( readseq.size() == readlen );

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
JunctionCounts::nucseq_probability(
	string const & nucseq,
	char const prefix
) const
{
	if ( nucseq.empty() ) return 1.0;

	Size const sz_alpha( alphabet_.size() );
	vector<Real> prefix_probs( sz_alpha, 0.0 ), probs( sz_alpha, 0.0 );

	Size const oldnuc( alphabet_.find( prefix ) );
	runtime_assert( oldnuc != string::npos );

	prefix_probs[ oldnuc ] = 1.0;

	Real totalprob(1.0);

	foreach_( char c, nucseq ) {
		// update current probs array
		Real cprob(0);
		for ( Size i=0; i<sz_alpha; ++i ) {
			probs[i] = 0.0;
			if ( degnuc_matches_nuc( c, alphabet_[i] ) ) {
				for ( Size j=0; j<sz_alpha; ++j ) {
					probs[i] += prefix_probs[j] * nuc_probs_[j][i]; // prob of i given j
				}
				cprob += probs[i];
			}
		}
		totalprob *= cprob;
		// update prefix_probs
		prefix_probs.swap( probs );
		runtime_assert( cprob>1e-6 ); // nuc_probs_ should be pretty flat, right?
		for ( Size i=0; i<sz_alpha; ++i ) prefix_probs[i] /= cprob;
	}

	return totalprob;
}

/// this is for human...
bool
dj_incompatible( string const & d_gene, string const & j_gene )
{
	return ( d_gene[4] != '1' && j_gene[4] == '1' );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
JunctionCounts::calc_ndn_nucseq_probability_given_j_gene(
	string const & ndn_nucseq,
	string const & j_gene,
	vector< Real > const & v_probs_by_prefix,
	vector< Real > const & j_probs_by_prefix // j_prefix is reverse-complemented wrt ndn_nucseq
) const
{
	Size const sz_alpha( alphabet_.size() );
	runtime_assert( v_probs_by_prefix.size() == sz_alpha );
	runtime_assert( j_probs_by_prefix.size() == sz_alpha );
	int const max_extra_d_trim( vdj_probs::max_extra_d_trim_for_probs );
	int const max_d_gap( vdj_probs::max_d_gap );

	Real prob(0);

	static UngappedAlignment al;

	int const ndn_len( ndn_nucseq.size() );
	runtime_assert( ndn_len>0 ); // handling that case outside

	// add prob that there was 'no D gene' in the sampling step
	if ( ndn_len <= max_vj_insert_ ) {
		string const v_insert( ndn_nucseq.substr(0,ndn_len/2) ),
			j_insert( reverse_complement( ndn_nucseq.substr(ndn_len/2) ) );
		Real prob_v(0.0);
		for ( Size v_prefix=0; v_prefix< sz_alpha; ++v_prefix ) {
			if ( v_probs_by_prefix[ v_prefix ]==0 ) continue;
			prob_v += v_probs_by_prefix[v_prefix] * nucseq_probability( v_insert, alphabet_[ v_prefix ] );
		}
		Real prob_j(0.0);
		for ( Size j_prefix=0; j_prefix< sz_alpha; ++j_prefix ) {
			if ( j_probs_by_prefix[ j_prefix ]==0 ) continue;
			prob_j += j_probs_by_prefix[j_prefix] * nucseq_probability( j_insert, alphabet_[ j_prefix ] );
		}
		prob += (1-d_success_rate_) * prob_v * prob_j * vj_insert_probs_[ ndn_len ];
		if ( verbose() ) cout << "prob_no_D: " << prob << ' ' << prob_v << ' ' << prob_j << ' ' <<
											 vj_insert_probs_[ndn_len] << endl;
	}

	// add prob that there was a D
	string vd_insertseq, jd_insertseq;
	if ( ndn_len >= (int) vdj_probs::min_d_trimmed_length ) {
		foreach_( string const & d_gene, d_genes_ ) {
			if ( dj_incompatible( d_gene, j_gene ) ) continue;
			Real prob_with_this_d(0), prob_d_gene( d_probs_.find(d_gene)->second ); //, prob_with_these_trims(0);
			string const & d_seq_inc_pnucs( d_fasta_inc_pnucs_.find(d_gene)->second );
			int const minlen( vdj_probs::min_d_trimmed_length ), d_gene_length( d_seq_inc_pnucs.size() - 2*max_pnucs_ );
			runtime_assert( (int) d_fasta_.find(d_gene)->second.size() == d_gene_length ); // sanity
			int best_score(0);
			for ( int best_score_gap=0; best_score_gap<= max_d_gap; ++best_score_gap ) {
				int const max_score( best_score_gap==0 ? 10000 : best_score-best_score_gap );
				if ( max_score<minlen ) break;
				for ( Size choose_tie=1; ;++choose_tie ) {
					Size num_ties;
					ungapped_align_degnucs_to_nucs( ndn_nucseq, d_seq_inc_pnucs, minlen, max_score, al, num_ties,
						-1000,                                         // mismatch_score: no mismatches
						max_vdj_insert_,                               // max_i_start
						int(ndn_nucseq.size()) - 1 - max_vdj_insert_,  // min_i_stop
						choose_tie );                                  // which tie to choose
					if ( best_score_gap==0 && choose_tie==1 ) best_score = al.score; // set this the first time through
					if ( choose_tie > num_ties || al.score < minlen ) break;
					// how many of the actual d gene's nucleotides are aligned?
					int const true_d_alignlen( min( al.j_stop(), d_gene_length+max_pnucs_ ) - max( al.j_start(), max_pnucs_ )+ 1);
					if ( true_d_alignlen >= minlen ) {
						vector< vector< Real > > const & trim_probs( d_trim_probs_.find(d_gene)->second );
						int const d_gene_length_inc_pnucs( d_seq_inc_pnucs.size() ),
							//len( al.length() ),
							max_total_extra_trim( true_d_alignlen - minlen ), // still want at least min_d_trimmed_length
							base_d0_trim( al.j_start() ),
							base_d1_trim( d_gene_length_inc_pnucs - 1 - al.j_stop() ),
							base_vd_insert( al.i_start ),
							base_dj_insert( ndn_len -1 - al.i_stop );
						for ( int extra_d0_trim=0; extra_d0_trim<= max_extra_d_trim-best_score_gap; ++extra_d0_trim ) {
							for ( int extra_d1_trim=0; extra_d1_trim <= max_extra_d_trim-best_score_gap &&
											extra_d0_trim+extra_d1_trim<=max_total_extra_trim; ++extra_d1_trim ) {
								if ( base_vd_insert + extra_d0_trim <= max_vdj_insert_ &&
									base_dj_insert + extra_d1_trim <= max_vdj_insert_ ) {
									vd_insertseq = ndn_nucseq.substr( 0, base_vd_insert+extra_d0_trim );
									jd_insertseq = reverse_complement( ndn_nucseq.substr( ndn_len - (extra_d1_trim+base_dj_insert) ) );
									Real prob_v(0); // this includes the V gene prob, the V trim prob, and the vd_insertseq nuc-prob
									for ( Size v_prefix=0; v_prefix< sz_alpha; ++v_prefix ) {
										if ( v_probs_by_prefix[v_prefix]==0 ) continue;
										prob_v += v_probs_by_prefix[v_prefix] * nucseq_probability( vd_insertseq, alphabet_[ v_prefix ] );
									}
									Real prob_j(0); // this includes the J gene prob, the J trim prob, and the jd_insertseq nuc-prob
									for ( Size j_prefix=0; j_prefix< sz_alpha; ++j_prefix ) {
										if ( j_probs_by_prefix[j_prefix]==0 ) continue;
										prob_j += j_probs_by_prefix[j_prefix] * nucseq_probability( jd_insertseq, alphabet_[ j_prefix ] );
									}

									Real const prob_with_these_trims
										( d_success_rate_ * prob_d_gene * prob_v * prob_j *
											vd_insert_probs_[ base_vd_insert+extra_d0_trim ] *
											trim_probs[ base_d0_trim+extra_d0_trim ][ base_d1_trim+extra_d1_trim ] *
											dj_insert_probs_[ base_dj_insert + extra_d1_trim ] );
									prob_with_this_d += prob_with_these_trims;
									if ( verbose() ) {
										cout << "trims_prob: " << prob_with_these_trims <<
											" align_len: " << al.length() <<
											" base_trims: " << base_d0_trim << ' ' << base_d1_trim <<
											" base_inserts: " << base_vd_insert << ' ' << base_dj_insert <<
											" extra_trimes: " << extra_d0_trim << ' ' << extra_d1_trim <<
											" prob_d: " << prob_d_gene <<
											" num_ties: " << num_ties <<
											" choose_tie: " << choose_tie <<
											" bsg: " << best_score_gap << endl;
									}
								}
							}
						}
					}
				} // choose_tie
			} // best_score_gap for this d-gene alignments
			prob += prob_with_this_d;
			if ( verbose() ) {
				cout << "prob_with_this_d: " << prob_with_this_d << ' ' << d_gene << ' ' << prob_d_gene << endl;
			}
		}
	}

	return prob;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// here we are assuming that this cdr3_nucseq is "error-corrected", ie that mismatches deep in v/j regions
// have been mutated to agree with the v/j sequence. So when we figure the trims we are insisting on exact matches!
//
//
Real
JunctionCounts::calc_cdr3_nucseq_probability(
	string const & cdr3_nucseq,
	strings const & v_ties,
	strings const & j_ties
) const
{
	//int const MAX_CDR3_NUCSEQ_LENGTH( 100 );
	int const max_extra_trim( vdj_probs::max_extra_vj_trim_for_probs );
	Size const sz_alpha( alphabet_.size() );

	if ( verbose() ) {
		cout << "calc_cdr3_nucseq_probability: " << cdr3_nucseq << " v_ties: ";
		foreach_( string v, v_ties ) cout << ' ' << v;
		cout << " j_ties: ";
		foreach_( string j, j_ties ) cout << ' ' << j;
		cout << endl;
	}

	// scratch space for the calculation
	static vector< vector< Real > > v_probs_for_ndn_start, j_probs_for_ndn_stop;
	static vector< bool > j1_gene_for_ndn_stop, j2_gene_for_ndn_stop;
	if ( v_probs_for_ndn_start.size() < cdr3_nucseq.size() ||
		v_probs_for_ndn_start.front().size() != sz_alpha ) {
		// dimension the static arrays
		Size const sz_big( cdr3_nucseq.size()*2 ); // silly
		v_probs_for_ndn_start.resize( sz_big );
		j_probs_for_ndn_stop.resize( sz_big );
		j1_gene_for_ndn_stop.resize( sz_big );
		j2_gene_for_ndn_stop.resize( sz_big );
		for ( Size i=0; i< sz_big; ++i ) {
			v_probs_for_ndn_start[i].resize( sz_alpha );
			j_probs_for_ndn_stop[i].resize( sz_alpha );
		}
	}

	int const cdr3_nucseq_len( cdr3_nucseq.size() );

	if ( v_ties.empty() || j_ties.empty() ) {
		return 0.0;
	}

	string single_j1_gene_for_d_choice, single_j2_gene_for_d_choice;
	foreach_( string const & g, j_ties ) {
		if      ( g[4] == '1' ) single_j1_gene_for_d_choice = g;
		else if ( g[4] == '2' ) single_j2_gene_for_d_choice = g;
		else {
			cerr << "funny j gene " << g << endl;
			exit(1);
		}
	}

	// whats the minimum that each v/j gene could have been trimmed back?
	vector<int> v_min_trims, v_ndn_starts;

	foreach_( string const & g, v_ties ) {
		string const & g_cdr3_nucseq( v_fasta_cdr3_region_inc_pnucs_.find(g)->second ); // includes pnucs
		int num_aligned(0);
		while ( num_aligned < (int) g_cdr3_nucseq.size() && num_aligned < cdr3_nucseq_len &&
			degnuc_matches_nuc( cdr3_nucseq[num_aligned], g_cdr3_nucseq[num_aligned] ) ) ++num_aligned;
		v_min_trims.push_back( g_cdr3_nucseq.size() - num_aligned );
		v_ndn_starts.push_back( num_aligned );
		if ( verbose() && v_min_trims.back()<max_pnucs_ ) {
			cout << "v_pnucs: ndn_start= " << v_ndn_starts.back() << ' ' << max_pnucs_-v_min_trims.back() << ' ' <<
				g << endl;
		}
	}

	vector<int> j_min_trims, j_ndn_stops;

	foreach_( string const & g, j_ties ) {
		string const & g_cdr3_nucseq( j_fasta_cdr3_region_inc_pnucs_.find(g)->second );
		int num_aligned(0); // Not allowing num_aligned to equal cdr3_nucseq_len, causes ndn_stop=-1 problems
		while ( num_aligned < (int) g_cdr3_nucseq.size() && num_aligned < cdr3_nucseq_len-1 &&
			degnuc_matches_nuc( cdr3_nucseq[cdr3_nucseq_len-1-num_aligned],
				g_cdr3_nucseq[g_cdr3_nucseq.size()-1-num_aligned] ) ) ++num_aligned;
		j_min_trims.push_back( g_cdr3_nucseq.size() - num_aligned );
		j_ndn_stops.push_back( cdr3_nucseq_len-1 - num_aligned );
		if ( verbose() && j_min_trims.back()<max_pnucs_ ) {
			cout << "j_pnucs: ndn_stop= " << j_ndn_stops.back() << ' ' << max_pnucs_-j_min_trims.back() << ' ' <<
				g << endl;
		}
	}


	int const max_ndn_start( max( v_ndn_starts ) ), min_ndn_start( min( v_ndn_starts ) );

	int v_trim(0), j_trim(0), extra_trim(0);
	for ( int ndn_start=max(0,min_ndn_start-max_extra_trim); ndn_start<= max_ndn_start; ++ndn_start ) {
		// which v-genes are compatible with this ndn_start?
		vector<Real> & v_probs_for_this_ndn_start( v_probs_for_ndn_start[ndn_start] );
		for ( Size i=0; i<sz_alpha; ++i ) v_probs_for_this_ndn_start[i] = 0.0;
		for ( Size gi=0; gi< v_ties.size(); ++gi ) {
			string const & g( v_ties[gi] );
			extra_trim = v_ndn_starts[gi] - ndn_start;
			if ( extra_trim <= max_extra_trim && extra_trim >= 0 ) {
				v_trim = v_min_trims[gi] + extra_trim;
				if ( v_trim <= max_trim_ + max_pnucs_ ) {
					string const & g_cdr3_nucseq_wpnucs( v_fasta_cdr3_region_inc_pnucs_.find(g)->second ); // includes pnucs
					string const & g_full_nucseq( v_fasta_.find(g)->second ); // no pnucs
					char const v_prefix( v_trim < max_pnucs_ ?
						g_cdr3_nucseq_wpnucs[ g_cdr3_nucseq_wpnucs.size()-1-v_trim ] :
						g_full_nucseq[ g_full_nucseq.size() - 1 - (v_trim-max_pnucs_) ] );
					Size const v_prefix_pos( alphabet_.find(v_prefix) );
					runtime_assert( v_prefix_pos != string::npos );
					v_probs_for_this_ndn_start[ v_prefix_pos ] +=
						v_probs_.find(g)->second * v_trim_probs_.find(g)->second[ v_trim ];
				}
			}
		}
		if ( verbose() ) {
			cout << "v_probs_for_ndn_start: " << ndn_start << ' ' <<
				v_probs_for_this_ndn_start[0] << ' ' <<
				v_probs_for_this_ndn_start[1] << ' ' <<
				v_probs_for_this_ndn_start[2] << ' ' <<
				v_probs_for_this_ndn_start[3] << endl;
		}
	}

	int const max_ndn_stop( max( j_ndn_stops ) ), min_ndn_stop( min( j_ndn_stops ) );

	for ( int ndn_stop=min_ndn_stop; ndn_stop <= max_ndn_stop+max_extra_trim && ndn_stop<cdr3_nucseq_len; ++ndn_stop ) {
		// which j-genes are compatible with this ndn_stop?
		runtime_assert( ndn_stop>=0 && ndn_stop < (int) j_probs_for_ndn_stop.size() );
		vector<Real> & j_probs_for_this_ndn_stop( j_probs_for_ndn_stop[ ndn_stop ] );
		for ( Size i=0; i<sz_alpha; ++i ) j_probs_for_this_ndn_stop[i] = 0.0;
		j1_gene_for_ndn_stop[ndn_stop] = false;
		j2_gene_for_ndn_stop[ndn_stop] = false;
		for ( Size gi=0; gi< j_ties.size(); ++gi ) {
			string const & g( j_ties[gi ] );
			extra_trim = ndn_stop - j_ndn_stops[gi];
			if ( extra_trim <= max_extra_trim && extra_trim >= 0 ) {
				j_trim = j_min_trims[gi] + extra_trim;
				if ( j_trim <= max_trim_+max_pnucs_ ) {
					string const & g_cdr3_nucseq_wpnucs( j_fasta_cdr3_region_inc_pnucs_.find(g)->second ); // includes pnucs
					string const & g_full_nucseq( j_fasta_.find(g)->second ); // no pnucs
					char const j_prefix( j_trim < max_pnucs_ ?
						reverse_complement_nucleotide( g_cdr3_nucseq_wpnucs[ j_trim ] ) :
						reverse_complement_nucleotide( g_full_nucseq[ j_trim-max_pnucs_ ] ) );
					Size const j_prefix_pos( alphabet_.find(j_prefix) );
					runtime_assert( j_prefix_pos != string::npos );
					j_probs_for_this_ndn_stop[ j_prefix_pos ] +=
						j_probs_.find(g)->second * j_trim_probs_.find(g)->second[ j_trim ];
					if      ( g[4] == '1' ) j1_gene_for_ndn_stop[ ndn_stop ] = true;
					else if ( g[4] == '2' ) j2_gene_for_ndn_stop[ ndn_stop ] = true;
				}
			}
		}
		if ( verbose() ) {
			cout << "j_probs_for_ndn_stop: " << ndn_stop << ' ' <<
				j_probs_for_this_ndn_stop[0] << ' ' <<
				j_probs_for_this_ndn_stop[1] << ' ' <<
				j_probs_for_this_ndn_stop[2] << ' ' <<
				j_probs_for_this_ndn_stop[3] << endl;
		}
	}

	Real prob(0);
	for ( int ndn_start=max(0,min_ndn_start-max_extra_trim); ndn_start<= max_ndn_start; ++ndn_start ) {
		for ( int ndn_stop=min_ndn_stop; ndn_stop <= max_ndn_stop+max_extra_trim && ndn_stop<cdr3_nucseq_len; ++ndn_stop){
			if ( ndn_start > ndn_stop+1 ) continue; // impossible
			if ( !( j1_gene_for_ndn_stop[ ndn_stop ] || j2_gene_for_ndn_stop[ ndn_stop ] ) ) continue;
			Real this_ndn_prob(0);
			if ( ndn_start == ndn_stop+1 ) {
				// 0 inserted nucleotides
				Real v_prob(0), j_prob(0);
				for ( Size i=0; i<sz_alpha; ++i ) {
					v_prob += v_probs_for_ndn_start[ ndn_start ][i];
					j_prob += j_probs_for_ndn_stop [ ndn_stop  ][i];
				}
				// no chance for a d gene in there
				this_ndn_prob = v_prob * j_prob * (1-d_success_rate_) * vj_insert_probs_[0];
			} else {
				string const & single_j_gene_for_d_choice( j2_gene_for_ndn_stop[ ndn_stop ] ?
					single_j2_gene_for_d_choice : single_j1_gene_for_d_choice );
				runtime_assert( !single_j_gene_for_d_choice.empty() );
				string const ndn_nucseq( cdr3_nucseq.substr( ndn_start, ndn_stop-ndn_start+1 ) );
				this_ndn_prob = calc_ndn_nucseq_probability_given_j_gene( ndn_nucseq, single_j_gene_for_d_choice,
					v_probs_for_ndn_start[ ndn_start ], j_probs_for_ndn_stop [ ndn_stop ] );
			}
			if ( verbose() ) {
				cout << "this_ndn_prob: " << ndn_start << ' ' << ndn_stop << " ndn_len: " << ndn_stop-ndn_start+1 <<
					' ' << this_ndn_prob << endl;
			}
			prob += this_ndn_prob;
		}
	}

	return prob;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Real
JunctionCounts::calc_cdr3_protseq_probability(
	string const & cdr3_protseq,
	strings const & v_ties,
	strings const & j_ties
) const
{
	if ( v_ties.empty() || j_ties.empty() ) return 0.0;

	strings cdr3_nucseqs;
	get_degenerate_coding_sequences( cdr3_protseq, cdr3_nucseqs );

	Real prob(0), single_prob(0);
	foreach_( string const & cdr3_nucseq, cdr3_nucseqs ) {
		single_prob = calc_cdr3_nucseq_probability( cdr3_nucseq, v_ties, j_ties );
		prob += single_prob;
		if ( verbose() ) {
			cout << "degenerate_nucseq: " << cdr3_protseq << ' ' << cdr3_nucseq << ' ' << single_prob << ' ' <<
				prob << endl;
		}
	}

	return prob;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
