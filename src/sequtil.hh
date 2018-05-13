// Simple utility functions involving nucleotide sequences
//

#ifndef INCLUDED_sequtil_hh
#define INCLUDED_sequtil_hh


#include "types.hh"
// #include "misc.hh"
// #include "io.hh"


void
get_degenerate_coding_sequences(
	string const & protseq,
	strings & nucseqs
)
{
	static map<char,strings> aa2degenerate_codons;

	if ( aa2degenerate_codons.empty() ) {
		aa2degenerate_codons['A'] = make_vector( string("gcn") );
		aa2degenerate_codons['C'] = make_vector( string("tgy") );
		aa2degenerate_codons['E'] = make_vector( string("gar") );
		aa2degenerate_codons['D'] = make_vector( string("gay") );
		aa2degenerate_codons['G'] = make_vector( string("ggn") );
		aa2degenerate_codons['F'] = make_vector( string("tty") );
		aa2degenerate_codons['I'] = make_vector( string("ath") );
		aa2degenerate_codons['H'] = make_vector( string("cay") );
		aa2degenerate_codons['K'] = make_vector( string("aar") );
		aa2degenerate_codons['M'] = make_vector( string("atg") );
		aa2degenerate_codons['L'] = make_vector( string("ttr") , string("ctn") );
		aa2degenerate_codons['N'] = make_vector( string("aay") );
		aa2degenerate_codons['Q'] = make_vector( string("car") );
		aa2degenerate_codons['P'] = make_vector( string("ccn") );
		aa2degenerate_codons['S'] = make_vector( string("tcn") , string("agy") );
		aa2degenerate_codons['R'] = make_vector( string("cgn") , string("agr") );
		aa2degenerate_codons['T'] = make_vector( string("acn") );
		aa2degenerate_codons['W'] = make_vector( string("tgg") );
		aa2degenerate_codons['V'] = make_vector( string("gtn") );
		aa2degenerate_codons['Y'] = make_vector( string("tay") );
		aa2degenerate_codons['*'] = make_vector( string("tar") , string("tga") ); // stop codon is '*'
	}

	nucseqs.clear();
	nucseqs.push_back( string() ); // empty string

	for ( Size i=0; i< protseq.size(); ++i ) {
		strings old_nucseqs;
		old_nucseqs.swap( nucseqs );
		runtime_assert( nucseqs.empty() );
		runtime_assert( aa2degenerate_codons.count(protseq[i] ) );
		strings const & codons( aa2degenerate_codons.find( protseq[i] )->second );
		foreach_( string const & old_nucseq, old_nucseqs ) {
			foreach_( string const & codon, codons ) {
				nucseqs.push_back( old_nucseq + codon );
			}
		}
	}
}



string
get_translation( string const & nucseq )
{
	static map<string,char> genetic_code;
	if ( genetic_code.empty() ) {
		genetic_code["aaa"] = 'K';
		genetic_code["aac"] = 'N';
		genetic_code["aag"] = 'K';
		genetic_code["aat"] = 'N';
		genetic_code["aca"] = 'T';
		genetic_code["acc"] = 'T';
		genetic_code["acg"] = 'T';
		genetic_code["act"] = 'T';
		genetic_code["aga"] = 'R';
		genetic_code["agc"] = 'S';
		genetic_code["agg"] = 'R';
		genetic_code["agt"] = 'S';
		genetic_code["ata"] = 'I';
		genetic_code["atc"] = 'I';
		genetic_code["atg"] = 'M';
		genetic_code["att"] = 'I';
		genetic_code["caa"] = 'Q';
		genetic_code["cac"] = 'H';
		genetic_code["cag"] = 'Q';
		genetic_code["cat"] = 'H';
		genetic_code["cca"] = 'P';
		genetic_code["ccc"] = 'P';
		genetic_code["ccg"] = 'P';
		genetic_code["cct"] = 'P';
		genetic_code["cga"] = 'R';
		genetic_code["cgc"] = 'R';
		genetic_code["cgg"] = 'R';
		genetic_code["cgt"] = 'R';
		genetic_code["cta"] = 'L';
		genetic_code["ctc"] = 'L';
		genetic_code["ctg"] = 'L';
		genetic_code["ctt"] = 'L';
		genetic_code["gaa"] = 'E';
		genetic_code["gac"] = 'D';
		genetic_code["gag"] = 'E';
		genetic_code["gat"] = 'D';
		genetic_code["gca"] = 'A';
		genetic_code["gcc"] = 'A';
		genetic_code["gcg"] = 'A';
		genetic_code["gct"] = 'A';
		genetic_code["gga"] = 'G';
		genetic_code["ggc"] = 'G';
		genetic_code["ggg"] = 'G';
		genetic_code["ggt"] = 'G';
		genetic_code["gta"] = 'V';
		genetic_code["gtc"] = 'V';
		genetic_code["gtg"] = 'V';
		genetic_code["gtt"] = 'V';
		genetic_code["taa"] = '*';
		genetic_code["tac"] = 'Y';
		genetic_code["tag"] = '*';
		genetic_code["tat"] = 'Y';
		genetic_code["tca"] = 'S';
		genetic_code["tcc"] = 'S';
		genetic_code["tcg"] = 'S';
		genetic_code["tct"] = 'S';
		genetic_code["tga"] = '*';
		genetic_code["tgc"] = 'C';
		genetic_code["tgg"] = 'W';
		genetic_code["tgt"] = 'C';
		genetic_code["tta"] = 'L';
		genetic_code["ttc"] = 'F';
		genetic_code["ttg"] = 'L';
		genetic_code["ttt"] = 'F';
	}

	runtime_assert( nucseq.size()%3==0 );
	Size const num_codons( nucseq.size()/3 );

	string protseq;
	for ( Size i=0; i<num_codons; ++i ) {
		protseq.push_back( genetic_code.find( nucseq.substr(3*i,3) )->second );
	}
	return protseq;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// here we consider k to match y, for example, since they both contain t... and h matches everything except 'g'
//
bool
degnuc_matches_degnuc( char const a, char const b )
{
	switch( a ) {
	case 'a':
		return (  b!='c'  &&  b!='b'  &&  b!='g'  &&  b!='k'  &&  b!='s'  &&  b!='t'  &&  b!='y'  );
	case 'c':
		return (  b!='a'  &&  b!='d'  &&  b!='g'  &&  b!='k'  &&  b!='r'  &&  b!='t'  &&  b!='w'  );
	case 'b':
		return (  b!='a'  );
	case 'd':
		return (  b!='c'  );
	case 'g':
		return (  b!='a'  &&  b!='c'  &&  b!='h'  &&  b!='m'  &&  b!='t'  &&  b!='w'  &&  b!='y'  );
	case 'h':
		return (  b!='g'  );
	case 'k':
		return (  b!='a'  &&  b!='c'  &&  b!='m'  );
	case 'm':
		return (  b!='g'  &&  b!='k'  &&  b!='t'  );
	case 'n':
		return true;
	case 's':
		return (  b!='a'  &&  b!='t'  &&  b!='w'  );
	case 'r':
		return (  b!='c'  &&  b!='t'  &&  b!='y'  );
	case 't':
		return (  b!='a'  &&  b!='c'  &&  b!='g'  &&  b!='m'  &&  b!='s'  &&  b!='r'  &&  b!='v'  );
	case 'w':
		return (  b!='c'  &&  b!='g'  &&  b!='s'  );
	case 'v':
		return (  b!='t'  );
	case 'y':
		return (  b!='a'  &&  b!='g'  &&  b!='r'  );
	default:
		cout << "ERROR degnuc_matches_degnuc: unrecognized anuc: " << a << ' ' << b << endl;
		cerr << "ERROR degnuc_matches_degnuc: unrecognized anuc: " << a << ' ' << b << endl;
		exit(1);
		return false;
	}
	return false;

}

bool
degnuc_matches_nuc( char const a, char const b )
{
	static bool const n_matches_everything( true ); //!option[ my_options::n_matches_nothing ] );
	switch ( a ) {
	case 'a':
		return (  b=='a'  );
	case 'c':
		return (  b=='c'  );
	case 'g':
		return (  b=='g'  );
	case 't':
		return (  b=='t'  );
	case 'n':
		return n_matches_everything;
	case 'b':
		return (  b!='a' );
	case 'd':
		return (  b!='c' );
	case 'h':
		return (  b!='g' );
	case 'k':
		return (  b=='g'  ||  b=='t'  );
	case 'm':
		return (  b=='a'  ||  b=='c'  );
	case 's':
		return (  b=='c'  ||  b=='g'  );
	case 'r':
		return (  b=='a'  ||  b=='g'  );
	case 'w':
		return (  b=='a'  ||  b=='t'  );
	case 'v':
		return (  b!='t' );
	case 'y':
		return (  b=='c'  ||  b=='t'  );
	default:
		cout << "ERROR degnuc_matches_nuc: unrecognized degnuc: " << a << ' ' << b << endl;
		cerr << "ERROR degnuc_matches_nuc: unrecognized degnuc: " << a << ' ' << b << endl;
		exit(1);
		return false;
	}
	return false;

}

struct UngappedAlignment {
	UngappedAlignment():
		score(0),offset(0),i_start(0),i_stop(0),mismatches(0),reversed(false){}

	void
	clear() {score = offset = i_start = i_stop = 0; mismatches = 0; reversed=false; }

	// these assume the definition of offset used in ungapped_align below...
	int
	j_start() const { return i_start + offset; }

	int
	j_stop() const { return i_stop + offset; }

	int
	length() const { return i_stop - i_start + 1 ; }

	int score;
	int offset;
	int i_start;
	int i_stop;
	int mismatches;
	bool reversed; // is the query sequence reverse-complemented?
};




void
ungapped_align(
	string const & a,
	string const & b,
	int const min_score,
	UngappedAlignment & best_align,
	int const mismatch_score = -3 // from: https://www.arabidopsis.org/Blast/BLASToptions.jsp
)
{
	// Allow 'x' in first (non-db) sequence; it matches nothing
	foreach_( char const bb, a ) { runtime_assert( bb=='a'|| bb=='c' || bb=='g' || bb=='t' || bb=='x' ); }
	foreach_( char const bb, b ) { runtime_assert( bb=='a'|| bb=='c' || bb=='g' || bb=='t' ); }

	int const alen( a.size() );
	int const blen( b.size() );

	// offset is how much a is shifted relative to b
	// ie a[0] aligns with b[offset]
	//    a[alen-1] aligns with b[alen-1+offset]
	//
	// to have any overlap, we need offset <= blen-1
	// and alen-1+offset >=0  ie  offset >= 1-alen
	//
	int const match_score(1);

	best_align.clear();
	best_align.score = min_score-1;
	best_align.reversed = false;
	//vector<int> scores( alen );

	for ( int offset = 1-alen; offset <= blen-1; ++offset ) {
		int const i_start( max(0,-1*offset) ), i_stop( min(alen-1,blen-1-offset) );
		int align_begin(-1), align_score(0);
		for ( int i=i_start; i<= i_stop; ++i ) {
			if ( a[i] == b[i+offset] ) { // match
				if ( align_score==0 ) align_begin = i; // not currently aligning
				align_score += match_score;
			} else { // mismatch
				// were we at the best alignment?
				if ( align_score > best_align.score ) {
					best_align.score = align_score;
					best_align.offset = offset;
					best_align.i_start = align_begin;
					best_align.i_stop = i-1;
				}
				align_score += mismatch_score;
				if ( align_score < 0 ) align_score = 0; // stop aligning
			}
		} // i=i_start,i_stop
		if ( align_score > best_align.score ) { // aligned all the way to the end
			best_align.score = align_score;
			best_align.offset = offset;
			best_align.i_start = align_begin;
			best_align.i_stop = i_stop;
		}
	} // offset

	//count mismatches, sanity check on score, alignment
	if ( best_align.score >= min_score ) {
		best_align.mismatches = 0;
		int rescore(0), istart( best_align.i_start ), istop( best_align.i_stop ),
			jstart( istart + best_align.offset ), jstop( istop + best_align.offset  );
		// best_align should be flanked by mismatches, or the sequence ends
		runtime_assert( istart==0 || jstart==0 || a[istart-1] != b[jstart-1] );
		runtime_assert( istop==(int)a.size()-1 || jstop==(int)b.size()-1 || a[istop+1] != b[jstop+1] );
		for ( int i=istart; i<= istop; ++i ) {
			if ( a[i] != b[i+best_align.offset] ) {
				++best_align.mismatches;
				rescore += mismatch_score;
			} else {
				rescore += match_score;
			}
		}
		runtime_assert( rescore == best_align.score );
	}

}

void
ungapped_align_degnucs_to_nucs(
	string const & a, // nucleotide sequence possibly containing degenerate nucleotides
	string const & b, // sequence of acgt's
	int const min_score,
	int const max_score,
	UngappedAlignment & best_align,
	Size & num_ties,
	int const mismatch_score = -3, // from: https://www.arabidopsis.org/Blast/BLASToptions.jsp
	int max_i_start = -1,
	int const min_i_stop = -1,
	Size const choose_tie = 1 // pick the first by default
)
{
	// the degnuc_matches_nuc function checks the degnucs in a, so just confirm b sequence
	foreach_( char const bb, b ) { runtime_assert( bb=='a'|| bb=='c' || bb=='g' || bb=='t' ); }

	int const alen( a.size() );
	int const blen( b.size() );

	if ( max_i_start == -1 ) max_i_start = alen; // signal that default value was used

	// offset is how much a is shifted relative to b
	// ie a[0] aligns with b[offset]
	//    a[alen-1] aligns with b[alen-1+offset]
	//
	// to have any overlap, we need offset <= blen-1
	// and alen-1+offset >=0  ie  offset >= 1-alen
	//
	int const match_score(1);

	best_align.clear();
	best_align.score = min_score-1;
	best_align.reversed = false;
	//vector<int> scores( alen );

	num_ties = 0;

	for ( int offset = 1-alen; offset <= blen-1; ++offset ) {
		int const i_start( max(0,-1*offset) ), i_stop( min(alen-1,blen-1-offset) );
		int align_begin(-1), align_score(0);
		for ( int i=i_start; i<= i_stop; ++i ) {
			if ( degnuc_matches_nuc( a[i], b[i+offset] ) ) { // match
				if ( align_score==0 ) align_begin = i; // not currently aligning
				align_score += match_score;
			} else { // mismatch
				// were we at the best alignment?
				if ( align_score >= best_align.score && align_score <= max_score && align_begin <= max_i_start &&
					i-1 >= min_i_stop ) {
					if ( align_score > best_align.score ) {
						num_ties=1;
					} else {
						++num_ties; // 2 or more
					}
					if ( num_ties==1 || num_ties == choose_tie ) { // always want to record the best score we've seen...
						best_align.score = align_score;
						best_align.offset = offset;
						best_align.i_start = align_begin;
						best_align.i_stop = i-1;
					}
				}
				align_score += mismatch_score;
				if ( align_score < 0 ) align_score = 0; // stop aligning
			}
		} // i=i_start,i_stop
		if ( align_score > best_align.score ) { // aligned all the way to the end
			best_align.score = align_score;
			best_align.offset = offset;
			best_align.i_start = align_begin;
			best_align.i_stop = i_stop;
		}
	} // offset

	//count mismatches, sanity check on score, alignment
	if ( best_align.score >= min_score ) {
		best_align.mismatches = 0;
		int rescore(0), istart( best_align.i_start ), istop( best_align.i_stop ),
			jstart( istart + best_align.offset ), jstop( istop + best_align.offset  );
		// best_align should be flanked by mismatches, or the sequence ends
		runtime_assert( istart==0 || jstart==0 || a[istart-1] != b[jstart-1] );
		runtime_assert( istop==(int)a.size()-1 || jstop==(int)b.size()-1 || a[istop+1] != b[jstop+1] );
		for ( int i=istart; i<= istop; ++i ) {
			if ( !degnuc_matches_nuc( a[i], b[i+best_align.offset] ) ) {
				++best_align.mismatches;
				rescore += mismatch_score;
			} else {
				rescore += match_score;
			}
		}
		runtime_assert( rescore == best_align.score );
	}

}

/// another version with an interface that parallels the original ungapped_align function
void
ungapped_align_degnucs_to_nucs(
	string const & a,
	string const & b,
	int const min_score,
	UngappedAlignment & best_align,
	int const mismatch_score = -3 // from: https://www.arabidopsis.org/Blast/BLASToptions.jsp
)
{
	Size num_ties; // unused
	Size const max_score( 10*a.size() ); //want this to be really big; relying on fact that match_score is 1, above
	ungapped_align_degnucs_to_nucs( a, b, min_score, max_score, best_align, num_ties, mismatch_score );
}



inline
char
reverse_complement_nucleotide( char const n )
{
	switch ( n ) {
	case 'a':
		return 't'; // a = a, t = t
	case 'c':
		return 'g'; // c = c, g = g
	case 'b':
		return 'v'; // b = cgt, v = acg
	case 'd':
		return 'h'; // d = agt, h = act
	case 'g':
		return 'c'; // g = g, c = c
	case 'h':
		return 'd'; // h = act, d = agt
	case 'k':
		return 'm'; // k = gt, m = ac
	case 'm':
		return 'k'; // m = ac, k = gt
	case 'n':
		return 'n'; // n = acgt, n = acgt
	case 's':
		return 's'; // s = cg, s = cg
	case 'r':
		return 'y'; // r = ag, y = ct
	case 't':
		return 'a'; // t = t, a = a
	case 'w':
		return 'w'; // w = at, w = at
	case 'v':
		return 'b'; // v = acg, b = cgt
	case 'y':
		return 'r'; // y = ct, r = ag
	case 'x':
		return 'x';
	default:
		runtime_assert( false );
		return 'x';
	}
	return 'x';
}

string
reverse_complement( string const & seq )
{
	string revseq;
	foreach_ ( char const c, seq ) {
		revseq.push_back( reverse_complement_nucleotide( c ) );
	}
	std::reverse( revseq.begin(), revseq.end() );
	return revseq;
}




Real
nucleotide_probability( char const x )
{
	switch ( x ) {

	case 'a':
		return 0.25;
	case 'c':
		return 0.25;
	case 'b':
		return 0.75;
	case 'd':
		return 0.75;
	case 'g':
		return 0.25;
	case 'h':
		return 0.75;
	case 'k':
		return 0.50;
	case 'm':
		return 0.50;
	case 'n':
		return 1.00;
	case 's':
		return 0.50;
	case 'r':
		return 0.50;
	case 't':
		return 0.25;
	case 'w':
		return 0.50;
	case 'v':
		return 0.75;
	case 'y':
		return 0.50;
	default:
		cout << "ERROR nucleotide_probability:: unrecognized nuc: " << x << endl;
		cerr << "ERROR nucleotide_probability:: unrecognized nuc: " << x << endl;
		exit(1);
		return 0.0;
	}
	return 0.0;
}



Real
nucleotide_probability_pwm(
	char const x,
	vector<Real> const & probs // 0=a, 1=c, 2=g, 3=t
)
{
	switch ( x ) {

	case 'a':
		return probs[0];
	case 'c':
		return probs[1];
	case 'b':
		return probs[1] + probs[2] + probs[3];
	case 'd':
		return probs[0] + probs[2] + probs[3];
	case 'g':
		return probs[2];
	case 'h':
		return probs[0] + probs[1] + probs[3];
	case 'k':
		return probs[2] + probs[3];
	case 'm':
		return probs[0] + probs[1];
	case 'n':
		return probs[0] + probs[1] + probs[2] + probs[3];
	case 's':
		return probs[1] + probs[2];
	case 'r':
		return probs[0] + probs[2];
	case 't':
		return probs[3];
	case 'w':
		return probs[0] + probs[3];
	case 'v':
		return probs[0] + probs[1] + probs[2];
	case 'y':
		return probs[1] + probs[3];
	default:
		cout << "ERROR nucleotide_probability_pwm:: unrecognized nuc: " << x << endl;
		cerr << "ERROR nucleotide_probability_pwm:: unrecognized nuc: " << x << endl;
		exit(1);
		return 0.0;
	}
	return 0.0;
}


Real
simple_nucseq_probability( string const & nucseq )
{
	Real prob(1.0);
	foreach_( char c, nucseq ) prob *= nucleotide_probability(c);
	return prob;
}


#endif
