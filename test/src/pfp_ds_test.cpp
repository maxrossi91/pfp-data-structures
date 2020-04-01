/* pfp-ds - prefix free parsing data structures
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file pfp_ds_test.cpp
   \brief pfp_ds_test.cpp define build and test prefix-free parsing datastructures.
   \author Massimiliano Rossi
   \date 20/03/2020
*/

#include<iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
#include<gsacak.h>
}
#include<malloc_count.h>

//**************************** From  Big-BWT ***********************************
// special symbols used by the construction algorithm:
//   they cannot appear in the input file
//   the 0 symbol is used in the final BWT file as the EOF char

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter
//******************************************************************************

template<typename T, typename S, typename lcp_t>
void LCP_array(S* s, const std::vector<T>& isa, const std::vector<T>& sa, size_t n, std::vector<lcp_t>& lcp){
    lcp[0]  = 0;

    T l = 0;
    for (size_t i = 0; i < n; ++i){
      // if i is the last character LCP is not defined
      T k = isa[i];
      if(k > 0){
        T j = sa[k-1];
        // I find the longest common prefix of the i-th suffix and the j-th suffix.
        while(s[i+l] == s[j+l]) l++;
        // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
        lcp[k] = l;
        if(l>0) l--;
      }
    }
}

// TODO: Extend it to integer alphabets
class dictionary{
public:
  std::vector<uint8_t> d;
  std::vector<uint_t> saD;
  std::vector<uint_t> isaD;
  std::vector<int_t> daD;
  std::vector<int_t> lcpD;
  sdsl::rmq_succinct_sct<> rmq_lcp_D;
  sdsl::bit_vector b_d; // Starting position of each phrase in D
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;
  std::vector<int_t> colex_daD;
  sdsl::rmq_succinct_sct<> rmq_colex_daD;
  sdsl::range_maximum_sct<>::type rMq_colex_daD;
  bool saD_flag = false;
  bool isaD_flag = false;
  bool daD_flag = false;
  bool lcpD_flag = false;
  bool rmq_lcp_D_flag = false;


  dictionary( std::vector<uint8_t>& d_,
              bool saD_flag_ = true,
              bool isaD_flag_ = true,
              bool daD_flag_ = true,
              bool lcpD_flag_ = true,
              bool rmq_lcp_D_flag_ = true ):
              d(d_)
  {
    build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_);

  }

  dictionary( std::string filename,
              bool saD_flag_ = true,
              bool isaD_flag_ = true,
              bool daD_flag_ = true,
              bool lcpD_flag_ = true,
              bool rmq_lcp_D_flag_ = true )
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".dict");
    read_file(tmp_filename.c_str(), d);

    build(saD_flag_, isaD_flag_, daD_flag_, lcpD_flag_, rmq_lcp_D_flag_);

  }

  inline size_t length_of_phrase(size_t id){
    assert(id > 0);
    return select_b_d(id+1)-select_b_d(id) - 1; // to remove the EndOfWord
  }

  inline size_t n_phrases(){
    return rank_b_d(d.size()-1);
  }

  void build(bool saD_flag_, bool isaD_flag_, bool daD_flag_, bool lcpD_flag_, bool rmq_lcp_D_flag_){

    // Building the bitvector with a 1 in each starting position of each phrase in D
    b_d.resize(d.size());
    for(size_t i = 0; i < 3*64; ++i) b_d[i] = false; // bug in resize
    b_d[0] = true; // Mark the first phrase
    for(int i = 1; i < d.size(); ++i )
      b_d[i] = (d[i-1]==EndOfWord);
    b_d[d.size()-1] = true; // This is necessary to get the length of the last phrase

    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);


    assert(!rmq_lcp_D_flag_ || (lcpD_flag || lcpD_flag_));

    // TODO: check if it has been already computed
    if(saD_flag_ && daD_flag_ && lcpD_flag_){
      saD.resize(d.size());
      lcpD.resize(d.size());
      daD.resize(d.size());
      // suffix array, LCP array, and Document array of the dictionary.
      verbose("Computing SA, LCP, and DA of dictionary");
      elapsed_time(
        gsacak(&d[0],&saD[0],&lcpD[0],&daD[0],d.size())
      );
    }else if(saD_flag_ && lcpD_flag_){
      saD.resize(d.size());
      lcpD.resize(d.size());
      // suffix array and LCP array of the dictionary.
      verbose("Computing SA, and LCP of dictionary");
      elapsed_time(
        gsacak(&d[0],&saD[0],&lcpD[0],nullptr,d.size())
      );
    } else if(saD_flag_ && daD_flag_){
      saD.resize(d.size());
      daD.resize(d.size());
      // suffix array and LCP array of the dictionary.
      verbose("Computing SA, and DA of dictionary");
      elapsed_time(
        gsacak(&d[0],&saD[0],nullptr,&daD[0],d.size())
      );
    } else if(saD_flag_){
      saD.resize(d.size());
      // suffix array and LCP array of the dictionary.
      verbose("Computing SA of dictionary");
      elapsed_time(
        gsacak(&d[0],&saD[0],nullptr,nullptr,d.size())
      );
    }

    assert(!isaD_flag_ || (saD_flag || saD_flag_) );
    if(isaD_flag_ && !isaD_flag){
      // inverse suffix array of the dictionary.
      verbose("Computing ISA of dictionary");
      elapsed_time(
        {
          isaD.resize(d.size());
          for(int i = 0; i < saD.size(); ++i){
            isaD[saD[i]] = i;
          }
        }
      );
    }

    assert(!rmq_lcp_D_flag_ || (lcpD_flag || lcpD_flag_));
    if(rmq_lcp_D_flag_ && ! rmq_lcp_D_flag){
      rmq_lcp_D_flag = true;

      verbose("Computing RMQ over LCP of dictionary");
      // Compute the LCP rank of D
      elapsed_time(
        rmq_lcp_D = sdsl::rmq_succinct_sct<>(&lcpD)
      );
    }

    // if(colex_daD_flag_){
      // co-lex document array of the dictionary.
      verbose("Computing co-lex DA of dictionary");
      elapsed_time(
        {  
          std::vector<uint_t>colex_id(n_phrases());
          std::vector<uint_t>inv_colex_id(n_phrases()); // I am using it as starting positions
          for(int i = 0, j = 0; i < d.size(); ++i )
            if(d[i+1]==EndOfWord){
              colex_id[j] = j;
              inv_colex_id[j++] = i;
            }

          colex_document_array_helper(inv_colex_id,colex_id,0,n_phrases());

          // computing inverse colex id
          for(int i = 0; i < colex_id.size(); ++i){
            inv_colex_id[colex_id[i]] = i;
          }
          colex_id.clear();

          colex_daD.resize(d.size());
          for(int i = 0; i < colex_daD.size(); ++i ){
            colex_daD[i]  = inv_colex_id[daD[i]];
          }
        }

        rmq_colex_daD = sdsl::rmq_succinct_sct<>(&colex_daD);
        rMq_colex_daD = sdsl::range_maximum_sct<>::type(&colex_daD);
      );


    // }

  }

  // I am assuming that the alphabet fits in an unsigned char
  void colex_document_array_helper(std::vector<uint_t>& sp, std::vector<uint_t>& id, size_t start, size_t end ){
    if((end - start < 2) || (start > end)) return;


    std::vector<uint32_t> count(256,0);
    for(size_t i = start; i < end; ++i ){
      count[ d[ sp[i] ] ] ++;
    }

    std::vector<uint32_t> psum(256,0);
    for(size_t i = 1; i < 256; ++i){
      psum[i] = psum[i-1] + count[i-1];
    }

    std::vector<uint_t> tmp(end-start+1,0);
    std::vector<uint_t> tmp_id(end-start+1,0);
    for(size_t i = start; i < end; ++i ){
      size_t index = psum[ d[sp[i] ] ] ++;
      tmp[index] = std::min(sp[i]-1, static_cast<uint_t>(d.size()-1));
      tmp_id[index] = id[i];
    }


    // Recursion
    size_t tmp_start = 0;
    for(size_t i = 0; i < 256; ++i ){
      for(size_t j = 0; j < count[i]; ++j){
        sp[start + j] = tmp[tmp_start];
        id[start + j] = tmp_id[tmp_start ++];
      }
      end = start + count[i];
      if(i > EndOfWord) colex_document_array_helper(sp,id,start,end);
      start = end;
    }

  }

};


// TODO: Extend it to non-integer alphabets
class parse{
public:
  std::vector<uint32_t> p;
  std::vector<uint_t> saP;
  std::vector<uint_t> isaP;
  std::vector<int_t> lcpP;
  sdsl::rmq_succinct_sct<> rmq_lcp_P;
  sdsl::bit_vector b_p; // Starting position of each phrase in D
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;
  bool saP_flag = false;
  bool isaP_flag = false;
  bool lcpP_flag = false;
  bool rmq_lcp_P_flag = false;

  size_t alphabet_size;

  parse(  std::vector<uint32_t>& p_,
          size_t alphabet_size_,
          bool saP_flag_ = true,
          bool isaP_flag_ = true,
          bool lcpP_flag_ = true,
          bool rmq_lcp_P_flag_ = true ):
          p(p_),
          alphabet_size(alphabet_size_)
  {
    build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);


  }

  parse(  std::string filename,
          size_t alphabet_size_,
          bool saP_flag_ = true,
          bool isaP_flag_ = true,
          bool lcpP_flag_ = true,
          bool rmq_lcp_P_flag_ = true ):
          alphabet_size(alphabet_size_)
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".parse");
    read_file(tmp_filename.c_str(), p);
    p.push_back(0); // this is the terminator for the sacak algorithm

    build(saP_flag_, isaP_flag_, lcpP_flag_, rmq_lcp_P_flag_);

  }

  void build(bool saP_flag_, bool isaP_flag_, bool lcpP_flag_, bool rmq_lcp_P_flag_){

    // TODO: check if it has been already computed
    if(saP_flag_){
      saP.resize(p.size());
      // suffix array of the parsing.
      verbose("Computing SA of the parsing");
      elapsed_time(
        sacak_int(&p[0],&saP[0],p.size(),alphabet_size);
      );
    }

    assert(!isaP_flag_ || (saP_flag || saP_flag_) );
    if(isaP_flag_ && !isaP_flag){
      // inverse suffix array of the parsing.
      verbose("Computing ISA of the parsing");
      elapsed_time(
        {
          isaP.resize(p.size());
          for(int i = 0; i < saP.size(); ++i){
            isaP[saP[i]] = i;
          }
        }
      )

    }

    if(lcpP_flag_){
      lcpP.resize(p.size());
      // LCP array of the parsing.
      verbose("Computing LCP of the parsing");
      elapsed_time(
        LCP_array(&p[0], isaP, saP, p.size(), lcpP);
      );
    }


    assert(!rmq_lcp_P_flag_ || (lcpP_flag || lcpP_flag_));
    if(rmq_lcp_P_flag_ && ! rmq_lcp_P_flag){
      rmq_lcp_P_flag = true;
      verbose("Computing RMQ over LCP of the parsing");
      // Compute the LCP rank of P
      elapsed_time(
        rmq_lcp_P = sdsl::rmq_succinct_sct<>(&lcpP);
      );
    }

  }

};


class pfp_wt {
public:
  using wt_bv = sdsl::bit_vector;
  struct wt_node {
    wt_node (wt_node * parent = nullptr)
      : parent(parent)
    { }

    wt_bv bit_vector;
    wt_bv::rank_1_type bv_rank_1;
    wt_bv::select_1_type bv_select_1;

    uint32_t phrase_id = 0;

    wt_node * parent = nullptr;
    std::unique_ptr<wt_node> left;
    std::unique_ptr<wt_node> right;

    bool is_leaf () const {
      return !left.get();
    }

    bool is_root () const {
      return parent == nullptr;
    }
  };

  pfp_wt()
    : root(new wt_node())
  { }

  pfp_wt(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse)
    : root(new wt_node()) {
    create_bwt_rec(*root, sorted_alphabet, parse);
  }

  void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {
    create_bwt_rec(*root, sorted_alphabet, parse);
  }

  // public interface
  void print_leafs() {
    print_leafs_rec(root);
    std::cout << std::endl;
  }

  // operator[] - get phrase ID on i-th position
  uint32_t operator[] (const size_t i) {
    assert(i < root->bit_vector.size());

    return get_phrase_id(*root, i); // WT[] is 0-based
  }

private:
  // root node of wavelet tree
  std::unique_ptr<wt_node> root;

  uint32_t get_phrase_id(const wt_node & node, const size_t i) {
    if (node.is_leaf()) {
      return node.phrase_id;
    }

    const auto rank_bit_1 = node.bv_rank_1(i + 1);
    if (node.bit_vector[i]) {
      return get_phrase_id(*node.right, rank_bit_1 - 1);
    }
    else {
      return get_phrase_id(*node.left, i - rank_bit_1);
    }
  }

  void create_bwt_rec(wt_node & w_structure, const std::vector<uint32_t> & alphabet, const std::vector<uint32_t> & parse) {
    // divide alphabet to two sets
    const uint32_t tres = alphabet.size() / 2;
    if (tres == 0) {
      // leaf
      w_structure.phrase_id = alphabet[0];

      return;
    }

    // create bit_vector
    w_structure.bit_vector.resize(parse.size());
    std::vector<uint32_t> parse_left;
    std::vector<uint32_t> parse_right;

    for (size_t i = 0; i < parse.size(); i++) {
      // TODO: direct access data structure
      const auto iter = std::find(alphabet.begin(), alphabet.end(), parse[i]);
      const auto idx = std::distance(alphabet.begin(), iter);

      if (idx < tres) {
        w_structure.bit_vector[i] = 0;
        parse_left.push_back(parse[i]);
      }
      else {
        w_structure.bit_vector[i] = 1;
        parse_right.push_back(parse[i]);
      }
    }

    // rank & select support
    w_structure.bv_rank_1 = wt_bv::rank_1_type(&w_structure.bit_vector);
    w_structure.bv_select_1 = wt_bv::select_1_type(&w_structure.bit_vector);

    w_structure.left = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));
    w_structure.right = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));

    create_bwt_rec(*(w_structure.left), std::vector<uint32_t>(alphabet.begin(), alphabet.begin() + tres), parse_left);
    create_bwt_rec(*(w_structure.right), std::vector<uint32_t>(alphabet.begin() + tres, alphabet.end()), parse_right);
  }

  void print_leafs_rec(const std::unique_ptr<wt_node> & node) {
    if (node->is_leaf()) {
      std::cout << node->phrase_id << " | ";
    }
    else {
      print_leafs_rec(node->left);
      print_leafs_rec(node->right);
    }
  }
};

class pf_parsing{
public:
  struct M_entry_t{
    uint_t len;
    uint_t left; // left and right are the extreemes of the range
    uint_t right;
  };


  dictionary dict;
  parse pars;
  std::vector<uint32_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  sdsl::bit_vector b_bwt;
  std::vector<M_entry_t> M;

  pfp_wt w_wt;

  pf_parsing( std::string filename, size_t w_):
              dict(filename),
              pars(filename,dict.n_phrases()+1),
              freq(1,0),
              w(w_)
  {

    // Uploading the frequency file
    uint32_t *occ;
    size_t d_words;
    std::string tmp_filename = filename + std::string(".occ");
    read_file<uint32_t> (tmp_filename.c_str(), occ, d_words);
    freq.insert(freq.end(),occ, occ + d_words);


    // Compute the length of the string;
    n = 0;
    for(int j = 0; j < pars.p.size()-1; ++j){
      // parse.p[j]: phrase_id
      assert(pars.p[j] != 0);
      n += dict.length_of_phrase(pars.p[j]) - w;
    }
    n += w - 1; // -1 is for the first dollar + w because n is the length including the last w markers


    verbose("Computing b_bwt and M of the parsing");
    elapsed_time({
      // Build the bitvector storing the position of the beginning of each phrase.
      b_bwt.resize(n);
      for(size_t i = 0; i < 3*64; ++i) b_bwt[i] = false; // bug in resize

      assert(dict.d[dict.saD[0]] == EndOfDict);
      size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
      size_t j = 0;
      while (i < dict.saD.size()) {
        size_t left = j;

        auto  sn = dict.saD[i];
        // Check if the suffix has length at least w and is not the complete phrase.
        auto phrase = dict.daD[i] + 1; assert(phrase > 0 && phrase < freq.size()); // + 1 because daD is 0-based
        size_t suffix_length = dict.select_b_d(dict.rank_b_d(sn+1) +1) - sn - 1;
        if(dict.b_d[sn] || suffix_length < w){
          ++i; // Skip
        }else{
          // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
          // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
          b_bwt[j++] = true; j += freq[phrase]-1; // the next bits are 0s
          i++;
          if (i < dict.saD.size()){
            auto new_sn = dict.saD[i];
            auto new_phrase = dict.daD[i] + 1; assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
            size_t new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn+1) +1) - new_sn - 1;

            while (i < dict.saD.size() && (dict.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length)){
                j += freq[new_phrase];
                ++i;

                if (i < dict.saD.size()){
                  new_sn = dict.saD[i];
                  new_phrase = dict.daD[i] + 1; assert(new_phrase > 0 && new_phrase < freq.size()); // + 1 because daD is 0-based
                  new_suffix_length = dict.select_b_d(dict.rank_b_d(new_sn+1) +1) - new_sn - 1;
                }
            }
          }

          // Computing M
          size_t right = j-1;
          M_entry_t m;
          m.len = suffix_length;
          m.left = dict.colex_daD[dict.rmq_colex_daD(left,right)];
          m.right = dict.colex_daD[dict.rMq_colex_daD(left,right)];

          M.push_back(m);
        }
      }
    }); /// End elapsed_time

    // Add W wavelet tree
    verbose("Computing W of BWT(P)");
    // create alphabet (phrases)
    std::vector<uint32_t> alphabet(dict.n_phrases());
    std::iota(alphabet.begin(), alphabet.end(), 1);

    // TODO: use existing co-lex sorted phrases
    auto co_lexi_dict_cmp = [&](const uint32_t i, const uint32_t j) {
      auto i_start = dict.select_b_d(i);
      auto i_end = i_start + dict.length_of_phrase(i) - 1;
      auto j_start = dict.select_b_d(j);
      auto j_end = j_start + dict.length_of_phrase(j) - 1;

      auto i_r_begin = dict.d.rend() - i_end - 1;
      auto i_r_end = dict.d.rend() - i_start;
      auto j_r_begin = dict.d.rend() - j_end - 1;
      auto j_r_end = dict.d.rend() - j_start;

      return std::lexicographical_compare(i_r_begin, i_r_end, j_r_begin, j_r_end);
    };
    std::sort(alphabet.begin(), alphabet.end(), co_lexi_dict_cmp);

    // create BWT(P)
    std::vector<uint32_t> bwt_p(pars.p.size(), 0);
    for (size_t i = 0; i < pars.p.size(); ++i)
    {
      if (pars.saP[i] > 0)
        bwt_p[i] = pars.p[pars.saP[i] - 1];
      else
        bwt_p[i] = pars.p[pars.saP.size() - 1 - 1];
    }

    w_wt.construct(alphabet, bwt_p);
  }
};


class pfp_lce_support{
protected:
  pf_parsing& pfp;

  sdsl::bit_vector b_p; // Starting sdsl::bit_vector::rank_1_type rank_b_p;rank_b_p::bit_vector::select_1_type select_b_p;
  sdsl::bit_vector::rank_1_type rank_b_p;
  sdsl::bit_vector::select_1_type select_b_p;

public:


  // This has to be changed using pfp_dictionary and pfp_parse
  pfp_lce_support(pf_parsing& pfp_):
                    pfp(pfp_),
                    b_p(pfp.n,0)
  { //

    // Build the bitvector storing the position of the beginning of each phrase.
    b_p.resize(pfp.n); // all should be initialized at false by sdsl
    for(size_t i = 0; i < 3*64; ++i) b_p[i] = false; // bug in resize
    b_p[0] = true; // phrase_0 becomes phrase 1
    size_t i = 0;
    for(int j = 0; j < pfp.pars.p.size()-1; ++j){
      // p[i]: phrase_id
      assert(pfp.pars.p[j] != 0);
      // phrase_length: select_b_d(p[i]+1)-select_b_d(p[i]);
      i += pfp.dict.length_of_phrase(pfp.pars.p[j]) - pfp.w;
      b_p[i] = true;
    }
//    b_p[n] = true; // We set the last phrase beginning in the last w characters of the string


    // Build rank and select on Sp
    rank_b_p = sdsl::bit_vector::rank_1_type(&b_p);
    select_b_p = sdsl::bit_vector::select_1_type(&b_p);
  }

  // return the longest common prefix of suffix i and j of T
  size_t lce(size_t i, size_t j){
    if( i == j )
      return pfp.n-i;

    // find the phrases of which i and j belongs to
    size_t p_i = rank_b_p(i+1);// - 1; // rank of the phrase in T that i belongs to. Span: [1..|P|]
    size_t p_j = rank_b_p(j+1);// - 1; // rank of the phrase in T that j belongs to. Span: [1..|P|]
    // NOTE: i+1 since rank return the number of 1s in v[0..i-1]


    // NOTE: the -1 is because the bitvector p start with a 1 but the first phrase is in position 0 in p
    auto id_p_i = pfp.pars.p[p_i - 1]; // phrase_id of the phrase that i belongs to.
    auto id_p_j = pfp.pars.p[p_j - 1]; // phrase_id of the phrase that j belongs to.

    auto tmp_s_pi = select_b_p(p_i);
    auto tmp_s_pj = select_b_p(p_j);
    // find the length of the suffixes iof id_p_i and id_p_j starting at i and j
    // Length of the phrrase - length of the prefix.
    size_t len_suff_i_in_p_i = pfp.dict.length_of_phrase(id_p_i) - (i - select_b_p(p_i) + 1);
    size_t len_suff_j_in_p_j = pfp.dict.length_of_phrase(id_p_j) - (j - select_b_p(p_j) + 1);
    size_t k = std::min(len_suff_i_in_p_i,len_suff_j_in_p_j);

     // find the occurrence of i and j in the concatenation of the phrases of D
     // Starting position of the phrase in D + length of the prefix.
     size_t occ_in_p_i_in_D = pfp.dict.select_b_d(id_p_i) + (i - select_b_p(p_i) + 1);
     size_t occ_in_p_j_in_D = pfp.dict.select_b_d(id_p_j) + (j - select_b_p(p_j) + 1);

     // compute the LCP of the suffix up to the phrase boundaries.
     auto lcp_ppi_ppj = 0;
     if(occ_in_p_i_in_D == occ_in_p_j_in_D){
       // if two positions map to the each position in the same phrase,
       // we have to ignore the lcp in D because the lcp is the length of the suffix of the phrase.
       lcp_ppi_ppj = k;
     }else{
       auto lcp_left = std::min(pfp.dict.isaD[occ_in_p_i_in_D],pfp.dict.isaD[occ_in_p_j_in_D])+1;
       auto lcp_right = max(pfp.dict.isaD[occ_in_p_i_in_D],pfp.dict.isaD[occ_in_p_j_in_D]);
       size_t lcp_ppi_ppj_i = pfp.dict.rmq_lcp_D(lcp_left,lcp_right);
       lcp_ppi_ppj = pfp.dict.lcpD[lcp_ppi_ppj_i];
     }

     if(lcp_ppi_ppj < k){ // if the lcp of the suffixes up to the phrase boundaries is shorter than the shortest suffix length
       return lcp_ppi_ppj;
     }else{ // otherwise, find the first phrase where the suffixes differ

       k -= pfp.w; // this align the lcp values
       // Check if one of those is the last should be dobe by the $ in P
       auto p_ip1_in_sa = pfp.pars.isaP[p_i]; // isaP is 0-based
       auto p_jp1_in_sa = pfp.pars.isaP[p_j];

       auto lcp_left = std::min(p_ip1_in_sa,p_jp1_in_sa)+1;
       auto lcp_right = max(p_ip1_in_sa,p_jp1_in_sa);
       size_t lcp_pi_pj_i = pfp.pars.rmq_lcp_P(lcp_left,lcp_right);
       auto lcp_pi_pj = pfp.pars.lcpP[lcp_pi_pj_i];

       // size_t l_com_phrases = 0
       // for l in range(lcp_pi_pj):
       // #               print(p_i + 1 + l)
       // l_com_phrases += len(self.D[self.P[p_i + 1 + l]]) - w

       size_t l_com_phrases = select_b_p(p_i + 1 + lcp_pi_pj) - select_b_p(p_i + 1); // Check if there is some error it might be here


       auto a = pfp.pars.p[p_i + lcp_pi_pj]; // p is 0-based
       auto b = pfp.pars.p[p_j + lcp_pi_pj];

       assert (a != 0 and b != 0);
       // if a == "#" or b == "#": # This cannot be possible due to the last $ in T'
       // return k + l_com_phrases
       // else:
       // Compute the lcp between phrases a and b
       auto a_in_sa = pfp.dict.isaD[pfp.dict.select_b_d(a)]; // position of the phrase a in saD
       auto b_in_sa = pfp.dict.isaD[pfp.dict.select_b_d(b)]; // position of the phrase b in saD

       lcp_left = std::min(a_in_sa,b_in_sa) + 1;
       lcp_right = max(a_in_sa,b_in_sa);

       size_t lcp_a_b_i = pfp.dict.rmq_lcp_D(lcp_left, lcp_right);
       auto lcp_a_b = pfp.dict.lcpD[lcp_a_b_i];

       return k + l_com_phrases + lcp_a_b;
     }



  }

};



// NOTE: The dictionary must contain the last w characters of each phrase
template< class d_t, class saD_t, class isaD_t, class daD_t, class pfD_t, class lcpD_t>
sdsl::bit_vector compute_B_BWT(
                  const d_t* d, // the collection of the dictionary of the prefix-free parsing
                  const std::vector<saD_t>& saD, // the inverse suffix array of D
                  const std::vector<isaD_t>& isaD, // the inverse suffix array of D
                  const std::vector<daD_t>& daD, // the document array of D
                  const std::vector<pfD_t>& pfD, // the frequencies of phrases of D in T
                  const std::vector<lcpD_t>& lcpD,
                  const sdsl::bit_vector::rank_1_type& rank_b_d,
                  const sdsl::bit_vector::select_1_type& select_b_d,
                  const sdsl::bit_vector& b_d,
                  size_t w,  // size of the window of the parsing
                  size_t n // length of the text (including terminators)
                ){

      // Build the bitvector storing the position of the beginning of each phrase.
      sdsl::bit_vector b_bwt(n,false); // all should be initialized at false by sdsl

      assert(d[saD[0]] == EndOfDict);
      size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
      size_t j = 0;
      while (i < saD.size()) {
        saD_t sn = saD[i];
        // Check if the suffix has length at least w and is not the complete phrase.
        daD_t phrase = daD[i] + 1; assert(phrase > 0 && phrase < pfD.size()); // + 1 because daD is 0-based
        size_t suffix_length = select_b_d(rank_b_d(sn+1) +1) - sn - 1;
        if(b_d[sn] || suffix_length < w){
          ++i; // Skip
        }else{
          // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
          // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
          b_bwt[j++] = true; j += pfD[phrase]-1; // the next bits are 0s
          i++;
          if (i < saD.size()){
            saD_t new_sn = saD[i];
            daD_t new_phrase = daD[i] + 1; assert(new_phrase > 0 && new_phrase < pfD.size()); // + 1 because daD is 0-based
            size_t new_suffix_length = select_b_d(rank_b_d(new_sn+1) +1) - new_sn - 1;

            while (i < saD.size() && (lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length)){
                j += pfD[new_phrase];
                ++i;

                if (i < saD.size()){
                  new_sn = saD[i];
                  new_phrase = daD[i] + 1; assert(new_phrase > 0 && new_phrase < pfD.size()); // + 1 because daD is 0-based
                  new_suffix_length = select_b_d(rank_b_d(new_sn+1) +1) - new_sn - 1;
                }
            }
          }
        }
      }

      return b_bwt;
}

class pfp_sa {
protected:
  pf_parsing& pfp;
public:
  pfp_sa(pf_parsing & pfp_)
    : pfp(pfp_)
  {
    
  }
};

// TODO:
//  - Rozdeleni alphabet je potreba kontrolovat na pritomnost v te konkretni casti alphabety
void create_W_simple() {
  std::vector<std::string> dict {"##GATTAC", "ACAT#", "AGATA##", "T#GATAC", "T#GATTAG"};
  std::vector<uint32_t> parse {1, 2, 4, 2, 5, 3};
  std::vector<uint32_t> indices {1, 2, 3, 4, 5};

  std::vector<uint8_t> dict2 = {'#', '#', 'G', 'A', 'T', 'T', 'A', 'C', EndOfWord,
                                'A', 'C', 'A', 'T', '#', EndOfWord,
                                'A', 'G', 'A', 'T', 'A', '#', '#', EndOfWord,
                                'T', '#', 'G', 'A', 'T', 'A', 'C', EndOfWord,
                                'T', '#', 'G', 'A', 'T', 'T', 'A', 'G', EndOfWord
                                };
  sdsl::bit_vector b_d = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0, 0
                          };
  sdsl::bit_vector::rank_1_type rank_b_d(&b_d);
  sdsl::bit_vector::select_1_type select_b_d(&b_d);
  const size_t n_phrases = 5;
  const std::vector<size_t> phrase_length {8, 5, 7, 7, 8};

  // co-lexi compare function
  auto co_lexi_dict_cmp = [&](const uint32_t i, const uint32_t j) {
    auto i_start = select_b_d(i);
    auto i_end = i_start + phrase_length[i - 1] - 1;
    auto j_start = select_b_d(j);
    auto j_end = j_start + phrase_length[j - 1] - 1;

    auto i_r_begin = dict2.rend() - i_end - 1;
    auto i_r_end = dict2.rend() - i_start;
    auto j_r_begin = dict2.rend() - j_end - 1;
    auto j_r_end = dict2.rend() - j_start;

    return std::lexicographical_compare(i_r_begin, i_r_end, j_r_begin, j_r_end);
  };
  std::sort(indices.begin(), indices.end(), co_lexi_dict_cmp);

  std::cout << "0: " << select_b_d(1) << " | 1: " << select_b_d(2) << std::endl;
  std::cout << "last: " << select_b_d(5) << std::endl;

  for (const auto i : indices)
    std::cout << i << " | ";
  std::cout << std::endl;

  // create wavelet tree based on co-lexi sorted phrases
  pfp_wt wt(indices, parse);
  wt.print_leafs();
  std::cout << "wt[0]" << wt[0] << std::endl;
  std::cout << "wt[1]" << wt[1] << std::endl;
  std::cout << "wt[2]" << wt[2] << std::endl;
  std::cout << "wt[3]" << wt[3] << std::endl;
  std::cout << "wt[4]" << wt[4] << std::endl;
  std::cout << "wt[5]" << wt[5] << std::endl;
}

int main(int argc, char const *argv[]) {

  create_W_simple();
  return 0;

  if(argc < 2)
    error("input file required");

  std::string filename = argv[1];

  size_t w = 10;
  pf_parsing pf(filename,w);

  verbose("Dictionary size (bytes):       " , sizeof(pf.dict.d[0]) * pf.dict.d.size());
  verbose("Dictionary SA size (bytes):    " , sizeof(pf.dict.saD[0]) * pf.dict.saD.size());
  verbose("Dictionary ISA size (bytes):   " , sizeof(pf.dict.isaD[0]) * pf.dict.isaD.size());
  verbose("Dictionary LCP size (bytes):   " , sizeof(pf.dict.lcpD[0]) * pf.dict.lcpD.size());
  verbose("Dictionary DA size (bytes):    " , sizeof(pf.dict.daD[0]) * pf.dict.daD.size());
  verbose("Dictionary co-DA size (bytes): " , sizeof(pf.dict.colex_daD[0]) * pf.dict.colex_daD.size());
  verbose("Parsing size (bytes):          " , sizeof(pf.pars.p[0]) * pf.pars.p.size());
  verbose("Parsing SA size (bytes):       " , sizeof(pf.pars.saP[0]) * pf.pars.saP.size());
  verbose("Parsing ISA size (bytes):      " , sizeof(pf.pars.isaP[0]) * pf.pars.isaP.size());
  verbose("Parsing LCP size (bytes):      " , sizeof(pf.pars.lcpP[0]) * pf.pars.lcpP.size());


  size_t n = pf.n;

  verbose("Providing LCE support");
  elapsed_time(
    pfp_lce_support lce_ds(pf)
  );

  verbose("Computing W");
  elapsed_time(
    pfp_sa pfp_sa(pf)
  );

  // Building b_bwt

#if 0
// TEST lca_ds
std::vector<unsigned char> text;
read_fasta_file(filename.c_str(), text);
verbose("Text size: " , text.size());

// build lcp of the Text
text.push_back(0);
uint_t t_len = text.size();
std::vector<uint_t> saT(t_len); //uint_t *saT = new uint_t[t_len];
std::vector<int_t> lcpT(t_len); //int_t *lcpT = new int_t[t_len];

verbose("Computing SA, ISA, and LCP of the text");
// time_t  start = time(NULL);
sacak(&text[0],&saT[0],t_len);

// Computing isaT
std::vector<uint_t> isaT(t_len); //uint_t *isaT = new uint_t[t_len];
for(int i = 0; i < t_len; ++i){
  isaT[saT[i]] = i;
}

LCP_array(&text[0], isaT, saT, t_len, lcpT);
sdsl::rmq_succinct_sct<> rmq_lcp_T = sdsl::rmq_succinct_sct<>(&lcpT);

/*
// Compute the Dictionary
std::vector<std::string> di;

std::string s;
for(size_t i=0;i<saD.size();i++){
  if(d[i]==EndOfWord) {
    di.push_back(s);
    s.clear();
  }else{
    s.append(1,d[i]);
  }
}
*/

verbose("Testing LCE ds");
// for(int i = 1; i < text.size()-1 ; ++i){
//   for(int j = i+1; j < text.size(); ++j){
//     auto a = lce_ds.lce(saT[i], saT[j]);
//     auto b = lcpT[rmq_lcp_T(i+1,j)];//lcpT[i+1];
//     assert( a == b);
//   }
// }
for(int i = 1; i < text.size()-1 ; ++i){
    auto a = lce_ds.lce(saT[i], saT[i+1]);
    auto b = lcpT[i+1];//lcpT[rmq_lcp_T(i+1,j)];//lcpT[i+1];
    assert( a == b);
}

#endif



  return 0;

}
