/* pfp - prefix free parsing wt W support
    Copyright (C) 2020 Ondřej Cvacho

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
   \file wt.hpp
   \brief wt.hpp define and build the prefix-free parsing wavelet tree support data structures.
   \author Ondřej Cvacho
   \date 03/04/2020
*/


#ifndef _PFP_WT_HH
#define _PFP_WT_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

class pfp_wt {
public:
  using size_type = size_t;

  pfp_wt() {};
  pfp_wt(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {};

  virtual ~pfp_wt() {};

  virtual void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) = 0;

  virtual uint32_t operator[] (const size_type i) = 0;

  virtual size_type size() const = 0;

  virtual size_type rank(const size_type i, const uint32_t c) const = 0;
  virtual size_type select(const size_type i, const uint32_t c) const = 0;
  virtual size_type range_count (const uint32_t t, const uint32_t b, const size_type i) const = 0;
  virtual size_type range_count_2d (const uint32_t t, const size_type i) const = 0;
  size_type range_select (const uint32_t t, const uint32_t b, const size_type r) const {
    size_type lo = 1;
    size_type hi = size();
    assert (r > 0 && r <= range_count(t, b, hi));

    size_type ans = 0;

    // find r-th "filled column" in the interval (t, b)
    // binary search using rank_count
    while (lo <= hi) {
      const size_type i = (lo + hi) >> 1;
      const size_type m = range_count(t, b, i);

      if (m == r) {
        ans = i - 1;
        hi = i - 1;
      }

      if (m < r)
        lo = i + 1;
      else
        hi = i - 1;
    }

    return ans;
  }


  virtual size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const = 0;
  //! Load from a stream.
  virtual void load(std::istream &in) = 0;
private:
};

class pfp_wt_custom : public pfp_wt {
public:
  using wt_bv = sdsl::bit_vector;
  struct wt_node {
  
    wt_node (wt_node * parent = nullptr)
      : parent(parent)
    { }

    wt_bv bit_vector;
    wt_bv::rank_1_type bv_rank_1;
    wt_bv::select_1_type bv_select_1;
    wt_bv::select_0_type bv_select_0;

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

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
    {
      sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
      size_type written_bytes = 0;

      written_bytes += bit_vector.serialize(out, child, "bit_vector");
      written_bytes += bv_rank_1.serialize(out, child, "bv_rank_1");
      written_bytes += bv_select_1.serialize(out, child, "bv_select_1");
      written_bytes += bv_select_0.serialize(out, child, "bv_select_0");
      written_bytes += sdsl::write_member(phrase_id, out, child, "phrase_id");
      written_bytes += sdsl::write_member(is_leaf(), out, child, "is_leaf()");
      if(!is_leaf()){
        written_bytes += (*left).serialize(out, child, "left");
        written_bytes += (*right).serialize(out, child, "right");
      }
      
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
    }

    //! Load from a stream.
    void load(std::istream &in)
    {
      bit_vector.load(in);
      bv_rank_1.load(in, &bit_vector);
      bv_select_1.load(in, &bit_vector);
      bv_select_0.load(in, &bit_vector);
      sdsl::load(phrase_id, in);
      bool leaf;
      sdsl::read_member(leaf, in);
      if(!leaf){
        left = std::unique_ptr<wt_node>(new wt_node(std::addressof(*this)));
        right = std::unique_ptr<wt_node>(new wt_node(std::addressof(*this)));
        left->load(in);
        right->load(in);
      }
    }
  };

  struct leaf_info {
    size_type alphabet_index;
    wt_node * leaf_link = nullptr;
  };

  pfp_wt_custom()
    : root(new wt_node())
  { }

  pfp_wt_custom(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse)
    : root(new wt_node()), alphabet(sorted_alphabet) {
    for (size_type i = 0; i < alphabet.size(); ++i) {
      leafs[alphabet[i]].alphabet_index = i;
    }

    create_bwt_rec(*root, alphabet.size(), 0, parse);
  }

  void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {
    alphabet = sorted_alphabet;
    for (size_type i = 0; i < alphabet.size(); ++i) {
      leafs[alphabet[i]].alphabet_index = i;
    }

    create_bwt_rec(*root, alphabet.size(), 0, parse);
  }

  void print_leafs() const {
    print_leafs_rec(root);
    std::cout << std::endl;
  }

  // operator[] - get phrase ID on i-th position
  uint32_t operator[] (const size_type i) override {
    // assert(i < root->bit_vector.size());

    return get_phrase_id(*root, i); // WT[] is 0-based
  }

  size_type rank(const size_type i, const uint32_t c) const override {
    // assert(i > 0 && i <= root->bit_vector.size());

    if (leafs.count(c) <= 0) {
      return 0;
    }

    wt_node * t_node = std::addressof(*root);
    size_type j = i;
    uint32_t alphabet_size = alphabet.size();
    uint32_t alphabet_start = 0;
    while (!t_node->is_leaf()) {
      const auto idx = leafs.at(c).alphabet_index - alphabet_start;
      const auto tres = alphabet_size / 2;
      const auto rank_bit_1 = t_node->bv_rank_1(j);

      if (idx < tres) {
        // left child
        j -= rank_bit_1;
        alphabet_size = tres;
        t_node = std::addressof(*t_node->left);
      }
      else {
        // right child
        j = rank_bit_1;
        alphabet_size -= tres;
        alphabet_start += tres;
        t_node = std::addressof(*t_node->right);
      }
    }

    return j;
  }

  size_type select(const size_type i, const uint32_t c) const override {
    // assert (i > 0 && i <= rank(size(), c));

    wt_node * t_node = leafs.at(c).leaf_link;
    size_type j = i;
    while (!t_node->is_root()) {
      if (std::addressof(*(t_node->parent->left)) == t_node) {
        j = t_node->parent->bv_select_0(j) + 1;
      }
      else {
        j = t_node->parent->bv_select_1(j) + 1;
      }
      t_node = t_node->parent;
    }

    return j - 1; // return index 0-based
  }

  size_type size() const override {
    return root->bit_vector.size();
  }

  size_type range_count (const uint32_t t, const uint32_t b, const size_type i) const override {
    // assert(i > 0 && i <= root->bit_vector.size());

    // t and b intervals of leafs
    // for each leaf from (t, b) - rank(i, c)
    const auto count_b = range_count_2d(alphabet[b], i);
  
    if (t >= 1)
      return count_b - range_count_2d(alphabet[t - 1], i);
    
    return count_b;
  }

  size_type range_count_2d (const uint32_t t, const size_type i) const override {
    // assert(i > 0 && i <= root->bit_vector.size());

    wt_node * node = std::addressof(*root);
    const size_type alphabet_index = leafs.at(t).alphabet_index;
    size_type alphabet_start = 0;
    size_type alphabet_size = alphabet.size();

    size_type count = 0;
    size_type j = i;
    size_type rank_j = 0;

    while (!node->is_leaf()) {
      const auto idx = alphabet_index - alphabet_start;
      const auto tres = alphabet_size / 2;
      rank_j = node->bv_rank_1(j);

      if (idx < tres) {
        j -= rank_j;
        alphabet_size = tres;
        node = std::addressof(*node->left);
      }
      else {
        // all 0s add towards count
        count += j - rank_j;

        j = rank_j;
        alphabet_size -= tres;
        alphabet_start += tres;
        node = std::addressof(*node->right);
      }
    }
    // + case when we arrived to leaf
    count += j;

    return count;
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += sdsl::serialize(alphabet, out, child, "alphabet");
    written_bytes += root->serialize(out, child, "root");
    
    std::vector<uint32_t> tmp_phrase_ids;
    std::vector<size_type> tmp_alphbet_indices;


    for (const auto &kv : leafs){
      tmp_phrase_ids.push_back(kv.first);
      tmp_alphbet_indices.push_back(kv.second.alphabet_index);
    }

    written_bytes += sdsl::serialize(tmp_phrase_ids, out, child, "phrase_ids");
    written_bytes += sdsl::serialize(tmp_alphbet_indices, out, child, "alphbet_indices");


    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    sdsl::load(alphabet, in);
    root->load(in);

    std::vector<uint32_t> tmp_phrase_ids;
    std::vector<size_type> tmp_alphbet_indices;

    sdsl::load(tmp_phrase_ids, in);
    sdsl::load(tmp_alphbet_indices, in);

    for(size_type i = 0; i < tmp_phrase_ids.size(); ++i){
      leaf_info tmp_info;
      tmp_info.alphabet_index = tmp_alphbet_indices[i];
      leafs[tmp_phrase_ids[i]] = tmp_info;
    }

    tmp_phrase_ids.clear();
    tmp_alphbet_indices.clear();

    // Fill leaves pointers
    std::function<void(wt_node &, std::map<uint32_t, leaf_info> &)> vis;
    vis = [&vis](wt_node &node, std::map<uint32_t, leaf_info> &l) -> void { 
      if(node.is_leaf()){
        // leaf
        l[node.phrase_id].leaf_link = std::addressof(node);
        return;
      } else {
        vis(*(node.left),l);
        vis(*(node.right),l);
      }
    };

    vis(*root,leafs);
  }

private:
  std::vector<uint32_t> alphabet;
  // root node of wavelet tree
  std::unique_ptr<wt_node> root;
  std::map<uint32_t, leaf_info> leafs;

  uint32_t get_phrase_id(const wt_node & node, const size_type i) const {
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

  void create_bwt_rec(wt_node & w_structure, const uint32_t alpha_size, const uint32_t alpha_start, const std::vector<uint32_t> & parse) {
    // divide alphabet to two sets
    const auto tres = alpha_size / 2;
    if (tres == 0) {
      // leaf
      w_structure.phrase_id = alphabet[alpha_start];
      leafs[w_structure.phrase_id].leaf_link = std::addressof(w_structure);

      return;
    }

    // create bit_vector
    w_structure.bit_vector.resize(parse.size());
    std::vector<uint32_t> parse_left;
    std::vector<uint32_t> parse_right;

    for (size_type i = 0; i < parse.size(); i++) {
      // TODO: direct access data structure
      const auto idx = leafs.at(parse[i]).alphabet_index - alpha_start;

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
    w_structure.bv_select_0 = wt_bv::select_0_type(&w_structure.bit_vector);

    w_structure.left = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));
    w_structure.right = std::unique_ptr<wt_node>(new wt_node(std::addressof(w_structure)));

    create_bwt_rec(*(w_structure.left), tres, alpha_start, parse_left);
    create_bwt_rec(*(w_structure.right), alpha_size - tres, alpha_start + tres, parse_right);
  }

  void print_leafs_rec(const std::unique_ptr<wt_node> & node) const {
    if (node->is_leaf()) {
      std::cout << node->phrase_id << " | ";
    }
    else {
      print_leafs_rec(node->left);
      print_leafs_rec(node->right);
    }
  }
};

class pfp_wt_sdsl : public pfp_wt {
public:
  pfp_wt_sdsl()
  { }

  pfp_wt_sdsl(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) {
    construct(sorted_alphabet, parse);
  }

  void construct(const std::vector<uint32_t> & sorted_alphabet, const std::vector<uint32_t> & parse) override {
    sdsl::int_vector<> parse_translate(parse.size(), 0);
    translate.resize(sorted_alphabet.size(), 0);
    i_translate.resize(sorted_alphabet.size(), 0);

    for (size_type i = 0; i < sorted_alphabet.size(); ++i) {
      translate[sorted_alphabet[i] - 1] = i;
      i_translate[i] = sorted_alphabet[i];
    }

    for (size_type i = 0; i < parse.size(); ++i) {
      parse_translate[i] = translate[parse[i] - 1];
    }

    sdsl::construct_im(wt_i, parse_translate);
  }

  // operator[] - get phrase ID on i-th position
  uint32_t operator[] (const size_type i) override {
    return i_translate[wt_i[i]];
  }

  size_type size() const override {
    return wt_i.size();
  }

  size_type rank(const size_type i, const uint32_t c) const override {
    wt_i.rank(i, translate[c - 1]);
  }

  size_type select(const size_type i, const uint32_t c) const override {
    wt_i.select(i, translate[c - 1]);
  }

  size_type range_count (const uint32_t t, const uint32_t b, const size_type i) const override {
    const auto count_b = range_count_2d(b, i);
  
    if (t >= 1)
      return count_b - range_count_2d(t - 1, i);
    
    return count_b;
  }

  size_type range_count_2d (const uint32_t t, const size_type i) const override {
    return wt_i.range_search_2d(0, i - 1, 0, t, false).first;
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const override {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += wt_i.serialize(out, child, "wt_i");
    written_bytes += sdsl::serialize(i_translate, out, child, "i_translate");
    written_bytes += sdsl::serialize(translate, out, child, "translate");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in) override {
    wt_i.load(in);
    sdsl::load(i_translate, in);
    sdsl::load(translate, in);
  }

private:
  sdsl::wt_int<> wt_i;
  std::vector<uint32_t> i_translate;
  std::vector<uint32_t> translate;
};

#endif /* end of include guard: _PFP_WT_HH */
