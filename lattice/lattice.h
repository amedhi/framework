/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2016-01-11 16:33:51
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-02-23 00:23:06
*----------------------------------------------------------------------------*/
#include <iostream>
#include <string>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../constants.h"
#include <Eigen/Dense>

namespace lattice {

/*---------------lattice types-----------------*/
enum class lattice_type {
  UNDEFINED, SQUARE, CHAIN
};

/*---------------Lattice site class-----------------*/
using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
//using BrvaisIdx = Eigen::Matrix<unsigned, 3, 1>;

class Site 
{
public:
  // ctors
  Site() {}
  Site(const int& type, const int& atomid, const Vector3i& bravindex, const Vector3d& coord, 
    const Vector3d& cell_coord);
  ~Site() {}
  // setter functions
  void static reset_count(void) { num_site=0; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void reset_coord(const Vector3d& v) { coord_=v; }
  void reset_cell_coord(const Vector3d& v) { cell_coord_=v; }
  void translate_by(const int& id_offset, const Vector3i& bravindex_offset, const Vector3d& coord_offset); 

  // getter functions
  unsigned id(void) const {return id_;}
  unsigned type(void) const {return type_;}
  unsigned atomid(void) const {return atomid_;}
  Vector3i bravindex(void) const { return bravindex_; }
  Vector3d coord(void) const {return coord_;}
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Site& site);
private:
  static unsigned num_site;
  unsigned id_ {0};
  unsigned type_ {0};
  unsigned atomid_ {0};
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3d coord_ {Vector3d(0.0, 0.0, 0.0)};
  Vector3d cell_coord_ {Vector3d(0.0, 0.0, 0.0)};
};

/*---------------Lattice bond class-----------------*/
class Bond 
{
public:
  // ctors
  Bond() {}
  Bond(const unsigned& type, const unsigned& ngb, const Vector3i& bravindex, const unsigned& src_id, 
    const Vector3i& src_offset, const unsigned& tgt_id, const Vector3i& tgt_offset);
  ~Bond() {}
  // setter functions
  void static reset_count(void) { num_bond=0; }
  void reset_id(const unsigned& id) { id_=id; }
  void reset_src_offset(const Vector3i& idx) { src_offset_=idx; }
  void reset_tgt_offset(const Vector3i& idx) { tgt_offset_=idx; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void shift_target_ids(const int& id_offset) { src_ += id_offset; tgt_ += id_offset; }
  void translate_by(const Vector3i& bravindex_offset) { bravindex_ += bravindex_offset; } 
  void connect(const int& src_id, const Vector3i& src_offset, const int& tgt_id, const Vector3i& tgt_offset);
  // getter functions
  unsigned id(void) const { return id_; }
  unsigned src_id(void) const { return src_; }
  unsigned tgt_id(void) const { return tgt_; }
  Vector3i bravindex(void) const { return bravindex_; }
  Vector3i src_offset(void) const { return src_offset_; }
  Vector3i tgt_offset(void) const { return tgt_offset_; }
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Bond& bond);
private:
  static unsigned num_bond;
  unsigned id_ {0};
  unsigned type_ {0};
  unsigned ngb_ {0};
  unsigned src_ {0}; 
  unsigned tgt_ {0}; 
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3i src_offset_ {Vector3i(0, 0, 0)};
  Vector3i tgt_offset_ {Vector3i(0, 0, 0)};
};

/*---------------Unitcell class-----------------*/
class Unitcell 
{
public:
  // ctors
  Unitcell() {}
  ~Unitcell() {}
  // setter functions
  void clear(void); 
  void set_basis(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3);
  int add_site(const Site& s) { sites.push_back(s); return sites.back().id(); }
  int add_site(const unsigned& type, const unsigned& atomid, const Vector3d& site_coord); 
  int add_site(const Site& s, const Vector3i& bravindex, const Vector3d& cell_coord);
  int add_bond(const Bond& b) { bonds.push_back(b); return bonds.back().id(); }
  int add_bond(const unsigned& type, const unsigned& ngb, const unsigned& src_id, const Vector3i& src_offset,
    const unsigned& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const Bond& bond, const Vector3i& src_offset, const Vector3i& tgt_offset);
  void reset_a1(const Vector3d& av1) { a1=av1; }
  void reset_a2(const Vector3d& av2) { a2=av2; }
  void reset_a3(const Vector3d& av3) { a2=av3; }
  void reset(const std::vector<Site>& new_sites, const std::vector<Bond>& new_bonds); 
  // getter functions
  Vector3d vector_a1(void) const { return a1; }
  Vector3d vector_a2(void) const { return a2; }
  Vector3d vector_a3(void) const { return a3; }
  Vector3i bravindex(void) const { return bravindex_; }
  Vector3d coord(void) const {return coord_;}
  unsigned num_site(void) const { return sites.size(); }
  unsigned num_bond(void) const { return bonds.size(); }
  void translate_by(const Vector3i& bravindex_offset, const int& cell_id_offset);
  void rotate_by(const Eigen::Matrix3d& matrix);
  Site site(const unsigned& i) const { return sites[i]; }
  Bond bond(const unsigned& i) const { return bonds[i]; }
private:
  int id {0};
  unsigned max_site_type_val {0};
  unsigned max_bond_type_val {0};
  unsigned max_neighb_val {0};
  Vector3d a1 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a2 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a3 {Vector3d(0.0, 0.0, 0.0)};
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3d coord_ {Vector3d(0.0, 0.0, 0.0)};
};

/*---------------spatial dimension type-----------------*/
enum class boundary_type {open, periodic, antiperiodic};

/*---------------Lattcie class defintion-----------------*/
class Lattice 
{
public:
  // ctors
  Lattice() {}
  Lattice(const input::Parameters& parms) { construct(parms); }
  ~Lattice() {}
  // setter functions
  int construct(const input::Parameters& parms);
  // getter functions

  // other methods 
  Vector3i boundary_wrap(const Vector3i& cell_idx) const;
  Vector3i get_next_bravindex(const Vector3i& current_index) const;
  Unitcell get_translated_cell(const Vector3i& bravindex_offset) const;
  int mapped_site_id(const unsigned& local_id, const Vector3i& bravindex) const;
  bool connect_bond(Bond& bond) const;
private:
  struct Extent {unsigned size; boundary_type bc; boundary_type periodicity;};
  enum Dimension {dim1, dim2, dim3};
  lattice_type id {lattice_type::UNDEFINED};
  std::string name {""};

  // unitcell & lattice dimensions
  Unitcell unitcell;
  Extent extent[3] = {Extent{1, boundary_type::open, boundary_type::open}, 
                      Extent{1, boundary_type::open, boundary_type::open},
                      Extent{1, boundary_type::open, boundary_type::open}
                     };

  // copy of initial lattice dimensions
  Extent copy_extent[3] {Extent{1, boundary_type::open, boundary_type::open}, 
                         Extent{1, boundary_type::open, boundary_type::open},
                         Extent{1, boundary_type::open, boundary_type::open}
                        };
  
  // number of unit cells in total and in one layer (for symmetrized lattice)
  unsigned num_total_cells {1};
  unsigned num_layer_cells {1};

  // helper functions
  int define_lattice(void); 
  int finalize_lattice(void); 
  boundary_type boundary_condition(std::string& bc) const;
  Eigen::Matrix3d rotation_matrix(const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime);
};



} // end namespace lattice









