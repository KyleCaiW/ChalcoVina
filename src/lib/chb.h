/*

   Copyright (c) 2025-2026

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Wenhao Cai <klocmetreb@hotmail.com>, 
		   China Pharmaceutical University 
           
*/

#ifndef VINA_CHB_H
#define VINA_CHB_H

#include "common.h"
#include "atom_constants.h"
#include <string>
#include <vector>

struct model;
struct atom;

// ==================== Data Structures ====================

// CHBMAP structure
struct CHBTable {
    std::string name;
    unsigned int dims[3];
    fl origin[3];
    fl inv_step[3];
    size_t stride1;
    std::vector<fl> data;

    CHBTable() : dims{0, 0, 0}, origin{0, 0, 0}, inv_step{0, 0, 0}, stride1(0) {}
};

// Intermolecular Chalcogen bond donor structure
struct chalcogen_donor {
    sz atom_index;
    vec coords;
    sz neighbor1;
    sz neighbor2;
    vec neighbor1_coords;
    vec neighbor2_coords;
};

// Intermolecular Chalcogen bond donor structure (indices only, for caching)
struct chalcogen_donor_indices {
    sz atom_index;
    sz neighbor1;
    sz neighbor2;
};

// Intermolecular Chalcogen bond acceptor structure
struct chalcogen_acceptor {
    sz atom_index;
    vec coords;
    sz ad_type;
};

// Intramolecular Chalcogen bond acceptor structure
struct intra_chb_pair {
    sz donor_atom;
    sz acceptor_atom;
    sz donor_neighbor1;
    sz donor_neighbor2;
    sz acceptor_type;
};

struct chalcogen_bond_detail {
    sz acceptor_type;
    fl distance;
    fl polar_angle_deg;
    fl azimuthal_angle_deg;
    fl raw_energy;
    fl weighted_energy;
};

// ==================== ChB Interface ====================

// ChB map lookup and management interface
class CHBTableManager {
public:
    // Load all ChB maps from a .map file
    static bool load_maps_from_file(const std::string& path);

    // Get energy and spherical gradients from the loaded ChB maps
    static fl get_map_energy_gradient(fl distance, fl polar_angle, fl azimuthal_angle, sz acceptor_type,
                                     fl& grad_distance, fl& grad_polar_angle, fl& grad_azimuthal_angle);

private:
    // Helper to get grid position for interpolation
    static bool get_grid_position(fl coord, fl origin, fl inv_step, unsigned int dim,
                              int& idx0, int& idx1, fl& weight);

    // Trilinear interpolation for energy and gradient
    static fl trilinear_interpolation_with_gradient(const CHBTable& map,
                                              fl r, fl theta, fl phi,
                                              fl& grad_r, fl& grad_theta, fl& grad_phi);
    
    static CHBTable* table_cache[AD_TYPE_SIZE];
};

// ==================== ChB Utility Functions ====================

// Atom type identification
bool is_divalent_sulfur(const atom& a);
bool is_chalcogen_acceptor(const atom& a);

// Functions to find chalcogen bond donors and acceptors
std::vector<chalcogen_donor> find_chalcogen_donors(const model& m);
std::vector<chalcogen_donor_indices> find_chalcogen_donor_indices(const model& m);
std::vector<chalcogen_donor> update_chalcogen_donors_from_indices(const model& m,
                                   const std::vector<chalcogen_donor_indices>& cached_indices);
std::vector<chalcogen_acceptor> find_chalcogen_acceptors(const model& m);

std::vector<chalcogen_bond_detail> analyze_chalcogen_bonds(const model& m,
    fl weight_O, fl weight_S, fl weight_N);

// Functions to calculate chalcogen bond energy and gradients
fl calculate_chalcogen_bond_energy(const model& m, fl weight_O, fl weight_S, fl weight_N, std::vector<vec>& gradients,
                                   const std::vector<chalcogen_acceptor>& cached_acceptors,
                                   const std::vector<chalcogen_donor_indices>& cached_donor_indices);

#endif
