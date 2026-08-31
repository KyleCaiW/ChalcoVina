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

#include "chb.h"
#include "model.h"
#include "atom_constants.h"
#include "common.h"
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <cstring>

// Constraints
const fl MIN_DISTANCE_SQ = 6.25; // 2.5*2.5
const fl MAX_DISTANCE_SQ = 20.25; // 4.5*4.5
const fl MAX_POLAR_ANGLE_DEG = 30.0;

// Static table storage
static std::unordered_map<std::string, CHBTable> chb_tables;

// Static array cache for O(1) lookup - indexed by acceptor AD type
CHBTable* CHBTableManager::table_cache[AD_TYPE_SIZE] = {nullptr};

bool is_divalent_sulfur(const atom& a) {
    return a.ad == AD_TYPE_S && a.bonds.size() == 2;
}

bool is_chalcogen_acceptor(const atom& a) {
    return (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA || a.ad == AD_TYPE_SA);
    //return (a.ad == AD_TYPE_OA);
}

// Helper function to map AD type to table cache index
inline int get_table_cache_index(sz ad_type) {
    switch(ad_type) {
        case AD_TYPE_OA: return 0;
        case AD_TYPE_NA: return 1;
        case AD_TYPE_SA: return 2;
        default: return -1;  // Not a chalcogen acceptor
    }
}

bool CHBTableManager::get_grid_position(fl coord, fl origin, fl inv_step, unsigned int dim,
                              int& idx0, int& idx1, fl& weight) {
    if (inv_step <= 0) {
        return false;
    }

    const fl normalized = (coord - origin) * inv_step;

    if (normalized < 0.0) {
        return false;
    }

    const fl max_norm = static_cast<fl>(dim - 1.0);
    if (normalized > max_norm + epsilon_fl) {
        return false;
    }

    idx0 = static_cast<int>(normalized);
    idx1 = idx0 + 1;
    weight = normalized - idx0;

    return true;
}

fl CHBTableManager::trilinear_interpolation_with_gradient(const CHBTable& table,
                                                      fl r, fl theta, fl phi,
                                                      fl& grad_r, fl& grad_theta, fl& grad_phi) {
    int r0, r1, t0, t1, p0, p1;
    fl wr, wt, wp;

    if (!CHBTableManager::get_grid_position(r,     table.origin[0], table.inv_step[0], table.dims[0], r0, r1, wr) ||
        !CHBTableManager::get_grid_position(theta, table.origin[1], table.inv_step[1], table.dims[1], t0, t1, wt) ||
        !CHBTableManager::get_grid_position(phi,   table.origin[2], table.inv_step[2], table.dims[2], p0, p1, wp))
    {
        grad_r = grad_theta = grad_phi = 0.0;
        return 0.0;
    }

    const size_t base_r0 = static_cast<size_t>(r0) * table.stride1;
    const size_t base_r1 = static_cast<size_t>(r1) * table.stride1;
    const size_t base_t0 = static_cast<size_t>(t0) * table.dims[2];
    const size_t base_t1 = static_cast<size_t>(t1) * table.dims[2];

    const fl* data_ptr = table.data.data();

    const fl f000 = data_ptr[base_r0 + base_t0 + p0];
    const fl f100 = data_ptr[base_r1 + base_t0 + p0];
    const fl f010 = data_ptr[base_r0 + base_t1 + p0];
    const fl f110 = data_ptr[base_r1 + base_t1 + p0];
    const fl f001 = data_ptr[base_r0 + base_t0 + p1];
    const fl f101 = data_ptr[base_r1 + base_t0 + p1];
    const fl f011 = data_ptr[base_r0 + base_t1 + p1];
    const fl f111 = data_ptr[base_r1 + base_t1 + p1];

    const fl mr = 1.0 - wr;
    const fl mt = 1.0 - wt;
    const fl mp = 1.0 - wp;

    const fl energy =
        f000 * mr * mt * mp +
        f100 * wr * mt * mp +
        f010 * mr * wt * mp +
        f110 * wr * wt * mp +
        f001 * mr * mt * wp +
        f101 * wr * mt * wp +
        f011 * mr * wt * wp +
        f111 * wr * wt * wp;

    const fl grad_r_local =
        f000 * (-1) * mt * mp +
        f100 *   1 * mt * mp +
        f010 * (-1) * wt * mp +
        f110 *   1 * wt * mp +
        f001 * (-1) * mt * wp +
        f101 *   1 * mt * wp +
        f011 * (-1) * wt * wp +
        f111 *   1 * wt * wp;

    const fl grad_theta_local =
        f000 * mr * (-1) * mp +
        f100 * wr * (-1) * mp +
        f010 * mr *   1 * mp +
        f110 * wr *   1 * mp +
        f001 * mr * (-1) * wp +
        f101 * wr * (-1) * wp +
        f011 * mr *   1 * wp +
        f111 * wr *   1 * wp;

    const fl grad_phi_local =
        f000 * mr * mt * (-1) +
        f100 * wr * mt * (-1) +
        f010 * mr * wt * (-1) +
        f110 * wr * wt * (-1) +
        f001 * mr * mt *   1 +
        f101 * wr * mt *   1 +
        f011 * mr * wt *   1 +
        f111 * wr * wt *   1;

    grad_r = grad_r_local * table.inv_step[0];
    grad_theta = grad_theta_local * table.inv_step[1];
    grad_phi = grad_phi_local * table.inv_step[2];

    return energy;
}

fl CHBTableManager::get_map_energy_gradient(fl distance, fl polar_angle_deg, fl azimuthal_angle_deg, sz acceptor_type,
                                       fl& grad_distance, fl& grad_polar_angle, fl& grad_azimuthal_angle) {
    int cache_index = get_table_cache_index(acceptor_type);
    if (cache_index < 0 || cache_index >= 3) {
        grad_distance = grad_polar_angle = grad_azimuthal_angle = 0.0;
        return 0.0;
    }

    CHBTable* table = table_cache[cache_index];
    if (!table) {
        grad_distance = grad_polar_angle = grad_azimuthal_angle = 0.0;
        return 0.0;
    }

    fl energy = trilinear_interpolation_with_gradient(*table, distance, polar_angle_deg, azimuthal_angle_deg,
                                         grad_distance, grad_polar_angle, grad_azimuthal_angle);
    
    // Convert angle gradients: dE/d(rad) = dE/d(deg) * 180/pi
    grad_polar_angle      = grad_polar_angle * 180.0 / M_PI;
    grad_azimuthal_angle  = grad_azimuthal_angle * 180.0 / M_PI;
    
    return energy;
}

// Read the ChbGrid.map database file
bool CHBTableManager::load_maps_from_file(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return false;
    }

    char magic[6];
    unsigned int num_tables;

    file.read(magic, 6);
    if (!file || std::string(magic, 6) != "CHBMAP") {
        return false;
    }

    file.read(reinterpret_cast<char*>(&num_tables), sizeof(unsigned int));
    if (!file) {
        return false;
    }
    
    struct TableIndexEntry {
        char name[32];
        uint64_t offset;
        uint64_t size;
    };

    if (num_tables == 0) {
        chb_tables.clear();
        return true;
    }

    std::vector<TableIndexEntry> index_entries(num_tables);
    file.read(reinterpret_cast<char*>(index_entries.data()), num_tables * sizeof(TableIndexEntry));

    if (!file) {
        return false;
    }

    chb_tables.clear();
    for (const auto& entry : index_entries) {
        CHBTable table;
        table.name = std::string(entry.name, strnlen(entry.name, 32));

        file.seekg(entry.offset);
        if (!file) {
            continue;
        }
        
        char header_buffer[60];
        file.read(header_buffer, 60);

        if (!file) {
            continue;
        }

        memcpy(table.dims,   header_buffer,      3 * sizeof(unsigned int));
        memcpy(table.origin, header_buffer + 12, 3 * sizeof(fl));
        
        fl step[3];
        memcpy(step,         header_buffer + 36, 3 * sizeof(fl));

        for (int i = 0; i < 3; ++i) {
            table.inv_step[i] = (step[i] > 1e-9) ? 1.0 / step[i] : 0.0;
        }

        table.stride1 = static_cast<size_t>(table.dims[1]) * table.dims[2];
        size_t num_points = static_cast<size_t>(table.dims[0]) * table.dims[1] * table.dims[2];

        if (num_points > 0) {
            table.data.resize(num_points);
            file.read(reinterpret_cast<char*>(table.data.data()), num_points * sizeof(fl));

            if (!file) {
                continue;
            }
        }

        chb_tables[table.name] = std::move(table);
    }

    // Map: S-O -> index 0, S-N -> index 1, S-S -> index 2
    auto it0 = chb_tables.find("S-O");
    table_cache[0] = (it0 != chb_tables.end()) ? &it0->second : nullptr;

    auto it1 = chb_tables.find("S-N");
    table_cache[1] = (it1 != chb_tables.end()) ? &it1->second : nullptr;

    auto it2 = chb_tables.find("S-S");
    table_cache[2] = (it2 != chb_tables.end()) ? &it2->second : nullptr;

    if (!table_cache[0]) std::cerr << "Warning: missing ChB map S-O\n";
    if (!table_cache[1]) std::cerr << "Warning: missing ChB map S-N\n";
    if (!table_cache[2]) std::cerr << "Warning: missing ChB map S-S\n";

    return true;
}

std::vector<chalcogen_donor> find_chalcogen_donors(const model& m) {
    std::vector<chalcogen_donor> donors;
    const atomv& atoms = m.get_atoms();

    for (sz i = 0; i < atoms.size(); ++i) {
        const atom& a = atoms[i];

        if (!m.is_atom_in_ligand(i)) {
            continue;
        }

        if (is_divalent_sulfur(a)) {
            chalcogen_donor donor;
            donor.atom_index = i;
            donor.coords = m.get_coords(i);
            donor.neighbor1 = a.bonds[0].connected_atom_index.i;
            donor.neighbor2 = a.bonds[1].connected_atom_index.i;
            donor.neighbor1_coords = m.get_coords(donor.neighbor1);
            donor.neighbor2_coords = m.get_coords(donor.neighbor2);
            donors.push_back(donor);
        }
    }

    return donors;
}

std::vector<chalcogen_donor_indices> find_chalcogen_donor_indices(const model& m) {
    std::vector<chalcogen_donor_indices> donor_indices;
    const atomv& atoms = m.get_atoms();

    for (sz i = 0; i < atoms.size(); ++i) {
        const atom& a = atoms[i];

        if (!m.is_atom_in_ligand(i)) {
            continue;
        }

        if (is_divalent_sulfur(a)) {
            chalcogen_donor_indices indices;
            indices.atom_index = i;
            indices.neighbor1 = a.bonds[0].connected_atom_index.i;
            indices.neighbor2 = a.bonds[1].connected_atom_index.i;
            donor_indices.push_back(indices);
        }
    }

    return donor_indices;
}

// update donors from indices
std::vector<chalcogen_donor> update_chalcogen_donors_from_indices(const model& m,
                                   const std::vector<chalcogen_donor_indices>& cached_indices) {
    std::vector<chalcogen_donor> donors;
    donors.reserve(cached_indices.size());

    for (const auto& indices : cached_indices) {
        chalcogen_donor donor;
        donor.atom_index = indices.atom_index;
        donor.neighbor1 = indices.neighbor1;
        donor.neighbor2 = indices.neighbor2;

        donor.coords = m.get_coords(indices.atom_index);
        donor.neighbor1_coords = m.get_coords(indices.neighbor1);
        donor.neighbor2_coords = m.get_coords(indices.neighbor2);

        donors.push_back(donor);
    }

    return donors;
}

std::vector<chalcogen_acceptor> find_chalcogen_acceptors(const model& m) {
    std::vector<chalcogen_acceptor> acceptors;
    const atomv& grid_atoms = m.get_grid_atoms();

    for (sz j = 0; j < grid_atoms.size(); ++j) {
        const atom& acceptor_atom = grid_atoms[j];

        if (is_chalcogen_acceptor(acceptor_atom)) {
            chalcogen_acceptor acceptor;
            acceptor.atom_index = j;
            acceptor.coords = acceptor_atom.coords;
            acceptor.ad_type = acceptor_atom.ad;
            acceptors.push_back(acceptor);
        }
    }

    return acceptors;
}

fl dihedral(const vec& p1, const vec& p2, const vec& p3, const vec& p4)
{
    vec b1 = p2 - p1;
    vec b2 = p3 - p2;
    vec b3 = p4 - p3;

    vec n1 = cross_product(b1, b2);
    vec n2 = cross_product(b2, b3);
    vec b2_norm = (1.0 / b2.norm()) * b2;
    vec m1 = cross_product(n1, b2_norm);

    fl x = n1 * n2;
    fl y = m1 * n2;

    return std::atan2(y, x);
}

// Main energy calculation functions
fl calculate_chalcogen_bond_energy(const model& m, fl weight_O, fl weight_S, fl weight_N, std::vector<vec>& gradients,
                                   const std::vector<chalcogen_acceptor>& cached_acceptors,
                                   const std::vector<chalcogen_donor_indices>& cached_donor_indices) {
    fl total_energy = 0.0;

    std::vector<chalcogen_donor> donors = update_chalcogen_donors_from_indices(m, cached_donor_indices);

    if (donors.empty() || cached_acceptors.empty()) {
        return 0.0;
    }

    for (const auto& donor : donors) {
        const vec& B = donor.coords;
        const vec& A = donor.neighbor1_coords;
        const vec& C = donor.neighbor2_coords;

        vec AB = B - A;
        vec CB = B - C;

        fl ab_norm = AB.norm();
        fl cb_norm = CB.norm();

        if (ab_norm < epsilon_fl || cb_norm < epsilon_fl) continue;

        vec AB_norm = (1.0 / ab_norm) * AB;
        vec CB_norm = (1.0 / cb_norm) * CB;

        for (const auto& acceptor : cached_acceptors) {
            // Select weight based on acceptor type
            fl weight;
            switch (acceptor.ad_type) {
                case AD_TYPE_OA: weight = weight_O; break;
                case AD_TYPE_SA: weight = weight_S; break;
                case AD_TYPE_NA: weight = weight_N; break;
                default: continue; // Skip unknown acceptor types
            }

            if (weight == 0.0) continue;

            const vec& D = acceptor.coords;
            vec BD = D - B;

            fl dist_sq = BD.norm_sqr();
            if (dist_sq > MAX_DISTANCE_SQ || dist_sq < MIN_DISTANCE_SQ) {
                continue;
            }

            fl distance = std::sqrt(dist_sq);
            vec BD_norm = (1.0 / distance) * BD;

            fl cos_theta_AB = AB_norm * BD_norm;
            fl cos_theta_CB = CB_norm * BD_norm;

            cos_theta_AB = std::max(static_cast<fl>(-1.0), std::min(static_cast<fl>(1.0), cos_theta_AB));
            cos_theta_CB = std::max(static_cast<fl>(-1.0), std::min(static_cast<fl>(1.0), cos_theta_CB));

            bool use_CB_side    = (cos_theta_CB >= cos_theta_AB);
            fl cos_polar = use_CB_side ? cos_theta_CB : cos_theta_AB;
            fl polar_angle_rad = std::acos(cos_polar);
            fl polar_angle_deg = polar_angle_rad * (180.0 / M_PI);
            
            if (polar_angle_deg > MAX_POLAR_ANGLE_DEG) {
                continue;
            }
            
            fl azimuthal_angle_rad;
            if (use_CB_side) {
                // D-B-C-A
                azimuthal_angle_rad = dihedral(D, B, C, A);
            } else {
                // D-B-A-C
                azimuthal_angle_rad = dihedral(D, B, A, C);
            }

            fl azimuthal_angle_deg = std::abs(azimuthal_angle_rad * (180.0 / M_PI)); 

            // Get energy and gradients from chalcogen bond map using geometric parameters
            fl dE_dr, dE_dtheta, dE_dphi;
            fl energy = CHBTableManager::get_map_energy_gradient(distance, polar_angle_deg, azimuthal_angle_deg, acceptor.ad_type,
                                                               dE_dr, dE_dtheta, dE_dphi);
            total_energy += energy * weight;

            // ====================================================================
            // Cartesian Gradient Calculation for atoms A, B, and C
            // ∇(E) = (∂E/∂r)∇(r) + (∂E/∂θ)∇(θ) + (∂E/∂φ)∇(φ)
            // ====================================================================
            vec grad_A(0, 0, 0), grad_B(0, 0, 0), grad_C(0, 0, 0);

            // --- Term 1: Radial Gradient (∂E/∂r * ∇(r)) ---
            // r = |D - B|
            vec grad_B_r = -1.0 * BD_norm;
            grad_B += dE_dr * grad_B_r;

            // --- Term 2: Polar Angle Gradient (∂E/∂θ * ∇(θ)) ---
            // ∂θ/∂cosθ = -1 / sinθ
            // ∇_x θ = (∂θ/∂cosθ) * ∇_x(cosθ)
            fl sin_polar_angle = std::sin(polar_angle_rad);
            if (sin_polar_angle > epsilon_fl) {
                fl inv_sin = 1.0 / sin_polar_angle;
                // coeff combines chain rule terms: dE/dTheta * dTheta/dCosTheta
                // dTheta/dCosTheta = -1 / sin(theta)
                fl common_coef = -dE_dtheta * inv_sin;

                // Identify the two vectors defining the angle:
                // u: vector from Neighbor (A or C) to Donor (B).
                // v: vector from Donor (B) to Acceptor (D).

                vec u, v;
                fl r_u, r_v;
                
                v = BD_norm;
                r_v = distance;

                // Determine which side to use based on theta_AB and theta_CB
                if (use_CB_side) {
                    u = CB_norm;
                    r_u = cb_norm;
                } else {
                    u = AB_norm;
                    r_u = ab_norm;
                }

                fl cos_theta = cos_polar;

                // Calculate vectors perpendicular to u and v in the plane of the angle
                vec perp_v_on_u = v - cos_theta * u;
                vec perp_u_on_v = u - cos_theta * v;
                
                vec term_neighbor = (1.0 / r_u) * perp_v_on_u;
                vec term_acceptor = (1.0 / r_v) * perp_u_on_v;

                // Gradient on Neighbor (A or C)
                // ∇_Neighbor E = - (dE/dTheta / sin) * (perp_v_on_u / r_u)
                // = common_coef * term_neighbor
                vec grad_term_neighbor = -1.0 * common_coef * term_neighbor;
                
                // Gradient on Acceptor (D) -- Implicitly needed for B's balance
                // ∇_Acceptor E = - (dE/dTheta / sin) * (perp_u_on_v / r_v)
                // = common_coef * term_acceptor
                vec grad_term_acceptor = common_coef * term_acceptor;

                // Gradient on Donor (B)
                // By translational invariance: ∇_B = -(∇_Neighbor + ∇_Acceptor)
                vec grad_term_B = -1.0 * (grad_term_neighbor + grad_term_acceptor);

                // Apply to gradients
                if (use_CB_side) {
                    grad_C += grad_term_neighbor;
                } else {
                    grad_A += grad_term_neighbor;
                }
                
                grad_B += grad_term_B;
            }

            // --- Term 3: Azimuthal Angle Gradient (dE/dphi * grad_phi) ---
            if (std::abs(azimuthal_angle_rad) > epsilon_fl && std::abs(dE_dphi) > epsilon_fl) {
                const fl sign_phi = (azimuthal_angle_rad >= 0) ? 1.0 : -1.0; // derivative of |phi|

                auto accumulate_dihedral_grad = [&](const vec& p1, const vec& p2, const vec& p3, const vec& p4,
                                                    vec& grad_p2, vec& grad_p3, vec& grad_p4) {
                    const vec b1 = p2 - p1;
                    const vec b2 = p3 - p2;
                    const vec b3 = p4 - p3;

                    const fl b2_len = b2.norm();
                    if (b2_len < epsilon_fl) return;

                    const vec n1 = cross_product(b1, b2);
                    const vec n2 = cross_product(b2, b3);
                    
                    const fl n1_norm_sq = n1.norm_sqr();
                    const fl n2_norm_sq = n2.norm_sqr();

                    if (n1_norm_sq < epsilon_fl || n2_norm_sq < epsilon_fl) return;

                    const vec b2n = (1.0 / b2_len) * b2;
                    const vec m1  = cross_product(n1, b2n);

                    const fl x = n1 * n2;
                    const fl y = m1 * n2;
                    const fl denom = x*x + y*y;
                    if (denom < epsilon_fl) return;

                    const fl dphi_dx = -y / denom;
                    const fl dphi_dy =  x / denom;

                    const vec gx_n1 = dphi_dx * n2;
                    const vec gx_n2 = dphi_dx * n1;

                    const vec gy_n1   = dphi_dy * cross_product(b2n, n2);
                    const vec gy_n2   = dphi_dy * m1;
                    const vec gy_b2n  = dphi_dy * cross_product(n2, n1);

                    const vec g_n1 = gx_n1 + gy_n1;
                    const vec g_n2 = gx_n2 + gy_n2;

                    const vec grad_b1 = cross_product(b2, g_n1);

                    vec grad_b2 = cross_product(g_n1, b1);
                    grad_b2 += cross_product(b3, g_n2);

                    const vec grad_b3 = cross_product(g_n2, b2);

                    const fl proj = b2n * gy_b2n;
                    grad_b2 += (1.0 / b2_len) * (gy_b2n - proj * b2n);

                    const vec grad_p1 = -1.0 * grad_b1;
                    const vec grad_p2_local = grad_b1 - grad_b2;
                    const vec grad_p3_local = grad_b2 - grad_b3;
                    const vec grad_p4_local = grad_b3;

                    const fl factor = dE_dphi * sign_phi;

                    grad_p2 += factor * grad_p2_local;
                    grad_p3 += factor * grad_p3_local;
                    grad_p4 += factor * grad_p4_local;
                };

                if (use_CB_side) {
                    // Dihedral defined as D-B-C-A
                    accumulate_dihedral_grad(D, B, C, A, grad_B, grad_C, grad_A);
                } else {
                    // Dihedral defined as D-B-A-C
                    accumulate_dihedral_grad(D, B, A, C, grad_B, grad_A, grad_C);
                }
            }

            // Add weighted gradients to global gradients vector
            gradients[donor.atom_index] += weight * grad_B;
            gradients[donor.neighbor1]  += weight * grad_A;
            gradients[donor.neighbor2]  += weight * grad_C;
        }
    }

    return total_energy;
}

std::vector<chalcogen_bond_detail> analyze_chalcogen_bonds(const model& m,
    fl weight_O, fl weight_S, fl weight_N) {
    std::vector<chalcogen_bond_detail> details;

    std::vector<chalcogen_donor_indices> cached_donor_indices = find_chalcogen_donor_indices(m);
    std::vector<chalcogen_donor> donors = update_chalcogen_donors_from_indices(m, cached_donor_indices);
    std::vector<chalcogen_acceptor> acceptors = find_chalcogen_acceptors(m);

    if (donors.empty() || acceptors.empty()) {
        return details;
    }

    for (const auto& donor : donors) {
        const vec& B = donor.coords;
        const vec& A = donor.neighbor1_coords;
        const vec& C = donor.neighbor2_coords;

        vec AB = B - A;
        vec CB = B - C;

        fl ab_norm = AB.norm();
        fl cb_norm = CB.norm();

        if (ab_norm < epsilon_fl || cb_norm < epsilon_fl) continue;

        vec AB_norm = (1.0 / ab_norm) * AB;
        vec CB_norm = (1.0 / cb_norm) * CB;

        for (const auto& acceptor : acceptors) {
            fl weight;
            switch (acceptor.ad_type) {
                case AD_TYPE_OA: weight = weight_O; break;
                case AD_TYPE_SA: weight = weight_S; break;
                case AD_TYPE_NA: weight = weight_N; break;
                default: continue;
            }

            if (weight == 0.0) continue;

            const vec& D = acceptor.coords;
            vec BD = D - B;

            fl dist_sq = BD.norm_sqr();
            if (dist_sq > MAX_DISTANCE_SQ || dist_sq < MIN_DISTANCE_SQ) {
                continue;
            }

            fl distance = std::sqrt(dist_sq);
            vec BD_norm = (1.0 / distance) * BD;

            fl cos_theta_AB = AB_norm * BD_norm;
            fl cos_theta_CB = CB_norm * BD_norm;

            cos_theta_AB = std::max(static_cast<fl>(-1.0), std::min(static_cast<fl>(1.0), cos_theta_AB));
            cos_theta_CB = std::max(static_cast<fl>(-1.0), std::min(static_cast<fl>(1.0), cos_theta_CB));

            bool use_CB_side = (cos_theta_CB >= cos_theta_AB);
            fl cos_polar = use_CB_side ? cos_theta_CB : cos_theta_AB;
            fl polar_angle_rad = std::acos(cos_polar);
            fl polar_angle_deg = polar_angle_rad * (180.0 / M_PI);

            if (polar_angle_deg > MAX_POLAR_ANGLE_DEG) {
                continue;
            }

            fl azimuthal_angle_rad;
            if (use_CB_side) {
                azimuthal_angle_rad = dihedral(D, B, C, A);
            } else {
                azimuthal_angle_rad = dihedral(D, B, A, C);
            }

            fl azimuthal_angle_deg = std::abs(azimuthal_angle_rad * (180.0 / M_PI));

            fl dE_dr, dE_dtheta, dE_dphi;
            fl raw_energy = CHBTableManager::get_map_energy_gradient(distance, polar_angle_deg,
                azimuthal_angle_deg, acceptor.ad_type, dE_dr, dE_dtheta, dE_dphi);

            fl weighted_energy = raw_energy * weight;

            chalcogen_bond_detail detail;
            detail.acceptor_type = acceptor.ad_type;
            detail.distance = distance;
            detail.polar_angle_deg = polar_angle_deg;
            detail.azimuthal_angle_deg = azimuthal_angle_deg;
            detail.raw_energy = raw_energy;
            detail.weighted_energy = weighted_energy;

            details.push_back(detail);
        }
    }

    return details;
}
