/****************************************************************************
 * Copyright (C) from 2009 to Present EPAM Systems.
 *
 * This file is part of Indigo toolkit.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ***************************************************************************/

#include "molecule/molecule_aromatic_stereo.h"

#include "base_c/bitarray.h"
#include "molecule/molecule.h"

using namespace indigo;

AromaticStereoValidator::AromaticStereoValidator(Molecule* molecule) : _mol(molecule), _dearomatizations_ready(false)
{
}

void AromaticStereoValidator::_ensureDearomatizations()
{
    if (_dearomatizations_ready)
        return;

    AromaticityOptions options;
    Dearomatizer dearomatizer(*_mol, nullptr, options);
    dearomatizer.enumerateDearomatizations(_dearomatizations);
    _dearomatizations_ready = true;
}

bool AromaticStereoValidator::_collectCenterCandidates(int atom_idx, CenterCandidates& center)
{
    if (_mol == nullptr || atom_idx < 0 || atom_idx >= _mol->vertexEnd() || !_mol->hasVertex(atom_idx) ||
        _mol->getAtomAromaticity(atom_idx) != ATOM_AROMATIC)
        return false;

    const Vertex& vertex = _mol->getVertex(atom_idx);

    // Ordinary Indigo stereocenter validation remains the authority for
    // element, charge, degree, and bond-order configurations. The aromatic
    // validator only supplies concrete Kekule assignments and proves them
    // against the whole aromatic system.
    if (vertex.degree() <= 2 || vertex.degree() > 4)
        return false;

    Array<int> vertices;
    Array<int> aromatic_bonds;
    vertices.push(atom_idx);

    for (int i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i))
    {
        const int edge_idx = vertex.neiEdge(i);
        vertices.push(vertex.neiVertex(i));

        if (_mol->getBondOrder(edge_idx) == BOND_AROMATIC)
            aromatic_bonds.push(edge_idx);
    }

    // Actual aromatic bond participation, not lowercase SMILES spelling,
    // defines the exceptional path. This also covers atoms that a saver must
    // spell uppercase with explicit aromatic bonds.
    if (aromatic_bonds.size() == 0)
        return false;

    Molecule local;
    Array<int> mapping;
    local.makeSubmolecule(*_mol, vertices, &mapping, SKIP_ALL);

    const int local_atom_idx = mapping[atom_idx];
    Array<int> local_aromatic_bonds;
    for (int i = 0; i < aromatic_bonds.size(); i++)
    {
        const Edge& edge = _mol->getEdge(aromatic_bonds[i]);
        const int neighbor_idx = edge.beg == atom_idx ? edge.end : edge.beg;
        const int local_edge_idx = local.findEdgeIndex(local_atom_idx, mapping[neighbor_idx]);
        if (local_edge_idx < 0)
            return false;
        local_aromatic_bonds.push(local_edge_idx);
    }

    const int combinations = 1 << aromatic_bonds.size();
    for (int mask = 0; mask < combinations; mask++)
    {
        for (int i = 0; i < local_aromatic_bonds.size(); i++)
        {
            const int bond_order = (mask & (1 << i)) != 0 ? BOND_DOUBLE : BOND_SINGLE;
            local.setBondOrder_Silent(local_aromatic_bonds[i], bond_order);
        }

        if (!local.isPossibleStereocenter(local_atom_idx))
            continue;

        Assignment assignment;
        for (int i = 0; i < aromatic_bonds.size(); i++)
        {
            const int bond_order = (mask & (1 << i)) != 0 ? BOND_DOUBLE : BOND_SINGLE;
            assignment.bond_orders.push_back({aromatic_bonds[i], bond_order});
        }
        center.assignments.push_back(assignment);
    }

    return !center.assignments.empty();
}

void AromaticStereoValidator::_unfixCandidateBonds(DearomatizationMatcher& matcher, const std::vector<int>& newly_fixed_bonds,
                                                    std::vector<int>& fixed_bond_orders)
{
    for (auto i = newly_fixed_bonds.rbegin(); i != newly_fixed_bonds.rend(); ++i)
    {
        matcher.unfixBond(*i);
        fixed_bond_orders[*i] = 0;
    }
}

bool AromaticStereoValidator::_areCentersJointlyPossible(int center_idx, DearomatizationMatcher& matcher, std::vector<int>& fixed_bond_orders)
{
    if (center_idx == static_cast<int>(_centers.size()))
        return true;

    for (const auto& assignment : _centers[center_idx].assignments)
    {
        std::vector<int> newly_fixed_bonds;
        bool valid = true;

        try
        {
            for (const auto& bond : assignment.bond_orders)
            {
                const int fixed_order = fixed_bond_orders[bond.edge_idx];
                if (fixed_order != 0)
                {
                    if (fixed_order != bond.bond_order)
                    {
                        valid = false;
                        break;
                    }
                    continue;
                }

                if (!matcher.isAbleToFixBond(bond.edge_idx, bond.bond_order) || !matcher.fixBond(bond.edge_idx, bond.bond_order))
                {
                    valid = false;
                    break;
                }

                fixed_bond_orders[bond.edge_idx] = bond.bond_order;
                newly_fixed_bonds.push_back(bond.edge_idx);
            }

            if (valid && _areCentersJointlyPossible(center_idx + 1, matcher, fixed_bond_orders))
                return true;
        }
        catch (...)
        {
            _unfixCandidateBonds(matcher, newly_fixed_bonds, fixed_bond_orders);
            throw;
        }

        _unfixCandidateBonds(matcher, newly_fixed_bonds, fixed_bond_orders);
    }

    return false;
}

namespace
{
    bool hasSerializableTetrahedralStereo(Molecule& mol)
    {
        for (int i = mol.stereocenters.begin(); i != mol.stereocenters.end(); i = mol.stereocenters.next(i))
        {
            const int atom_idx = mol.stereocenters.getAtomIndex(i);
            if (mol.stereocenters.isTetrahydral(atom_idx) && mol.stereocenters.getType(atom_idx) >= MoleculeStereocenters::ATOM_AND)
                return true;
        }
        return false;
    }
} // namespace

void AromaticStereoValidator::suppressIncompatibleAromatization(Molecule& mol, const byte* proposed_aromatic_bonds,
                                                                std::vector<bool>& suppressed_bonds)
{
    suppressed_bonds.clear();

    if (!hasSerializableTetrahedralStereo(mol))
        return;

    suppressed_bonds.assign(mol.edgeEnd(), false);

    std::vector<bool> will_be_aromatic(mol.edgeEnd(), false);
    std::vector<bool> newly_proposed_aromatic(mol.edgeEnd(), false);
    bool has_newly_proposed_aromatic_bond = false;

    for (int e_idx = mol.edgeBegin(); e_idx < mol.edgeEnd(); e_idx = mol.edgeNext(e_idx))
    {
        const bool proposed = bitGetBit(proposed_aromatic_bonds, e_idx) != 0;
        const bool already_aromatic = mol.getBondOrder(e_idx) == BOND_AROMATIC;
        will_be_aromatic[e_idx] = proposed || already_aromatic;
        newly_proposed_aromatic[e_idx] = proposed && !already_aromatic;
        has_newly_proposed_aromatic_bond |= newly_proposed_aromatic[e_idx];
    }

    if (!has_newly_proposed_aromatic_bond)
        return;

    // Components are built from the aromatic graph that would exist after
    // applying the proposal. Existing aromatic bonds are included so extending
    // an already-aromatic system cannot bypass stereo compatibility checks.
    std::vector<int> vertex_component(mol.vertexEnd(), -1);
    std::vector<std::vector<int>> component_vertices;

    for (int v_idx = mol.vertexBegin(); v_idx < mol.vertexEnd(); v_idx = mol.vertexNext(v_idx))
    {
        if (vertex_component[v_idx] != -1)
            continue;

        bool has_aromatic_neighbor = false;
        const Vertex& vertex = mol.getVertex(v_idx);
        for (int i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i))
            if (will_be_aromatic[vertex.neiEdge(i)])
            {
                has_aromatic_neighbor = true;
                break;
            }

        if (!has_aromatic_neighbor)
            continue;

        const int component_idx = static_cast<int>(component_vertices.size());
        component_vertices.emplace_back();
        std::vector<int> pending;
        pending.push_back(v_idx);
        vertex_component[v_idx] = component_idx;

        while (!pending.empty())
        {
            const int current = pending.back();
            pending.pop_back();
            component_vertices[component_idx].push_back(current);

            const Vertex& current_vertex = mol.getVertex(current);
            for (int i = current_vertex.neiBegin(); i != current_vertex.neiEnd(); i = current_vertex.neiNext(i))
            {
                const int edge_idx = current_vertex.neiEdge(i);
                if (!will_be_aromatic[edge_idx])
                    continue;

                const int neighbor = current_vertex.neiVertex(i);
                if (vertex_component[neighbor] == -1)
                {
                    vertex_component[neighbor] = component_idx;
                    pending.push_back(neighbor);
                }
            }
        }
    }

    const int component_count = static_cast<int>(component_vertices.size());
    std::vector<bool> component_has_proposal(component_count, false);
    std::vector<std::vector<int>> component_stereocenters(component_count);

    for (int e_idx = mol.edgeBegin(); e_idx < mol.edgeEnd(); e_idx = mol.edgeNext(e_idx))
        if (newly_proposed_aromatic[e_idx])
        {
            const Edge& edge = mol.getEdge(e_idx);
            const int component_idx = vertex_component[edge.beg];
            if (component_idx >= 0)
                component_has_proposal[component_idx] = true;
        }

    for (int i = mol.stereocenters.begin(); i != mol.stereocenters.end(); i = mol.stereocenters.next(i))
    {
        const int atom_idx = mol.stereocenters.getAtomIndex(i);
        const int component_idx = atom_idx < static_cast<int>(vertex_component.size()) ? vertex_component[atom_idx] : -1;
        if (component_idx < 0 || !component_has_proposal[component_idx])
            continue;
        if (!mol.stereocenters.isTetrahydral(atom_idx) || mol.stereocenters.getType(atom_idx) < MoleculeStereocenters::ATOM_AND)
            continue;

        component_stereocenters[component_idx].push_back(atom_idx);
    }

    std::vector<bool> blocked_components(component_count, false);

    for (int component_idx = 0; component_idx < component_count; component_idx++)
    {
        if (component_stereocenters[component_idx].empty())
            continue;

        // Validate on the aromatic component plus every immediate neighbor of
        // its atoms. This preserves complete center degree/connectivity without
        // cloning the whole molecule.
        Array<int> local_vertices;
        std::vector<bool> included(mol.vertexEnd(), false);

        for (int atom_idx : component_vertices[component_idx])
        {
            if (!included[atom_idx])
            {
                included[atom_idx] = true;
                local_vertices.push(atom_idx);
            }

            const Vertex& vertex = mol.getVertex(atom_idx);
            for (int i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i))
            {
                const int neighbor = vertex.neiVertex(i);
                if (!included[neighbor])
                {
                    included[neighbor] = true;
                    local_vertices.push(neighbor);
                }
            }
        }

        Molecule local;
        Array<int> mapping;
        local.makeSubmolecule(mol, local_vertices, &mapping, SKIP_ALL);

        for (int e_idx = mol.edgeBegin(); e_idx < mol.edgeEnd(); e_idx = mol.edgeNext(e_idx))
        {
            if (!newly_proposed_aromatic[e_idx])
                continue;

            const Edge& edge = mol.getEdge(e_idx);
            if (vertex_component[edge.beg] != component_idx)
                continue;

            const int local_edge_idx = local.findEdgeIndex(mapping[edge.beg], mapping[edge.end]);
            if (local_edge_idx < 0)
                throw Exception("internal: proposed aromatic bond is absent from validation submolecule");
            local.setBondOrder(local_edge_idx, BOND_AROMATIC, true);
        }

        AromaticStereoValidator validator(&local);
        for (int atom_idx : component_stereocenters[component_idx])
            if (!validator.addStereocenter(mapping[atom_idx]))
            {
                blocked_components[component_idx] = true;
                break;
            }
    }

    for (int e_idx = mol.edgeBegin(); e_idx < mol.edgeEnd(); e_idx = mol.edgeNext(e_idx))
    {
        if (!newly_proposed_aromatic[e_idx])
            continue;

        const Edge& edge = mol.getEdge(e_idx);
        const int component_idx = vertex_component[edge.beg];
        if (component_idx >= 0 && blocked_components[component_idx])
            suppressed_bonds[e_idx] = true;
    }
}

bool AromaticStereoValidator::addStereocenter(int atom_idx)
{
    CenterCandidates center;
    if (!_collectCenterCandidates(atom_idx, center))
        return false;

    _ensureDearomatizations();
    _centers.push_back(center);

    DearomatizationMatcher matcher(_dearomatizations, *_mol, nullptr);
    std::vector<int> fixed_bond_orders(_mol->edgeEnd(), 0);
    const bool possible = _areCentersJointlyPossible(0, matcher, fixed_bond_orders);

    if (!possible)
        _centers.pop_back();

    return possible;
}
