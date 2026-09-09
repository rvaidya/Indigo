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
