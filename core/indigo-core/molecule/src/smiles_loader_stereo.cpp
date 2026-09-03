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

#include <cctype>
#include <memory>
#include <regex>
#include <unordered_set>

#include "base_cpp/scanner.h"
#include "graph/cycle_basis.h"
#include "molecule/elements.h"
#include "molecule/molecule.h"
#include "molecule/molecule_arom.h"
#include "molecule/molecule_dearom.h"
#include "molecule/molecule_stereocenters.h"
#include "molecule/query_molecule.h"
#include "molecule/smiles_loader.h"

using namespace indigo;

namespace
{
    struct AromaticStereoValidationContext
    {
        DearomatizationsStorage dearomatizations;
        std::unique_ptr<DearomatizationMatcher> matcher;
    };

    DearomatizationMatcher& getDearomatizationMatcher(Molecule& mol, AromaticStereoValidationContext& context)
    {
        if (context.matcher == nullptr)
        {
            AromaticityOptions options;
            Dearomatizer dearomatizer(mol, nullptr, options);
            dearomatizer.enumerateDearomatizations(context.dearomatizations);
            context.matcher = std::make_unique<DearomatizationMatcher>(context.dearomatizations, mol, nullptr);
        }
        return *context.matcher;
    }

    bool isPossibleAromaticStereocenter(Molecule* mol, int atom_idx, AromaticStereoValidationContext& context)
    {
        if (mol == nullptr || atom_idx < 0 || atom_idx >= mol->vertexEnd() || !mol->hasVertex(atom_idx) ||
            mol->getAtomAromaticity(atom_idx) != ATOM_AROMATIC)
            return false;

        const Vertex& vertex = mol->getVertex(atom_idx);
        if (vertex.degree() <= 2 || vertex.degree() > 4)
            return false;

        Array<int> vertices;
        Array<int> aromatic_bonds;
        vertices.push(atom_idx);

        for (int i = vertex.neiBegin(); i != vertex.neiEnd(); i = vertex.neiNext(i))
        {
            const int edge_idx = vertex.neiEdge(i);
            const int bond_order = mol->getBondOrder(edge_idx);

            if (bond_order == BOND_TRIPLE)
                return false;

            vertices.push(vertex.neiVertex(i));
            if (bond_order == BOND_AROMATIC)
                aromatic_bonds.push(edge_idx);
        }

        if (aromatic_bonds.size() == 0)
            return false;

        Molecule local;
        Array<int> mapping;
        local.makeSubmolecule(*mol, vertices, &mapping, SKIP_ALL);

        const int local_atom_idx = mapping[atom_idx];
        Array<int> local_aromatic_bonds;
        for (int i = 0; i < aromatic_bonds.size(); i++)
        {
            const Edge& edge = mol->getEdge(aromatic_bonds[i]);
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

            DearomatizationMatcher& matcher = getDearomatizationMatcher(*mol, context);
            Array<int> fixed_bonds;
            bool valid = true;

            try
            {
                for (int i = 0; i < aromatic_bonds.size(); i++)
                {
                    const int edge_idx = aromatic_bonds[i];
                    const int bond_order = (mask & (1 << i)) != 0 ? BOND_DOUBLE : BOND_SINGLE;

                    if (!matcher.isAbleToFixBond(edge_idx, bond_order))
                    {
                        valid = false;
                        break;
                    }

                    if (!matcher.fixBond(edge_idx, bond_order))
                    {
                        valid = false;
                        break;
                    }
                    fixed_bonds.push(edge_idx);
                }
            }
            catch (...)
            {
                for (int i = 0; i < fixed_bonds.size(); i++)
                    matcher.unfixBond(fixed_bonds[i]);
                throw;
            }

            for (int i = 0; i < fixed_bonds.size(); i++)
                matcher.unfixBond(fixed_bonds[i]);

            if (valid)
                return true;
        }

        return false;
    }
} // namespace

void SmilesLoader::_calcStereocenters()
{
    int i, j;
    AromaticStereoValidationContext aromatic_stereo_context;

    for (i = 0; i < _atoms.size(); i++)
    {
        if (_atoms[i].chirality == 0)
            continue;

        if (_bmol->getVertex(i).degree() == 2) // allene stereo center
        {
            int subst[4];
            int subst2[4];
            int left, right;
            bool pure_h[4];

            if (!MoleculeAlleneStereo::possibleCenter(*_bmol, i, left, right, subst, pure_h))
            {
                if (!stereochemistry_options.ignore_errors)
                    throw Error("chirality on atom %d makes no sense", i);
                continue;
            }

            int parity = 3 - _atoms[i].chirality;

            for (j = 0; j < 4; j++)
                if (subst[j] == -1)
                    subst2[j] = -1;
                else
                    subst2[j] = subst[j];

            // Daylight doc says: Hydrogens attached to substituted allene-like atoms
            // are taken to be immediately following that atom
            if (subst2[1] == -1)
                subst2[1] = left;
            if (subst2[3] == -1)
                subst2[3] = right;

            if (subst2[1] < subst2[0])
            {
                std::swap(subst2[1], subst2[0]);
                parity = 3 - parity;
            }

            if (subst2[3] < subst2[2])
            {
                std::swap(subst2[3], subst2[2]);
                parity = 3 - parity;
            }

            // move hydrogens from [0] and [2] to [1] and [3] respectively
            if (pure_h[0])
            {
                if (subst[1] == -1)
                    throw Error("unexpected: subst[1] = -1");
                std::swap(subst[0], subst[1]);
                parity = 3 - parity;
            }
            if (pure_h[2])
            {
                if (subst[3] == -1)
                    throw Error("unexpected: subst[3] = -1");
                std::swap(subst[2], subst[3]);
                parity = 3 - parity;
            }

            _bmol->allene_stereo.add(i, left, right, subst, parity);
        }
        else // ordinary tetrahedral stereo center
        {
            int pyramid[4] = {-1, -1, -1, -1};
            int counter = 0;
            int h_index = -1;

            if (_atoms[i].parent != -1)
                pyramid[counter++] = _atoms[i].parent;

            if (_atoms[i].neighbors.size() == 3)
            {
                h_index = counter;
                pyramid[counter++] = -1;
            }

            for (j = _atoms[i].neighbors.begin(); j != _atoms[i].neighbors.end(); j = _atoms[i].neighbors.next(j))
            {
                int nei = _atoms[i].neighbors.at(j);

                if (counter >= 4)
                {
                    if (!stereochemistry_options.ignore_errors)
                        throw Error("too many bonds for chiral atom %d", i);
                    break;
                }

                if (nei != _atoms[i].parent)
                    pyramid[counter++] = nei;
            }

            if (j != _atoms[i].neighbors.end())
                continue;

            if (counter < 3)
            {
                if (!stereochemistry_options.ignore_errors)
                    throw Error("only %d bonds for chiral atom %d", counter, i);
                continue;
            }

            if (counter == 4)
            {
                j = pyramid[0];
                pyramid[0] = pyramid[1];
                pyramid[1] = pyramid[2];
                pyramid[2] = pyramid[3];
                pyramid[3] = j;

                if (h_index == 0)
                    h_index = 3;
                else if (h_index > 0)
                    h_index--;
            }

            if (h_index >= 0)
            {
                if (counter != 4)
                {
                    if (!stereochemistry_options.ignore_errors)
                        throw Error("implicit hydrogen not allowed with %d neighbor atoms", counter - 1);
                    continue;
                }

                bool parity = true;

                for (j = h_index; j < 3; j++)
                {
                    std::swap(pyramid[j], pyramid[j + 1]);
                    parity = !parity;
                }

                if (!parity)
                    std::swap(pyramid[0], pyramid[1]);
            }

            if (_atoms[i].chirality == 2)
                std::swap(pyramid[0], pyramid[1]);

            bool possible_stereocenter = _bmol->isPossibleStereocenter(i);
            if (!possible_stereocenter && _atoms[i].aromatic)
                possible_stereocenter = isPossibleAromaticStereocenter(_mol, i, aromatic_stereo_context);

            if (!possible_stereocenter)
            {
                if (!stereochemistry_options.ignore_errors)
                    throw Error("chirality not possible on atom #%d", i);
                continue;
            }

            _bmol->addStereocenters(i, MoleculeStereocenters::ATOM_ABS, 0, pyramid);
        }
    }
}

void SmilesLoader::_calcCisTrans()
{
    QS_DEF(Array<int>, dirs);
    int i;

    dirs.clear();

    for (i = 0; i < _bonds.size(); i++)
        dirs.push(_bonds[i].dir);

    // there could be bonds added to stereocenters
    for (; i < _bmol->edgeEnd(); i++)
        dirs.push(0);

    _bmol->buildFromSmilesCisTrans(dirs.ptr());
    if (_qmol != 0)
    {
        for (i = 0; i < _bonds.size(); i++)
            if (_bmol->cis_trans.getParity(i) != 0)
                _qmol->setBondStereoCare(i, true);
    }
}

void SmilesLoader::_validateStereoCenters()
{
    AromaticStereoValidationContext aromatic_stereo_context;
    for (int i = _bmol->stereocenters.begin(); i < _bmol->stereocenters.end(); i = _bmol->stereocenters.next(i))
    {
        auto atom_idx = _bmol->stereocenters.getAtomIndex(i);
        bool possible_stereocenter = _bmol->isPossibleStereocenter(atom_idx);
        if (!possible_stereocenter && _bmol->stereocenters.isTetrahydral(atom_idx))
            possible_stereocenter = isPossibleAromaticStereocenter(_mol, atom_idx, aromatic_stereo_context);

        if (possible_stereocenter || _bmol->stereocenters.isAtropisomeric(atom_idx))
            continue;
        if (!stereochemistry_options.ignore_errors)
            throw Error("atom %d is not a stereocenter", atom_idx);
    }
}

void SmilesLoader::_addExplicitHForStereo()
{
    for (int i = 0; i < _atoms.size(); i++)
    {
        if ((_atoms[i].chirality > 0) && (_bmol->getVertex(i).degree() == 2) && (_atoms[i].hydrogens == 1))
        {
            _AtomDesc& atom = _atoms.push(_neipool);
            _BondDesc* bond = &_bonds.push();

            atom.label = ELEM_H;
            int exp_h_idx = _mol->addAtom(atom.label);

            bond->beg = i;
            bond->end = _atoms.size() - 1;
            bond->type = BOND_SINGLE;
            bond->index = _mol->addBond_Silent(bond->beg, bond->end, bond->type);

            _atoms[i].neighbors.add(exp_h_idx);
            _atoms[exp_h_idx].neighbors.add(i);
            _atoms[exp_h_idx].parent = i;

            _atoms[i].hydrogens = 0;
        }
    }
}

void SmilesLoader::_addLigandsForStereo()
{
    bool add_explicit_h = false;
    int num_ligands = 0;

    for (int i = 0; i < _atoms.size(); i++)
    {
        if ((_atoms[i].chirality > 0) && (_bmol->getVertex(i).degree() < 3) && !_isAlleneLike(i))
        {
            if (_atoms[i].hydrogens == 1)
            {
                add_explicit_h = true;
                num_ligands = 3 - _bmol->getVertex(i).degree() - _atoms[i].hydrogens;
            }
            else
                num_ligands = 3 - _bmol->getVertex(i).degree();

            for (int j = 0; j < num_ligands; j++)
            {
                _AtomDesc& atom = _atoms.push(_neipool);
                _BondDesc* bond = &_bonds.push();
                std::unique_ptr<QueryMolecule::Atom> qatom;

                if (add_explicit_h)
                    qatom.reset(QueryMolecule::Atom::nicht(new QueryMolecule::Atom(QueryMolecule::ATOM_NUMBER, ELEM_H)));
                else
                    qatom.reset(QueryMolecule::Atom::oder(QueryMolecule::Atom::nicht(new QueryMolecule::Atom(QueryMolecule::ATOM_NUMBER, ELEM_H)),
                                                          new QueryMolecule::Atom(QueryMolecule::ATOM_NUMBER, ELEM_H)));

                std::unique_ptr<QueryMolecule::Bond> qbond = std::make_unique<QueryMolecule::Bond>(QueryMolecule::BOND_ORDER, BOND_SINGLE);

                atom.star_atom = true;
                int any_atom_idx = _qmol->addAtom(qatom.release());

                bond->beg = i;
                bond->end = _atoms.size() - 1;
                bond->type = BOND_SINGLE;
                bond->dir = 0;
                bond->topology = 0;
                bond->index = _qmol->addBond(i, any_atom_idx, qbond.release());

                _atoms[i].neighbors.add(any_atom_idx);
                _atoms[any_atom_idx].neighbors.add(i);
                _atoms[any_atom_idx].parent = i;
            }

            if (_atoms[i].hydrogens == 1)
            {
                _AtomDesc& atom = _atoms.push(_neipool);
                _BondDesc* bond = &_bonds.push();

                std::unique_ptr<QueryMolecule::Atom> qatom = std::make_unique<QueryMolecule::Atom>(QueryMolecule::ATOM_NUMBER, ELEM_H);
                std::unique_ptr<QueryMolecule::Bond> qbond = std::make_unique<QueryMolecule::Bond>(QueryMolecule::BOND_ORDER, BOND_SINGLE);

                atom.label = ELEM_H;
                int exp_h_idx = _qmol->addAtom(qatom.release());

                bond->beg = i;
                bond->end = _atoms.size() - 1;
                bond->type = BOND_SINGLE;
                bond->dir = 0;
                bond->topology = 0;
                bond->index = _qmol->addBond(i, exp_h_idx, qbond.release());

                _atoms[i].neighbors.add(exp_h_idx);
                _atoms[exp_h_idx].neighbors.add(i);
                _atoms[exp_h_idx].parent = i;

                _atoms[i].hydrogens = 0;
                _qmol->getAtom(i).removeConstraints(QueryMolecule::ATOM_TOTAL_H);
            }
        }
    }
}

bool SmilesLoader::_isAlleneLike(int i)
{
    if (_bmol->getVertex(i).degree() == 2)
    {
        int subst[4];
        int left, right;
        bool pure_h[4];

        if (MoleculeAlleneStereo::possibleCenter(*_bmol, i, left, right, subst, pure_h))
            return true;
    }
    return false;
}
