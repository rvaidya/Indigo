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

#ifndef __molecule_aromatic_stereo_h__
#define __molecule_aromatic_stereo_h__

#include <vector>

#include "molecule/molecule_dearom.h"

namespace indigo
{
    class Molecule;

    // Validates tetrahedral stereocenters that participate in aromatic bond
    // systems by resolving their incident aromatic bonds to concrete Kekule
    // orders and proving that all accepted centers can coexist in one valid
    // whole-system dearomatization.
    //
    // The chemistry authority remains Molecule::isPossibleStereocenter();
    // this class adds no element/charge whitelist of its own.
    class AromaticStereoValidator
    {
    public:
        explicit AromaticStereoValidator(Molecule* molecule);

        bool addStereocenter(int atom_idx);

    private:
        struct BondOrder
        {
            int edge_idx;
            int bond_order;
        };

        struct Assignment
        {
            std::vector<BondOrder> bond_orders;
        };

        struct CenterCandidates
        {
            std::vector<Assignment> assignments;
        };

        void _ensureDearomatizations();
        bool _collectCenterCandidates(int atom_idx, CenterCandidates& center);
        void _unfixCandidateBonds(DearomatizationMatcher& matcher, const std::vector<int>& newly_fixed_bonds,
                                  std::vector<int>& fixed_bond_orders);
        bool _areCentersJointlyPossible(int center_idx, DearomatizationMatcher& matcher, std::vector<int>& fixed_bond_orders);

        Molecule* _mol;
        DearomatizationsStorage _dearomatizations;
        bool _dearomatizations_ready;
        std::vector<CenterCandidates> _centers;
    };
} // namespace indigo

#endif
