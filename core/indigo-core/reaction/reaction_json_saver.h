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

#ifndef __reaction_json_saver__
#define __reaction_json_saver__

#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "base_cpp/exception.h"

namespace indigo
{

    class Output;
    class BaseReaction;
    class BaseMolecule;
    class Vec2f;
    class ReactionJsonSaver
    {
    public:
        explicit ReactionJsonSaver(Output& output);
        ~ReactionJsonSaver();

        void saveReaction(BaseReaction& rxn);

        DECL_ERROR;

    protected:
        Output& _output;

    private:
        ReactionJsonSaver(const ReactionJsonSaver&); // no implicit copy
        static void _getBounds(BaseMolecule& mol, Vec2f& min_vec, Vec2f& max_vec, double scale);
    };

} // namespace indigo

#endif
