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

#include "graph/automorphism_search.h"

#include "graph/embedding_enumerator.h"
#include "graph/graph_decomposer.h"

#include <limits>
#include <chrono>
#include <cstdint>
#include <cstdio>
using namespace indigo;

IMPL_ERROR(AutomorphismSearch, "automorphism search");
IMPL_TIMEOUT_EXCEPTION(AutomorphismSearch, "automorphism search");

CP_DEF(AutomorphismSearch);

namespace
{
    struct ComponentMatchContext
    {
        AutomorphismSearch* search;
        const Array<int>* sub_to_graph;
        const Array<int>* super_to_graph;
    };

    struct SearchProfile
    {
        uint64_t dispatches = 0;
        uint64_t refine_calls = 0;
        uint64_t refine_ns = 0;
        uint64_t target_calls = 0;
        uint64_t target_ns = 0;
        uint64_t prune_calls = 0;
        uint64_t prune_ns = 0;
        uint64_t compare_calls = 0;
        uint64_t compare_ns = 0;
        uint64_t automorphism_checks = 0;
        int components = 0;
        int classes = 0;
        int seeded = 0;
        int generators = 0;
    };

    SearchProfile search_profile;
    using ProfileClock = std::chrono::steady_clock;

}


AutomorphismSearch::AutomorphismSearch()
    : CP_INIT, TL_CP_GET(_call_stack), TL_CP_GET(_lab), TL_CP_GET(_ptn), TL_CP_GET(_graph), TL_CP_GET(_mapping), TL_CP_GET(_inv_mapping), TL_CP_GET(_degree),
      TL_CP_GET(_tcells), TL_CP_GET(_fix), TL_CP_GET(_mcr), TL_CP_GET(_moved_vertices), TL_CP_GET(_seeded_component_automorphisms),
      TL_CP_GET(_generators_by_vertex), TL_CP_GET(_fixed_generator_count), TL_CP_GET(_active), TL_CP_GET(_workperm), TL_CP_GET(_workperm2), TL_CP_GET(_bucket),
      TL_CP_GET(_count), TL_CP_GET(_firstlab), TL_CP_GET(_canonlab), TL_CP_GET(_orbits), TL_CP_GET(_fixedpts), TL_CP_GET(_work_active_cells),
      TL_CP_GET(_edge_ranks_in_refine), TL_CP_GET(_edge_rank_cache), TL_CP_GET(_edge_counts), TL_CP_GET(_edge_count_touched),
      TL_CP_GET(_generator_seen_epoch), TL_CP_GET(_long_prune_candidates), TL_CP_GET(_lab_position), TL_CP_GET(_cell_starts), TL_CP_GET(_ptn_change_stack)
{
    getcanon = true;
    compare_vertex_degree_first = true;
    refine_reverse_degree = false;
    refine_by_sorted_neighbourhood = false;
    worksize = 100;

    context = 0;
    cb_vertex_rank = 0;
    cb_vertex_cmp = 0;
    cb_check_automorphism = 0;
    cb_compare_mapped = 0;
    cb_automorphism = 0;
    cb_edge_rank = 0;
    context_automorphism = 0;
    _given_graph = 0;
    _long_prune_epoch = 0;
    ignored_vertices = 0;

    _cancellation_handler = getCancellationHandler();

    _call_stack.clear();
}

AutomorphismSearch::~AutomorphismSearch()
{
}

void AutomorphismSearch::_prepareGraph(Graph& graph)
{
    QS_DEF(Array<int>, buckets);
    QS_DEF(Array<int>, ranks);
    int i;

    ranks.clear();
    buckets.clear();

    _graph.clear();
    _mapping.clear();
    _degree.clear();

    _ptn.clear();

    _inv_mapping.clear_resize(graph.vertexEnd());
    _degree.clear_resize(graph.vertexEnd());
    _degree.zerofill();

    if (cb_vertex_rank != 0 && cb_vertex_cmp != 0)
        throw Error("both vertex rank and vertex compare callbacks specified");

    _given_graph = &graph;

    if (cb_vertex_cmp != 0)
    {
        for (i = graph.vertexBegin(); i != graph.vertexEnd(); i = graph.vertexNext(i))
            if (ignored_vertices == 0 || !ignored_vertices[i])
                _mapping.push(i);

        for (i = graph.edgeBegin(); i != graph.edgeEnd(); i = graph.edgeNext(i))
        {
            const Edge& edge = graph.getEdge(i);

            if (ignored_vertices == 0 || (!ignored_vertices[edge.beg] && !ignored_vertices[edge.end]))
            {
                _degree[edge.beg]++;
                _degree[edge.end]++;
            }
        }

        _mapping.qsort(_cmp_vertices, this);

        int rank = 0;

        for (i = 0; i < _mapping.size(); i++)
        {
            if (i > 0 && _cmp_vertices(_mapping[i], _mapping[i - 1], this) != 0)
                rank++;

            ranks.push(rank);
        }
    }
    else
    {
        for (i = graph.vertexBegin(); i != graph.vertexEnd(); i = graph.vertexNext(i))
        {
            if (ignored_vertices != 0 && ignored_vertices[i])
                continue;

            int rank = 0;

            if (cb_vertex_rank != 0)
                rank = cb_vertex_rank(graph, i, context);

            ranks.push(rank);
            _mapping.push(i);
        }
    }

    for (i = 0; i < _mapping.size(); i++)
    {
        _graph.addVertex();
        _inv_mapping[_mapping[i]] = i;
        _ptn.push(AUTOMORPHISM_INFINITY);

        while (buckets.size() <= ranks[i])
            buckets.push(0);

        buckets[ranks[i]]++;
    }

    for (i = graph.edgeBegin(); i != graph.edgeEnd(); i = graph.edgeNext(i))
    {
        const Edge& edge = graph.getEdge(i);

        if (ignored_vertices != 0 && (ignored_vertices[edge.beg] || ignored_vertices[edge.end]))
            continue;

        int beg = _inv_mapping[edge.beg];
        int end = _inv_mapping[edge.end];

        _graph.addEdge(beg, end);
    }

    int start = 0;

    for (i = 0; i < buckets.size(); i++)
    {
        if (buckets[i] == 0)
            continue;

        int end = start + buckets[i];

        buckets[i] = start;

        _ptn[end - 1] = 0;

        start = end;
    }

    _n = _graph.vertexCount();

    _lab.clear_resize(_n);

    for (i = 0; i < _n; i++)
        _lab[buckets[ranks[i]]++] = i;
}

bool AutomorphismSearch::_componentVertexMatch(Graph& subgraph, Graph& supergraph, const int* /*core_sub*/, int sub_idx, int super_idx, void* userdata)
{
    ComponentMatchContext& match = *(ComponentMatchContext*)userdata;
    AutomorphismSearch& self = *match.search;

    if (subgraph.getVertex(sub_idx).degree() != supergraph.getVertex(super_idx).degree())
        return false;

    int graph_sub = match.sub_to_graph->at(sub_idx);
    int graph_super = match.super_to_graph->at(super_idx);
    if (graph_sub < 0 || graph_super < 0 || !self._graph.hasVertex(graph_sub) || !self._graph.hasVertex(graph_super))
        throw Error("internal: incomplete disconnected component vertex mapping");

    int original_sub = self._mapping[graph_sub];
    int original_super = self._mapping[graph_super];

    if (self.cb_vertex_cmp != 0 && self.cb_vertex_cmp(*self._given_graph, original_sub, original_super, self.context) != 0)
        return false;

    if (self.cb_vertex_rank != 0 &&
        self.cb_vertex_rank(*self._given_graph, original_sub, self.context) != self.cb_vertex_rank(*self._given_graph, original_super, self.context))
        return false;

    return true;
}

bool AutomorphismSearch::_componentEdgeMatch(Graph& subgraph, Graph& supergraph, int sub_idx, int super_idx, void* userdata)
{
    ComponentMatchContext& match = *(ComponentMatchContext*)userdata;
    AutomorphismSearch& self = *match.search;

    if (self.cb_edge_rank == 0)
        return true;

    const Edge& sub_edge = subgraph.getEdge(sub_idx);
    const Edge& super_edge = supergraph.getEdge(super_idx);

    int graph_sub_beg = match.sub_to_graph->at(sub_edge.beg);
    int graph_sub_end = match.sub_to_graph->at(sub_edge.end);
    int graph_super_beg = match.super_to_graph->at(super_edge.beg);
    int graph_super_end = match.super_to_graph->at(super_edge.end);

    if (graph_sub_beg < 0 || graph_sub_end < 0 || graph_super_beg < 0 || graph_super_end < 0)
        throw Error("internal: incomplete disconnected component edge mapping");

    int original_sub = self._given_graph->findEdgeIndex(self._mapping[graph_sub_beg], self._mapping[graph_sub_end]);
    int original_super = self._given_graph->findEdgeIndex(self._mapping[graph_super_beg], self._mapping[graph_super_end]);

    if (original_sub < 0 || original_super < 0)
        return false;

    return self.cb_edge_rank(*self._given_graph, original_sub, self.context) == self.cb_edge_rank(*self._given_graph, original_super, self.context);
}

bool AutomorphismSearch::_trySeedComponentSwap(const GraphDecomposer& decomposer, int component1, int component2, Array<int>& permutation,
                                               ReusableObjArray<Array<int>>* automorphisms, int max_automorphisms)
{
    if (decomposer.getComponentVerticesCount(component1) != decomposer.getComponentVerticesCount(component2) ||
        decomposer.getComponentEdgesCount(component1) != decomposer.getComponentEdgesCount(component2))
        return false;

    QS_DEF(Array<int>, component1_vertices);
    QS_DEF(Array<int>, component2_vertices);
    component1_vertices.clear();
    component2_vertices.clear();

    for (int i = _graph.vertexBegin(); i != _graph.vertexEnd(); i = _graph.vertexNext(i))
    {
        int component = decomposer.getComponent(i);
        if (component == component1)
            component1_vertices.push(i);
        if (component == component2)
            component2_vertices.push(i);
    }

    Graph component1_graph;
    Graph component2_graph;
    QS_DEF(Array<int>, graph_to_component1);
    QS_DEF(Array<int>, graph_to_component2);
    component1_graph.makeSubgraph(_graph, component1_vertices, &graph_to_component1);
    component2_graph.makeSubgraph(_graph, component2_vertices, &graph_to_component2);

    QS_DEF(Array<int>, component1_to_graph);
    QS_DEF(Array<int>, component2_to_graph);
    component1_to_graph.clear_resize(component1_graph.vertexEnd());
    component2_to_graph.clear_resize(component2_graph.vertexEnd());
    component1_to_graph.fffill();
    component2_to_graph.fffill();

    for (int i = 0; i < component1_vertices.size(); i++)
    {
        int graph_vertex = component1_vertices[i];
        int component_vertex = graph_to_component1[graph_vertex];
        if (component_vertex < 0 || !component1_graph.hasVertex(component_vertex) || component1_to_graph[component_vertex] != -1)
            throw Error("internal: invalid disconnected component subgraph mapping");
        component1_to_graph[component_vertex] = graph_vertex;
    }
    for (int i = 0; i < component2_vertices.size(); i++)
    {
        int graph_vertex = component2_vertices[i];
        int component_vertex = graph_to_component2[graph_vertex];
        if (component_vertex < 0 || !component2_graph.hasVertex(component_vertex) || component2_to_graph[component_vertex] != -1)
            throw Error("internal: invalid disconnected component supergraph mapping");
        component2_to_graph[component_vertex] = graph_vertex;
    }

    for (int i = component1_graph.vertexBegin(); i != component1_graph.vertexEnd(); i = component1_graph.vertexNext(i))
        if (component1_to_graph[i] < 0)
            throw Error("internal: incomplete disconnected component subgraph mapping");
    for (int i = component2_graph.vertexBegin(); i != component2_graph.vertexEnd(); i = component2_graph.vertexNext(i))
        if (component2_to_graph[i] < 0)
            throw Error("internal: incomplete disconnected component supergraph mapping");

    ComponentMatchContext match = {this, &component1_to_graph, &component2_to_graph};
    EmbeddingEnumerator enumerator(component2_graph);
    enumerator.setSubgraph(component1_graph);
    enumerator.userdata = &match;
    enumerator.cb_match_vertex = _componentVertexMatch;
    enumerator.cb_match_edge = _componentEdgeMatch;

    QS_DEF(Array<int>, candidate_permutation);
    candidate_permutation.clear_resize(_n);

    enumerator.processStart();
    while (enumerator.processNext())
    {
        for (int i = 0; i < _n; i++)
            candidate_permutation[i] = i;

        // The embedding maps query-local vertices to target-local vertices.
        const int* mapping = enumerator.getSubgraphMapping();
        bool valid_mapping = true;
        bool identity_mapping = true;

        for (int i = component1_graph.vertexBegin(); i != component1_graph.vertexEnd(); i = component1_graph.vertexNext(i))
        {
            int mapped = mapping[i];
            if (mapped < 0 || !component2_graph.hasVertex(mapped))
            {
                valid_mapping = false;
                break;
            }

            int graph_vertex1 = component1_to_graph[i];
            int graph_vertex2 = component2_to_graph[mapped];
            if (graph_vertex1 < 0 || graph_vertex1 >= _n || graph_vertex2 < 0 || graph_vertex2 >= _n)
                throw Error("internal: incomplete disconnected component automorphism mapping");

            if (graph_vertex1 != graph_vertex2)
                identity_mapping = false;
            candidate_permutation[graph_vertex1] = graph_vertex2;
            if (component1 != component2)
                candidate_permutation[graph_vertex2] = graph_vertex1;
        }

        if (!valid_mapping || (component1 == component2 && identity_mapping))
            continue;
        if (!_isPermutation(candidate_permutation))
        {
            if (component1 == component2)
                continue;
            throw Error("internal: disconnected component swap is not a permutation");
        }
        if (!_isAutomorphism(candidate_permutation))
            continue;

        if (automorphisms != nullptr)
        {
            automorphisms->push().copy(candidate_permutation);
            if (automorphisms->size() == 1)
                permutation.copy(candidate_permutation);
            if (automorphisms->size() >= max_automorphisms)
                return true;
            continue;
        }
        permutation.copy(candidate_permutation);
        return true;
    }

    return automorphisms != nullptr && automorphisms->size() > 0;
}



bool AutomorphismSearch::_isPermutation(const Array<int>& permutation)
{
    if (permutation.size() != _n)
        return false;

    _workperm2.zerofill();
    for (int i = 0; i < _n; i++)
    {
        int mapped = permutation[i];
        if (mapped < 0 || mapped >= _n || _workperm2[mapped] != 0)
            return false;
        _workperm2[mapped] = 1;
    }
    return true;
}


void AutomorphismSearch::_storeAutomorphism(const Array<int>& permutation, bool save_for_orbits)
{
    if (_fix.size() == worksize)
    {
        _removeGenerator(_fix.size() - 1);
        _fix.pop();
        _mcr.pop();
        _moved_vertices.pop();
    }

    int generator_index = _fix.size();
    _buildFixMcr(permutation, _fix.push(), _mcr.push(), _moved_vertices.push());
    _fixed_generator_count.push(0);
    _registerGenerator(generator_index);
    if (save_for_orbits)
    {
        _seeded_component_automorphisms.push().copy(permutation);
    }
}

void AutomorphismSearch::_registerGenerator(int generator_index)
{
    const Array<int>& moved_vertices = _moved_vertices[generator_index];
    for (int i = 0; i < moved_vertices.size(); i++)
    {
        int vertex = moved_vertices[i];
        _generators_by_vertex[vertex].push(generator_index);
        _fixed_generator_count[generator_index] += _fixedpts[vertex];
    }
}

void AutomorphismSearch::_removeGenerator(int generator_index)
{
    const Array<int>& moved_vertices = _moved_vertices[generator_index];
    for (int i = 0; i < moved_vertices.size(); i++)
    {
        Array<int>& generators = _generators_by_vertex[moved_vertices[i]];
        for (int j = 0; j < generators.size(); j++)
            if (generators[j] == generator_index)
            {
                generators.remove(j);
                break;
            }
    }
    _fixed_generator_count.pop();
}

void AutomorphismSearch::_setFixedPoint(int vertex, int value)
{
    if (_fixedpts[vertex] == value)
        return;

    int delta = value == 0 ? -1 : 1;
    Array<int>& generators = _generators_by_vertex[vertex];
    for (int i = 0; i < generators.size(); i++)
        _fixed_generator_count[generators[i]] += delta;
    _fixedpts[vertex] = value;
}
void AutomorphismSearch::_insertCellStart(int start)
{
    int position = 0;
    while (position < _cell_starts.size() && _cell_starts[position] < start)
        position++;
    if (position < _cell_starts.size() && _cell_starts[position] == start)
        return;

    _cell_starts.push(0);
    for (int i = _cell_starts.size() - 1; i > position; i--)
        _cell_starts[i] = _cell_starts[i - 1];
    _cell_starts[position] = start;
}

void AutomorphismSearch::_removeCellStart(int start)
{
    for (int position = 0; position < _cell_starts.size(); position++)
        if (_cell_starts[position] == start)
        {
            _cell_starts.remove(position);
            return;
        }
    throw Error("internal: cell start %d is missing", start);
}

void AutomorphismSearch::_seedDisconnectedComponentAutomorphisms()
{
    GraphDecomposer decomposer(_graph);
    int components_count = decomposer.decompose();
    if (components_count < 2)
        return;

    QS_DEF(Array<int>, representatives);
    QS_DEF(Array<int>, component_representatives);
    QS_DEF(Array<int>, permutation);
    ReusableObjArray<Array<int>> component_automorphisms;
    representatives.clear();
    component_representatives.clear_resize(components_count);
    component_representatives.fffill();

    for (int component = 0; component < components_count; component++)
    {
        int representative = -1;
        for (int i = 0; i < representatives.size(); i++)
        {
            int candidate = representatives[i];
            if (decomposer.getComponentVerticesCount(candidate) != decomposer.getComponentVerticesCount(component) ||
                decomposer.getComponentEdgesCount(candidate) != decomposer.getComponentEdgesCount(component))
                continue;
            if (_trySeedComponentSwap(decomposer, candidate, component, permutation))
            {
                representative = candidate;
                break;
            }
        }
        if (representative == -1)
        {
            representative = component;
            representatives.push(component);
        }
        component_representatives[component] = representative;
    }
    search_profile.components = components_count;
    search_profile.classes = representatives.size();

    for (int component = 0; component < components_count; component++)
    {

        component_automorphisms.clear();
        if (_trySeedComponentSwap(decomposer, component, component, permutation, &component_automorphisms, 4))
            for (int i = 0; i < component_automorphisms.size(); i++)
                _storeAutomorphism(component_automorphisms[i], true);
    }

    for (int component = 0; component < components_count; component++)
    {
        int representative = component_representatives[component];
        int next = -1;
        int class_size = 0;
        for (int candidate = 0; candidate < components_count; candidate++)
            if (component_representatives[candidate] == representative)
            {
                class_size++;
                if (candidate > component && next == -1)
                    next = candidate;
            }
        if (next == -1 && component != representative && class_size > 2)
            next = representative;
        if (next != -1 && _trySeedComponentSwap(decomposer, component, next, permutation))
            _storeAutomorphism(permutation, true);
    }
    search_profile.seeded = _seeded_component_automorphisms.size();
    search_profile.generators = _fix.size();

}

void AutomorphismSearch::_activateSeededComponentAutomorphisms()
{
    int pending = _seeded_component_automorphisms.size();
    search_profile.seeded = pending;
    for (int i = 0; i < pending; i++)
    {
        const Array<int>& permutation = _seeded_component_automorphisms[i];
        _joinOrbits(permutation);
        _handleAutomorphism(permutation);
    }
    _seeded_component_automorphisms.clear();
}

int AutomorphismSearch::_cmp_vertices(int idx1, int idx2, void* context)
{
    const AutomorphismSearch* self = (const AutomorphismSearch*)context;

    int degree_diff = self->_degree[idx1] - self->_degree[idx2];

    if (self->compare_vertex_degree_first)
        if (degree_diff != 0)
            return degree_diff;
    if (self->cb_vertex_cmp == 0)
        return degree_diff;

    int ret_cb_vertex_cmp = self->cb_vertex_cmp(*self->_given_graph, idx1, idx2, self->context);
    if (ret_cb_vertex_cmp != 0)
        return ret_cb_vertex_cmp;

    if (!self->compare_vertex_degree_first)
        return degree_diff;

    return 0;
}

void AutomorphismSearch::getCanonicalNumbering(Array<int>& numbering)
{
    int i;

    numbering.clear();

    for (i = 0; i < _mapping.size(); i++)
        numbering.push(_mapping[_canonlab[i]]);
}

void AutomorphismSearch::getOrbits(Array<int>& orbits) const
{
    orbits.clear_resize(_given_graph->vertexEnd());
    orbits.fffill();

    for (int i = 0; i < _mapping.size(); i++)
        orbits[_mapping[i]] = _orbits[i];
}

void AutomorphismSearch::getCanonicallyOrderedOrbits(Array<int>& orbits) const
{
    // Each vertex in the orbit has its canonical number.
    // Canonical orbit index is the minimal canonical index of the vertices
    // from this orbit
    QS_DEF(Array<int>, min_vertex_in_orbit);
    min_vertex_in_orbit.clear_resize(_given_graph->vertexEnd());
    min_vertex_in_orbit.fffill();

    for (int i = 0; i < _mapping.size(); i++)
    {
        int vertex = _canonlab[i];
        int orbit = _orbits[vertex];

        if (min_vertex_in_orbit[orbit] == -1 || min_vertex_in_orbit[orbit] > i)
            min_vertex_in_orbit[orbit] = i;
    }

    orbits.clear_resize(_given_graph->vertexEnd());
    orbits.fffill();

    for (int i = 0; i < _mapping.size(); i++)
        orbits[_mapping[i]] = min_vertex_in_orbit[_orbits[i]];
}

void AutomorphismSearch::process(Graph& graph)
{
    _prepareGraph(graph);
    search_profile = SearchProfile();
    _lab_position.clear_resize(_n);
    _cell_starts.clear();
    _ptn_change_stack.clear();
    for (int position = 0; position < _n; position++)
    {
        _lab_position[_lab[position]] = position;
        if (position == 0 || _ptn[position - 1] == 0)
            _cell_starts.push(position);
    }
    _edge_counts.clear_resize(_graph.vertexEnd());
    _edge_counts.zerofill();
    _edge_count_touched.clear();
    _edge_rank_cache.clear_resize(_graph.edgeEnd());
    _edge_rank_cache.zerofill();
    if (cb_edge_rank != 0)
        for (int edge_idx = _graph.edgeBegin(); edge_idx != _graph.edgeEnd(); edge_idx = _graph.edgeNext(edge_idx))
        {
            const Edge& edge = _graph.getEdge(edge_idx);
            int mapped_beg = _mapping[edge.beg];
            int mapped_end = _mapping[edge.end];
            int mapped_edge_idx = _given_graph->findEdgeIndex(mapped_beg, mapped_end);
            if (mapped_edge_idx == -1)
                throw Error("Internal error: edge must exists");
            _edge_rank_cache[edge_idx] = cb_edge_rank(*_given_graph, mapped_edge_idx, context);
        }

    _active.clear_resize(_n);
    _workperm.clear_resize(_n);
    _workperm2.clear_resize(_n);
    _firstlab.clear_resize(_n);
    _canonlab.clear_resize(_n);
    _fixedpts.clear_resize(_n);
    _count.clear_resize(_n);
    _orbits.clear_resize(_n);
    _fix.clear();
    _mcr.clear();
    _moved_vertices.clear();
    _seeded_component_automorphisms.clear();
    _generators_by_vertex.clear();
    for (int i = 0; i < _n; i++)
        _generators_by_vertex.push().clear();
    _fixed_generator_count.clear();
    _generator_seen_epoch.clear();
    _long_prune_candidates.clear();
    _long_prune_epoch = 0;

    if (_n == 0)
        return;

    _fixedpts.zerofill();
    _needshortprune = false;
    _orbits_num = _n;

    {
        int i, numcells = 0;

        _ptn[_n - 1] = 0;

        for (i = 0; i < _n; i++)
            if (_ptn[i] != 0)
                _ptn[i] = AUTOMORPHISM_INFINITY;
            else
                numcells++;

        _active.zerofill();

        for (i = 0; i < _n; i++)
        {
            _active[i] = 1;
            while (_ptn[i])
                i++;
        }

        for (i = 0; i < _n; ++i)
            _orbits[i] = i;

        if (getcanon)
        {
            _seedDisconnectedComponentAutomorphisms();
        }
        _Call& call = _call_stack.push();
        call.level = 1;
        call.numcells = numcells;
        call.place = _INITIAL;
    }

    int retval = -1;

    while (_call_stack.size() > 0)
    {
        search_profile.dispatches++;
        _Call call = _call_stack.top();

        if (call.place == _INITIAL)
        {
            retval = _firstNode(call.level, call.numcells);
            if (retval >= 0)
                _call_stack.pop();
        }
        else if (call.place == _FIRST_LOOP)
        {
            int tv = -1;
            if (retval != -1)
            {
                // handle the value returned from _FIRST_TO_FIRST or _FIRST_TO_OTHER
                tv = _tcells[call.level][call.k];

                if (tv == call.tv1)
                    _gca_first = call.level;

                _setFixedPoint(tv, 0);

                if (retval < call.level)
                {
                    _call_stack.pop();
                    continue; // break the _FIRST_LOOP and keep the retval;
                }

                if (_needshortprune)
                {
                    _needshortprune = false;
                    call.k = _shortPrune(_tcells[call.level], _mcr.top(), call.k);
                }
                _recover(call.level);
                // advance the _FIRST_LOOP counter
                call.k++;
            }

            for (; call.k < _tcells[call.level].size(); call.k++)
            {
                tv = _tcells[call.level][call.k];

                if (_orbits[tv] == tv) // not equivalent to the previous child?
                    break;
            }

            if (call.k == _tcells[call.level].size())
            {
                // return from _FIRST_LOOP
                retval = call.level - 1;
                _call_stack.pop();
                continue;
            }

            _call_stack.top() = call;

            _breakout(call.level + 1, call.tc, tv);
            _cosetindex = tv;
            _setFixedPoint(tv, 1);

            _Call& newcall = _call_stack.push();
            newcall.level = call.level + 1;
            newcall.numcells = call.numcells + 1;
            if (tv == call.tv1)
                newcall.place = _FIRST_TO_FIRST;
            else
                newcall.place = _FIRST_TO_OTHER;

            // discard the old return value
            retval = -1;
        }
        else if (call.place == _FIRST_TO_FIRST)
        {
            retval = _firstNode(call.level, call.numcells);

            if (retval >= 0)
                // _FIRST_LOOP did not happen; pass the return value to the caller
                _call_stack.pop();
        }
        else if (call.place == _FIRST_TO_OTHER || call.place == _OTHER_TO_OTHER)
        {
            retval = _otherNode(call.level, call.numcells);

            if (retval >= 0)
                // _OTHER_LOOP did not happen; pass the return value to the caller
                _call_stack.pop();
        }
        else if (call.place == _OTHER_LOOP)
        {
            int tv;

            if (retval != -1)
            {
                // handle the value returned from _OTHER_TO_OTHER
                tv = _tcells[call.level][call.k];

                _setFixedPoint(tv, 0);

                if (retval < call.level)
                {
                    _call_stack.pop();
                    continue; // break the _OTHER_LOOP and keep the retval;
                }

                // use stored automorphism data to prune target cell
                if (_needshortprune)
                {
                    _needshortprune = false;
                    call.k = _shortPrune(_tcells[call.level], _mcr.top(), call.k);
                }

                if (tv == call.tv1)
                    call.k = _longPrune(_tcells[call.level], call.k);

                _recover(call.level);
                // advance the _OTHER_LOOP counter
                call.k++;
            }

            if (call.k == _tcells[call.level].size())
            {
                // return from _OTHER_LOOP
                retval = call.level - 1;
                _call_stack.pop();
                continue;
            }

            _call_stack.top() = call;

            tv = _tcells[call.level][call.k];

            _breakout(call.level + 1, call.tc, tv);
            _setFixedPoint(tv, 1);

            _Call& newcall = _call_stack.push();
            newcall.level = call.level + 1;
            newcall.numcells = call.numcells + 1;
            newcall.place = _OTHER_TO_OTHER;

            // discard the old return value
            retval = -1;
        }
        else
            throw Error("internal: bad command %d", call.place);
    }
}

int AutomorphismSearch::_firstNode(int level, int numcells)
{
    _refine(level, numcells);

    _tcells.resize(level + 1);

    if (numcells == _n) // found first leaf?
    {
        _gca_first = level;

        _firstlab.copy(_lab);

        if (getcanon)
        {
            _canonlevel = _gca_canon = level;
            _canonlab.copy(_lab);
            _activateSeededComponentAutomorphisms();
        }

        return level - 1;
    }

    // locate new target cell
    int tc = _targetcell(level, _tcells[level]);
    int tv1 = _tcells[level][0];

    _call_stack.pop();

    // use the elements of the target cell to produce the children
    _Call& call = _call_stack.push();
    call.level = level;
    call.k = 0;
    call.tc = tc;
    call.tv1 = tv1;
    call.numcells = numcells;
    call.place = _FIRST_LOOP;

    return -1;
}

int AutomorphismSearch::_otherNode(int level, int numcells)
{
    _refine(level, numcells);
    _tcells.resize(level + 1);
    int rtnlevel = _processNode(level, numcells);

    if (rtnlevel < level) // keep returning if necessary
        return rtnlevel;

    int tc = _targetcell(level, _tcells[level]);

    if (_needshortprune)
    {
        _needshortprune = false;
        _shortPrune(_tcells[level], _mcr.top(), 0);
    }

    int tv1 = _tcells[level][0];

    _call_stack.pop();

    // use the elements of the target cell to produce the children
    _Call& call = _call_stack.push();
    call.level = level;
    call.k = 0;
    call.tc = tc;
    call.tv1 = tv1;
    call.numcells = numcells;
    call.place = _OTHER_LOOP;

    return -1;
}

void AutomorphismSearch::_recover(int level)
{

    while (_ptn_change_stack.size() > 0)
    {
        int boundary = _ptn_change_stack.top();
        if (_ptn[boundary] <= level)
            break;
        _ptn[boundary] = AUTOMORPHISM_INFINITY;
        _removeCellStart(boundary + 1);
        _ptn_change_stack.pop();
    }

    if (getcanon)
    {
        if (level < _gca_canon)
            _gca_canon = level;

        if (level < _gca_first)
            throw Error("internal error?");
    }
}

void AutomorphismSearch::_breakout(int level, int tc, int tv)
{
    _active.zerofill();

    _active[tc] = 1;

    int i = tc;
    int prev = tv;

    do
    {
        int next = _lab[i];

        _lab[i++] = prev;
        _lab_position[prev] = i - 1;
        prev = next;
    } while (prev != tv);
    _ptn[tc] = level;
    _ptn_change_stack.push(tc);
    _insertCellStart(tc + 1);
}

int AutomorphismSearch::_shortPrune(Array<int>& tcell, Array<int>& mcr, int idx)
{
    int i, j;
    int ret = idx;

    for (i = j = 0; i < tcell.size(); i++)
        if (mcr[tcell[i]])
            tcell[j++] = tcell[i];
        else if (idx >= i)
            ret--;

    tcell.resize(j);
    return ret;
}

int AutomorphismSearch::_longPrune(Array<int>& tcell, int idx)
{
    auto profile_start = ProfileClock::now();
    search_profile.prune_calls++;
    int i, j, k;
    int ret = idx;

    if (_generator_seen_epoch.size() < _fix.size())
    {
        int old_size = _generator_seen_epoch.size();
        _generator_seen_epoch.resize(_fix.size());
        for (int generator = old_size; generator < _generator_seen_epoch.size(); generator++)
            _generator_seen_epoch[generator] = 0;
    }
    if (_long_prune_epoch == std::numeric_limits<int>::max())
    {
        _generator_seen_epoch.zerofill();
        _long_prune_epoch = 0;
    }
    int epoch = ++_long_prune_epoch;
    _long_prune_candidates.clear();
    for (int position = 0; position < tcell.size(); position++)
    {
        const Array<int>& generators = _generators_by_vertex[tcell[position]];
        for (int generator_index = 0; generator_index < generators.size(); generator_index++)
        {
            int generator = generators[generator_index];
            if (_generator_seen_epoch[generator] != epoch)
            {
                _generator_seen_epoch[generator] = epoch;
                _long_prune_candidates.push(generator);
            }
        }
    }

    for (int candidate_index = 0; candidate_index < _long_prune_candidates.size(); candidate_index++)
    {
        k = _long_prune_candidates[candidate_index];
        if (_fixed_generator_count[k] != 0)
            continue;

        bool affects_target_cell = false;
        for (int position = 0; position < tcell.size(); position++)
            if (!_mcr[k][tcell[position]])
            {
                affects_target_cell = true;
                break;
            }
        if (!affects_target_cell)
            continue;
        for (i = j = 0; i < tcell.size(); i++)
            if (_mcr[k][tcell[i]])
                tcell[j++] = tcell[i];
            else if (idx >= i)
                ret--;

        tcell.resize(j);
        idx = ret;
    }

    search_profile.prune_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(ProfileClock::now() - profile_start).count();
    return ret;
}

int AutomorphismSearch::_processNode(int level, int numcells)
{
    int i;

    // no idea what this nauty's if() means.
    // if (_eqlev_first != level && (!getcanon || _comp_canon < 0))
    //   code = 4;

    if (numcells != _n) // discrete partition?
        return level;


    if (_cancellation_handler != nullptr && _cancellation_handler->isCancelled())
    {
        std::fprintf(stderr,
                     "PROFILE n=%d components=%d classes=%d generators=%d seeded=%d dispatch=%llu refine=%llu/%llu ms target=%llu/%llu ms prune=%llu/%llu ms compare=%llu/%llu ms checks=%llu\\n",
                     _n, search_profile.components, search_profile.classes, search_profile.generators, search_profile.seeded,
                     (unsigned long long)search_profile.dispatches, (unsigned long long)search_profile.refine_calls,
                     (unsigned long long)(search_profile.refine_ns / 1000000),
                     (unsigned long long)search_profile.target_calls, (unsigned long long)(search_profile.target_ns / 1000000),
                     (unsigned long long)search_profile.prune_calls, (unsigned long long)(search_profile.prune_ns / 1000000),
                     (unsigned long long)search_profile.compare_calls, (unsigned long long)(search_profile.compare_ns / 1000000),
                     (unsigned long long)search_profile.automorphism_checks);
        throw TimeoutException("%s", _cancellation_handler->cancelledRequestMessage());
    }

    for (i = 0; i < _n; i++)
        _workperm[_firstlab[i]] = _lab[i];

    if (_isAutomorphism(_workperm))
    {
        // _lab is equivalent to firstlab
        _storeAutomorphism(_workperm, false);
        _joinOrbits(_workperm);
        _handleAutomorphism(_workperm);

        return _gca_first;
    }

    if (getcanon)
    {
        /*if (_comp_canon == 0)
        {
           // again, strange nauty's if()
           if (level < _canonlevel)
              _comp_canon = 1;
           else
              _comp_canon = _compareCanon();
        }*/

        int comp_canon = _compareCanon();

        if (comp_canon == 0)
        {
            // _lab is equivalent to canonlab
            for (i = 0; i < _n; i++)
                _workperm[_canonlab[i]] = _lab[i];

            _storeAutomorphism(_workperm, false);

            int norb = _orbits_num;

            _joinOrbits(_workperm);

            if (norb != _orbits_num)
            {
                _handleAutomorphism(_workperm);
                if (_orbits[_cosetindex] < _cosetindex)
                    return _gca_first;
            }
            if (_gca_canon != _gca_first)
                _needshortprune = true;
            return _gca_canon;
        }
        else if (comp_canon > 0)
        {
            // _lab is better than canonlab
            _canonlab.copy(_lab);
            _canonlevel = _gca_canon = level;
        }
    }

    return level - 1;
}

void AutomorphismSearch::_joinOrbits(const Array<int>& perm)
{
    int i, j1, j2;

    for (i = 0; i < _n; i++)
    {
        j1 = _orbits[i];

        while (_orbits[j1] != j1)
            j1 = _orbits[j1];

        j2 = _orbits[perm[i]];

        while (_orbits[j2] != j2)
            j2 = _orbits[j2];

        if (j1 < j2)
            _orbits[j2] = j1;
        else if (j1 > j2)
            _orbits[j1] = j2;
    }

    _orbits_num = 0;

    for (i = 0; i < _n; i++)
    {
        _orbits[i] = _orbits[_orbits[i]];
        if (_orbits[i] == i)
            _orbits_num++;
    }
}

bool AutomorphismSearch::_isAutomorphism(Array<int>& perm)
{
    search_profile.automorphism_checks++;
    for (int i = _graph.edgeBegin(); i != _graph.edgeEnd(); i = _graph.edgeNext(i))
    {
        const Edge& edge = _graph.getEdge(i);

        if (!_graph.haveEdge(perm[edge.beg], perm[edge.end]))
            return false;
    }

    if (cb_check_automorphism != 0)
    {
        QS_DEF(Array<int>, perm_mapping);

        perm_mapping.clear_resize(_given_graph->vertexEnd());
        perm_mapping.fffill();

        for (int i = 0; i < _n; i++)
            perm_mapping[_mapping[i]] = _mapping[perm[i]];

        return cb_check_automorphism(*_given_graph, perm_mapping, context);
    }

    return true;
}

// lab vs. canonlab
int AutomorphismSearch::_compareCanon()
{
    int i;
    auto profile_start = ProfileClock::now();
    search_profile.compare_calls++;

    QS_DEF(Array<int>, map);
    QS_DEF(Array<int>, canon_map);

    map.clear_resize(_n);
    canon_map.clear_resize(_n);

    for (i = 0; i < _n; i++)
    {
        map[i] = _mapping[_lab[i]];
        canon_map[i] = _mapping[_canonlab[i]];
    }

    if (cb_compare_mapped == 0)
        throw Error("cb_compare_mapped = 0");
    int result = cb_compare_mapped(*_given_graph, map, canon_map, context);
    search_profile.compare_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(ProfileClock::now() - profile_start).count();
    return result;
}

void AutomorphismSearch::_buildFixMcr(const Array<int>& perm, Array<int>& fix, Array<int>& mcr, Array<int>& moved_vertices)
{
    int i;

    fix.clear_resize(_n);
    mcr.clear_resize(_n);
    moved_vertices.clear();
    fix.zerofill();
    mcr.zerofill();

    _workperm2.zerofill();

    for (i = 0; i < _n; ++i)
    {
        if (perm[i] != i)
            moved_vertices.push(i);

        if (perm[i] == i)
        {
            fix[i] = 1;
            mcr[i] = 1;
        }
        else if (_workperm2[i] == 0)
        {
            int l = i;

            do
            {
                _workperm2[l] = 1;
                l = perm[l];
            } while (l != i);

            mcr[i] = 1;
        }
    }
}

int AutomorphismSearch::_targetcell(int /*level*/, Array<int>& cell)
{
    int i = 0, j, k;
    auto profile_start = ProfileClock::now();
    search_profile.target_calls++;
    int ibest = -1, jbest = -1, bestdegree = -1;
    for (int cell_index = 0; cell_index < _cell_starts.size(); cell_index++)
    {
        i = _cell_starts[cell_index];
        j = cell_index + 1 < _cell_starts.size() ? _cell_starts[cell_index + 1] - 1 : _n - 1;
        if (i == j)
            continue;

        int degree = _degree[_mapping[_lab[i]]];

        // Prefer zero-degree cells, then the smallest non-singleton cell.
        int cell_size = j - i;
        int best_size = jbest - ibest;
        bool prefer = ibest == -1 || (degree == 0 && bestdegree != 0) ||
                      ((degree == 0) == (bestdegree == 0) && (ibest == -1 || cell_size < best_size));
        if (prefer)
        {
            jbest = j;
            ibest = i;
            bestdegree = degree;
        }
    }


    if (ibest == -1)
        throw Error("(intenal error) target cell cannot be found");

    i = ibest;
    j = jbest;

    cell.clear();

    int imin = 0;

    for (k = i; k <= j; k++)
    {
        cell.push(_lab[k]);

        if (cell.size() > 0 && cell[cell.size() - 1] < cell[imin])
            imin = cell.size() - 1;
    }
    if (imin > 0)
        cell.swap(0, imin);

    search_profile.target_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(ProfileClock::now() - profile_start).count();
    return i;
}

void AutomorphismSearch::_refine(int level, int& numcells)
{
    auto profile_start = ProfileClock::now();
    search_profile.refine_calls++;
    if (refine_by_sorted_neighbourhood)
        _refineBySortingNeighbourhood(level, numcells);
    else
        _refineOriginal(level, numcells);
    search_profile.refine_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(ProfileClock::now() - profile_start).count();
}

void AutomorphismSearch::_refineOriginal(int level, int& numcells)
{
    int hint = 0;

    int split1 = -1;

    while (numcells < _n)
    {
        int split2;

        if (_active[hint])
            split1 = hint;
        else
        {
            int i;
            for (i = 0; i < _n; i++)
            {
                split1 = (split1 + 1) % _n;
                if (_active[split1])
                    break;
            }
            if (i == _n)
                break;
        }

        _active[split1] = 0;

        for (split2 = split1; _ptn[split2] > level; split2++)
            ;

        _edge_ranks_in_refine.clear();

        _refineByCell(split1, split2, level, numcells, hint, -1);
        // Check if there are exists edge with different ranks
        // Last element in _edge_ranks_in_refine array contains positive value that
        // means that such edge rank exists. But it is nessesary to refine by all ranks
        // except one because cells have already been refined by edges without ranks.
        for (int i = 0; i < _edge_ranks_in_refine.size() - 1; i++)
            if (_edge_ranks_in_refine[i] != 0)
                _refineByCell(split1, split2, level, numcells, hint, i);
    }
}

void AutomorphismSearch::_refineBySortingNeighbourhood(int level, int& numcells)
{
    // This refine procedure works like refining by sorting neighbourhood ranks

    while (true)
    {
        // Collect active cells
        _work_active_cells.clear();
        for (int i = 0; i < _n; i++)
        {
            int split1;
            if (_active[i])
            {
                split1 = i;

                int split2;
                for (split2 = split1; _ptn[split2] > level; split2++)
                    ;

                int(&split_cell)[2] = _work_active_cells.push();
                split_cell[0] = split1;
                split_cell[1] = split2;

                _active[i] = 0;
            }
        }

        if (_work_active_cells.size() == 0)
            break;

        // Refine all cells by collected active cells
        for (int i = 0; i < _work_active_cells.size(); i++)
        {
            int(&split_cell)[2] = _work_active_cells[i];

            int split1 = split_cell[0], split2 = split_cell[1];

            int dummy_hint;
            _refineByCell(split1, split2, level, numcells, dummy_hint, -1);

            if (numcells == _n)
                return;
        }
    }
}


void AutomorphismSearch::_refineByCell(int split1, int split2, int level, int& numcells, int& hint, int target_edge_rank)
{
    int i, j;
    Array<int>& edge_counts = _edge_counts;
    _edge_count_touched.clear();
    for (j = split1; j <= split2; j++)
    {
        int splitter_vertex = _lab[j];
        const Vertex& splitter = _graph.getVertex(splitter_vertex);
        for (int nei = splitter.neiBegin(); nei != splitter.neiEnd(); nei = splitter.neiNext(nei))
        {
            int neighbor = splitter.neiVertex(nei);
            if (cb_edge_rank != 0)
            {
                int edge_rank = _edge_rank_cache[splitter.neiEdge(nei)];
                if (target_edge_rank == -1)
                {
                    while (_edge_ranks_in_refine.size() <= edge_rank)
                        _edge_ranks_in_refine.push(0);
                    _edge_ranks_in_refine[edge_rank]++;
                }
                else if (target_edge_rank != edge_rank)
                    continue;
            }
            if (edge_counts[neighbor]++ == 0)
                _edge_count_touched.push(neighbor);
        }
    }

    if (split1 == split2) // trivial splitting cell
    {
        int cell1, cell2;
        QS_DEF(Array<int>, affected_cell_starts);
        affected_cell_starts.clear();

        for (int touched = 0; touched < _edge_count_touched.size(); touched++)
        {
            int position = _lab_position[_edge_count_touched[touched]];
            int low = 0, high = _cell_starts.size();
            while (low < high)
            {
                int middle = low + (high - low) / 2;
                if (_cell_starts[middle] <= position)
                    low = middle + 1;
                else
                    high = middle;
            }
            int start = _cell_starts[low - 1];
            if (_ptn[start] <= level)
                continue;
            bool already_added = false;
            for (int i = 0; i < affected_cell_starts.size(); i++)
                if (affected_cell_starts[i] == start)
                {
                    already_added = true;
                    break;
                }
            if (!already_added)
                affected_cell_starts.push(start);
        }

        for (int i = 1; i < affected_cell_starts.size(); i++)
        {
            int start = affected_cell_starts[i];
            int j = i;
            while (j > 0 && affected_cell_starts[j - 1] > start)
            {
                affected_cell_starts[j] = affected_cell_starts[j - 1];
                j--;
            }
            affected_cell_starts[j] = start;
        }

        for (int affected = 0; affected < affected_cell_starts.size(); affected++)
        {
            cell1 = affected_cell_starts[affected];
            for (cell2 = cell1; _ptn[cell2] > level; cell2++)
                ;

            int c1 = cell1, c2 = cell2;
            while (c1 <= c2)
            {
                if (edge_counts[_lab[c1]] != 0)
                    c1++;
                else
                {
                    std::swap(_lab[c1], _lab[c2]);
                    _lab_position[_lab[c1]] = c1;
                    _lab_position[_lab[c2]] = c2;
                    c2--;
                }
            }

            if (c2 >= cell1 && c1 <= cell2)
            {
                if (_ptn[c2] > level)
                {
                    _ptn[c2] = level;
                    _ptn_change_stack.push(c2);
                    _insertCellStart(c2 + 1);
                }
                if (c1 <= cell2)
                    _insertCellStart(c1);
                numcells++;

                if (_active[cell1] || (c2 - cell1 >= cell2 - c1 && !refine_by_sorted_neighbourhood))
                {
                    _active[c1] = 1;
                    if (c1 == cell2)
                        hint = c1;
                }
                else
                {
                    _active[cell1] = 1;
                    if (c2 == cell1)
                        hint = cell1;
                }
            }
        }
    }
    else // nontrivial splitting cell
    {
        int cell1, cell2;

        for (cell1 = 0; cell1 < _n; cell1 = cell2 + 1)
        {
            for (cell2 = cell1; _ptn[cell2] > level; ++cell2)
                ;
            if (cell1 == cell2)
                continue;

            int bmin = _n;

            _bucket.clear();

            for (i = cell1; i <= cell2; i++)
            {
                int cnt = edge_counts[_lab[i]];
                while (_bucket.size() <= cnt)
                    _bucket.push(0);

                _bucket[cnt]++;

                if (cnt < bmin)
                    bmin = cnt;

                _count[i] = cnt;
            }

            if (bmin == _bucket.size() - 1)
                continue;

            if (refine_reverse_degree)
            {
                for (i = cell1; i <= cell2; i++)
                {
                    _count[i] = _bucket.size() - _count[i] - 1;
                }
                for (i = _bucket.size() - 1; i >= bmin; i--)
                {
                    int dest = _bucket.size() - i - 1;
                    if (dest < i)
                        std::swap(_bucket[i], _bucket[dest]);
                }
                _bucket.resize(_bucket.size() - bmin);
                bmin = 0;
            }

            int c1 = cell1, c2;
            int maxcell = -1, maxpos = -1;
            int last_c1 = -1;

            for (i = bmin; i < _bucket.size(); i++)
            {
                if (_bucket[i] == 0)
                    continue;

                int subcell_size = _bucket[i];
                c2 = c1 + subcell_size;
                _bucket[i] = c1;
                last_c1 = c1;

                if (c2 - c1 > maxcell)
                {
                    maxcell = c2 - c1;
                    maxpos = c1;
                }

                if (c1 != cell1)
                {
                    _active[c1] = 1;
                    if (c2 - c1 == 1)
                        hint = c1;
                    numcells++;
                }
                if (c2 <= cell2 && _ptn[c2 - 1] > level)
                {
                    _ptn[c2 - 1] = level;
                    _ptn_change_stack.push(c2 - 1);
                    _insertCellStart(c2);
                }

                c1 = c2;
            }

            for (i = cell1; i <= cell2; i++)
                _workperm2[_bucket[_count[i]]++] = _lab[i];

            for (i = cell1; i <= cell2; i++)
            {
                _lab[i] = _workperm2[i];
                _lab_position[_lab[i]] = i;
            }

            if (_active[cell1] == 0)
            {
                _active[cell1] = 1;

                // When sorting by neighbourhood is is allowed only to exclude
                // the last created subcell. For ordinary refine greatest cell is excluded.
                if (!refine_by_sorted_neighbourhood)
                    _active[maxpos] = 0;
                else
                    _active[last_c1] = 0;
            }
        }
    }
    for (i = 0; i < _edge_count_touched.size(); i++)
        edge_counts[_edge_count_touched[i]] = 0;
}

void AutomorphismSearch::_handleAutomorphism(const Array<int>& perm)
{
    if (cb_automorphism != 0)
    {
        QS_DEF(Array<int>, perm2);
        int i;

        perm2.clear_resize(_given_graph->vertexEnd());
        perm2.fffill();

        for (i = 0; i < _n; i++)
            perm2[_mapping[i]] = _mapping[perm[i]];

        cb_automorphism(perm2.ptr(), context_automorphism);
    }
}
