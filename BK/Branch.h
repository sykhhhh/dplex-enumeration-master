#pragma once
#ifndef BRANCH_H
#define BRANCH_H

#include "Graph.h"
double fuck_time;
/**
 * @brief partial set in BK algorithm
 * make sure after pushing v, it is still a DPlex
 */
class Plex
{
public:
    Graph *g;
    vector<ui> *vertex_index_map; // (*vertex_index_map)[v] is the original index of v in G
    ui sz;
    vector<ui> vertices_id;
    vector<ui> vertices_original_id;
    vector<ui> cnt_neighbor_in;  // |S \cap N_G^{in}|
    vector<ui> cnt_neighbor_out; // |S \cap N_G^{out}|
    int k, l;
    Plex(Graph &_g, vector<ui> &v_map, vector<ui> &S,
         int paramK, int paramL) : sz(0), vertex_index_map(&v_map), g(&_g), k(paramK), l(paramL)
    {
        ui n = g->n;
        assert(vertex_index_map->size() == n);
        cnt_neighbor_in.resize(n);
        cnt_neighbor_out.resize(n);
        vertices_id = S;
        vertices_original_id.resize(S.size());
        for (ui i = 0; i < S.size(); i++)
        {
            ui v = S[i];
            assert(v < n);
            vertices_original_id[i] = v_map[v];
            for (ui u : g->vertices[v].out_neighbors)
            {
                assert(u < n);
                cnt_neighbor_in[u]++;
            }
            for (ui u : g->vertices[v].in_neighbors)
            {
                assert(u < n);
                cnt_neighbor_out[u]++;
            }
        }
        sz = S.size();
    }
    ~Plex() {}
    Plex &operator=(const Plex &other)
    {
        g = other.g;
        vertex_index_map = other.vertex_index_map;
        sz = other.sz;
        vertices_id = other.vertices_id;
        vertices_original_id = other.vertices_original_id;
        cnt_neighbor_in = other.cnt_neighbor_in;
        cnt_neighbor_out = other.cnt_neighbor_out;
        k = other.k;
        l = other.l;
        return *this;
    }
    inline void push_vertex(ui v, vector<ui> &C, vector<ui> &X)
    {
        assert(v < vertex_index_map->size());
        assert(sz == vertices_id.size());
        vertices_id.push_back(v);
        vertices_original_id.push_back((*vertex_index_map)[v]);
        sz++;
        assert(v < g->n);
        // update info
        Timer t;
        for (ui w : g->vertices[v].out_neighbors)
        {
            assert(w < g->n);
            cnt_neighbor_in[w]++;
        }
        for (ui w : g->vertices[v].in_neighbors)
        {
            assert(w < g->n);
            cnt_neighbor_out[w]++;
        }
        fuck_time+=t.get_time();
        // update sets: C and X
        // reduce C
        {
            ui size_C = 0;
            for (ui i = 0; i < C.size(); i++)
            {
                ui w = C[i];
                if (cnt_neighbor_in[w] + l < sz + 1 || cnt_neighbor_out[w] + k < sz + 1)
                    continue;
                C[size_C++] = w;
            }
            C.resize(size_C);
            sort(C.begin(), C.end());
        }
        // reduce X
        {
            ui size_X = 0;
            for (ui i = 0; i < X.size(); i++)
            {
                ui w = X[i];
                if (cnt_neighbor_in[w] + l < sz + 1 || cnt_neighbor_out[w] + k < sz + 1)
                    continue;
                X[size_X++] = w;
            }
            X.resize(size_X);
            sort(X.begin(), X.end());
        }
        // reduce C and X
        for (ui u : vertices_id)
        {
            if (cnt_neighbor_in[u] + l == vertices_id.size())
            {
                if (u == v || !g->exist_edge(v, u)) // remove all non-in-neighbors of u
                {
                    auto copy = C;
                    C.resize(0);
                    intersect_set(copy, g->vertices[u].in_neighbors, C);
                    copy = X;
                    X.resize(0);
                    intersect_set(copy, g->vertices[u].in_neighbors, X);
                }
            }
            if (cnt_neighbor_out[u] + k == vertices_id.size()) // remove all non-out-neighbors of u
            {
                if (u == v || !g->exist_edge(u, v))
                {
                    auto copy = C;
                    C.resize(0);
                    intersect_set(copy, g->vertices[u].out_neighbors, C);
                    copy = X;
                    X.resize(0);
                    intersect_set(copy, g->vertices[u].out_neighbors, X);
                }
            }
        }
    }
    inline ui size() const
    {
        return sz;
    }
    bool is_DPlex()
    {
        unordered_map<int, int> din, dout;
        auto S = vertices_id;
        unique_vector(S);
        for (ui u : S)
        {
            for (ui v : S)
            {
                if (v >= u)
                    break;
                if (g->exist_edge(u, v))
                    dout[u]++, din[v]++;
                if (g->exist_edge(v, u))
                    dout[v]++, din[u]++;
            }
        }
        for (ui u : S)
        {
            if (dout[u] + k < S.size() || din[u] + l < S.size())
                return false;
        }
        return true;
    }
};

/**
 * @brief conduct recursive BK algorithm
 */
class Branch
{
    Graph &g;
    Output &out;
    int k, l;
    int max_plex_size;

public:
    ll dfs_cnt;
    double bk_time;
    double copy_time;
    double reduce_time;
    double induce_time;
    map<int, int> cnt;
    Branch(Graph &_g, Output &_o, int paramK, int paramL) : g(_g), out(_o), k(paramK), l(paramL),
                                                            bk_time(0), copy_time(0), reduce_time(0),
                                                            induce_time(0)
    {
        dfs_cnt = 0;
    }
    ~Branch() {}

    void print_results()
    {
        print_module_time("BK", bk_time);
        print_module_time("copy", copy_time);
        print_module_time("reduceCX", reduce_time);
        print_module_time("induce", induce_time);
        puts("");
    }

    void run()
    {
        ui n = g.n;
        vector<ui> C(n);
        for (ui i = 0; i < n; i++)
            C[i] = i;
        vector<ui> X;
        vector<ui> S;
        generate_set_for_BK_invocation(S, C, X);
        print_results();
    }

    /**
     * @brief remove u in C\cup X if S+u is not a DPlex; then invoke BK
     */
    void generate_set_for_BK_invocation(vector<ui> &S, vector<ui> &C, vector<ui> &X)
    {
        // invoke BK
        // step-1: induce subgraph
        Timer induce_timer;
        vector<int> vis(g.n, -1);
        Graph subg(S, C, X, vis, g);
        vector<ui> vertex_map(subg.n);
        for (ui i = 0; i < subg.n; i++)
            vertex_map[i] = subg.vertices[i].origin_id;
        Plex S_plex(subg, vertex_map, S, k, l);
        induce_time += induce_timer.get_time();

        Timer t;
        bk(S_plex, C, X);
        bk_time += t.get_time();
    }

    /**
     * @brief recursive BK algorithm
     * @param S partial set, which is always a DPlex
     * @param C candidate set; for u in C, S+u is a DPlex
     * @param X excluding set; for u in X, S+u is a DPlex
     * we use vertices in C to extend S, and use X to judge whether S is maximal
     */
    void bk(Plex &S, vector<ui> &C, vector<ui> &X)
    {
        dfs_cnt++;
        // if(dfs_cnt%1000000==0) cerr<<out.counter<<endl;
        if (!C.size())
        {
            if (!X.size()) // is maximal
            {
                out.dump(S.vertices_original_id, S.sz);
            }
            return;
        }
        else if (C.size() == 1)
        {
            ui u = C[0];
            C.pop_back();
            S.push_vertex(u, C, X);
            bk(S, C, X);
            return;
        }

        ui total_size = S.sz + C.size();
        if (total_size < k)
        {
            return;
        }
        while (C.size())
        {
            ui u = *C.rbegin();
            C.pop_back();
            // generate a new branch: must include u
            {
                Timer t;
                auto S_new = S;
                auto C_new = C;
                auto X_new = X;
                copy_time += t.get_time();
                Timer tt;
                S_new.push_vertex(u, C_new, X_new);
                reduce_time += tt.get_time();
                // enter the recursive sub-branch
                bk(S_new, C_new, X_new);
            }
            // exclude u
            X.push_back(u);
        }
    }
};

#endif