#ifndef GRAPH_H
#define GRAPH_H

#include "LinearHeap.h"

class Vertex
{
public:
    ui id;
    ui origin_id;
    ui din;
    ui dout;
    vector<ui> out_neighbors;
    vector<ui> in_neighbors;
    Vertex() : din(0), dout(0)
    {
    }
    ~Vertex()
    {
    }
    Vertex &operator=(const Vertex &other)
    {
        id = other.id;
        origin_id = other.origin_id;
        din = other.din;
        dout = other.dout;
        out_neighbors = other.out_neighbors;
        in_neighbors = other.in_neighbors;
        return *this;
    }
    bool operator==(const Vertex &other) const
    {
        return id == other.id;
    }
    bool operator<(const Vertex &other) const
    {
        return id < other.id;
    }
};

class Graph
{
public:
    ui n, m;
    vector<Vertex> vertices;
    Graph() : n(0), m(0)
    {
    }
    /**
     * @brief induce the subgraph G[S]
     * @param vis tool-array: make sure vis[0...n-1]=-1 before and after this function
     */
    Graph(vector<ui> &S, vector<int> &vis, const Graph &origin_g)
    {
        sort(S.begin(), S.end());
        n = S.size();
        assert(vis.size() >= n);
        // next, vis[u] is the new id of u in subgraph
        for (ui i = 0; i < n; i++)
        {
            assert(vis[S[i]] == -1);
            assert(S[i] < vis.size());
            vis[S[i]] = i;
        }
        vertices.resize(n);
        for (ui u_origin : S)
        {
            ui u = vis[u_origin];
            auto &vertex = vertices[u];
            vertex.id = u;
            auto &origin_vertex = origin_g.vertices[u_origin];
            vertex.origin_id = origin_vertex.origin_id;
            for (ui v_origin : origin_vertex.out_neighbors)
            {
                int v = vis[v_origin];
                if (v == -1)
                    continue;
                assert(v < S.size());
                assert(u != v);
                vertex.out_neighbors.push_back(v);
                vertices[v].in_neighbors.push_back(u);
            }
        }
        m = 0;
        for (auto &vertex : vertices)
        {
            vertex.in_neighbors.shrink_to_fit();
            vertex.out_neighbors.shrink_to_fit();
            vertex.din = vertex.in_neighbors.size();
            vertex.dout = vertex.out_neighbors.size();
            m += vertex.din;
            for (ui u : vertex.out_neighbors)
                assert(u != vertex.id);
            for (ui u : vertex.in_neighbors)
                assert(u != vertex.id);
        }
        // clear the tool-array
        for (ui i = 0; i < n; i++)
        {
            vis[S[i]] = -1;
        }
    }
    /**
     * @brief induce the subgraph G[S + X]; and ignore the edges between X
     * @param vis tool-array: make sure vis[0...n-1]=-1 before and after this function
     *
     * @return each vertex in S,C,X will be re-mapped to the new index of subgraph
     */
    Graph(vector<ui> &S, vector<ui> &C, vector<ui> &X, vector<int> &vis, const Graph &origin_g)
    {
        vector<ui> V = C;
        V.insert(V.end(), S.begin(), S.end());
        n = V.size() + X.size();
        unique_vector(V);
        assert( n== V.size() + X.size());
        for(ui u:X)
            assert(!has(V, u));
        vertices.resize(n);
        ui idx = 0;
        // re-map vertices to {0,...,n-1}
        for (ui v : V)
        {
            assert(vis[v] == -1);
            vis[v] = idx++;
        }
        for (ui v : X)
        {
            assert(vis[v] == -1);
            vis[v] = idx++;
        }
        // out-neighbors of S+C
        for (ui v_origin : V)
        {
            ui v = vis[v_origin];
            auto &vertex = vertices[v];
            auto &origin_vertex = origin_g.vertices[v_origin];
            vertex.id = v;
            vertex.origin_id = origin_vertex.origin_id;
            for (ui u_origin : origin_vertex.out_neighbors)
            {
                ui u = vis[u_origin];
                if (u == -1)
                    continue;
                assert(u < idx);
                vertex.out_neighbors.push_back(u);
                vertices[u].in_neighbors.push_back(v);
            }
        }
        // out-neighbors of X in S+C
        for (ui v_origin : X)
        {
            ui v = vis[v_origin];
            assert(v >= V.size());
            auto &vertex = vertices[v];
            auto &origin_vertex = origin_g.vertices[v_origin];
            vertex.id = v;
            vertex.origin_id = origin_vertex.origin_id;
            for (ui u_origin : origin_vertex.out_neighbors)
            {
                ui u = vis[u_origin];
                if (u == -1)
                    continue;
                if (u >= V.size()) // u in X, then we ignore (v,u)
                    continue;
                vertex.out_neighbors.push_back(u);
                vertices[u].in_neighbors.push_back(v);
            }
        }
        m = 0;
        for (auto &vertex : vertices)
        {
            sort(vertex.out_neighbors.begin(), vertex.out_neighbors.end());
            sort(vertex.in_neighbors.begin(), vertex.in_neighbors.end());
            vertex.din = vertex.in_neighbors.size();
            vertex.dout = vertex.out_neighbors.size();
            m += vertex.din;
        }
        // re-map and clear the tool-array
        for (ui &v : S)
        {
            ui new_id = vis[v];
            vis[v] = -1;
            v = new_id;
        }
        for (ui &v : C)
        {
            ui new_id = vis[v];
            vis[v] = -1;
            v = new_id;
        }
        for (ui &v : X)
        {
            ui new_id = vis[v];
            vis[v] = -1;
            v = new_id;
        }
    }
    ~Graph()
    {
    }
    Graph &operator=(const Graph &other)
    {
        n = other.n;
        m = other.m;
        vertices = other.vertices;
        return *this;
    }
    /**
     * @brief flip each edge, i.e., E' = {(u,v) | (v,u) \in E}
     */
    void flip()
    {
        for (auto &vertex : vertices)
        {
            swap(vertex.din, vertex.dout);
            swap(vertex.in_neighbors, vertex.out_neighbors);
        }
    }
    /**
     * @return whether (a, b) is in E
     */
    bool exist_edge(ui a, ui b)
    {
        return has(vertices[a].out_neighbors, b);
    }
    /**
     * @brief read edges from file where the file format can be ".txt" ".bin" ".out"
     *
     * ".txt" ".out" are text format
     * ".bin"  is binary format
     */
    void readFromFile(string file_path)
    {
        string suffix = get_file_name_suffix(file_path);
        if (suffix == "bin")
        {
            FILE *in = fopen(file_path.c_str(), "rb");
            if (in == nullptr)
            {
                printf("Failed to open %s \n", file_path.c_str());
                exit(1);
            }
            ui size_int;
            fread(&size_int, sizeof(ui), 1, in);
            if (size_int != sizeof(ui))
            {
                printf("sizeof int is different: graph_file(%u), machine(%u)\n", size_int, (int)sizeof(ui));
                exit(1);
            }
            fread(&n, sizeof(ui), 1, in);
            fread(&m, sizeof(ui), 1, in);
            cout << "File: " << get_file_name_without_suffix(file_path) << " n= " << n << " m= " << m << endl;
            ui *d = new ui[n]; // d[u] is the number of u's out-neighbors
            ui *pstart = new ui[n + 1];
            ui *edge_to = new ui[m];
            fread(d, sizeof(ui), n, in);
            fread(edge_to, sizeof(ui), m, in);
            pstart[0] = 0;
            for (ui i = 1; i <= n; i++)
                pstart[i] = pstart[i - 1] + d[i - 1];
            vertices.resize(n);
            for (ui u = 0; u < n; u++)
            {
                auto &vertex = vertices[u];
                vertex.out_neighbors.resize(d[u]);
                for (ui i = pstart[u], p = 0; i < pstart[u + 1]; i++, p++)
                {
                    ui j = edge_to[i];
                    vertex.out_neighbors[p] = j;
                    vertices[j].in_neighbors.push_back(u);
                    // printf("%d %d\n",u+1,j+1);
                }
            }
            for (ui u = 0; u < n; u++)
            {
                auto &vertex = vertices[u];
                vertex.id = vertex.origin_id = u;
                vertex.dout = vertex.out_neighbors.size();
                vertex.din = vertex.in_neighbors.size();
                vertex.in_neighbors.shrink_to_fit();
            }
            delete[] d;
            delete[] pstart;
            delete[] edge_to;
        }
        else // default graph file format: n m \n edges
        {
            ifstream in(file_path);
            if (!in.is_open())
            {
                printf("Failed to open %s \n", file_path.c_str());
                fflush(stdout);
                exit(1);
            }
            in >> n >> m;
            cout << "File: " << get_file_name_without_suffix(file_path) << " n= " << n << " m= " << m << endl;
            vector<pii> edges(m);
            ui idx = 0;
            for (ui i = 0; i < m; i++)
            {
                ui a, b;
                in >> a >> b;
                if (a == b)
                    continue;
                assert(a < n && b < n);
                edges[idx++] = {a, b};
            }
            edges.resize(idx);
            unique_pii(edges, n);
            m = edges.size();
            vertices.resize(n);
            for (auto &h : edges)
            {
                ui a = h.x, b = h.y;
                vertices[a].out_neighbors.push_back(b);
                vertices[b].in_neighbors.push_back(a);
            }
            for (ui u = 0; u < n; u++)
            {
                auto &vertex = vertices[u];
                vertex.id = vertex.origin_id = u;
                vertex.dout = vertex.out_neighbors.size();
                vertex.din = vertex.in_neighbors.size();
                vertex.out_neighbors.shrink_to_fit();
                vertex.in_neighbors.shrink_to_fit();
            }
            in.close();
        }
        printf("Graph init ok\n");
        fflush(stdout);
    }
    /**
     * @brief v_1, v_2, ..., v_n is the pseudo-degeneracy order, which will be dumped to vertex_id[0...n-1]
     * @return degeneracy of G
     */
    ui compute_degeneracy_order(int k, int l, vector<ui> &vertex_id)
    {
        assert(!vertex_id.size());
        assert(k <= l);
        vector<ui> din(n), dout(n);
        for (ui i = 0; i < n; i++)
            din[i] = vertices[i].din, dout[i] = vertices[i].dout;
        vector<ui> pd(n); // pseudo-degree
        for (ui i = 0; i < n; i++)
            pd[i] = min(din[i] + l, dout[i] + k);
        LinearHeap heap(n + l, n, pd);
        vertex_id.clear();
        vector<bool> popped(n);
        ui ret = 0;
        while (heap.sz)
        {
            ui u = heap.get_min_node();
            heap.delete_node(u);
            vertex_id.push_back(u);
            popped[u] = 1;
            ret = max(ret, pd[u]);
            for (ui v : vertices[u].out_neighbors)
            {
                if (popped[v])
                    continue;
                din[v]--;
                if (pd[v] > din[v] + l)
                {
                    pd[v] = din[v] + l;
                    heap.decrease(pd[v], v);
                }
            }
            for (ui v : vertices[u].in_neighbors)
            {
                if (popped[v])
                    continue;
                dout[v]--;
                if (pd[v] > dout[v] + k)
                {
                    pd[v] = dout[v] + k;
                    heap.decrease(pd[v], v);
                }
            }
        }
        ret = min(ret, n);
        return ret;
    }

    /**
     * @brief given a permutation of {0,...,n-1}, resort V
     * @param vertex_id v in new graph corresponds to vertex_id[v] in old graph
     */
    void re_sort(vector<ui> &vertex_id)
    {
        assert(vertex_id.size() == n);
        vector<ui> id_in_new(n); // v in old graph corresponds to id_in_new[v] in new graph
        for (ui i = 0; i < n; i++)
        {
            id_in_new[vertex_id[i]] = i;
        }
        for (auto &s : vertices)
        {
            s.id = id_in_new[s.id];
            for (auto &v : s.out_neighbors)
            {
                v = id_in_new[v];
            }
            sort(s.out_neighbors.begin(), s.out_neighbors.end());
            for (auto &v : s.in_neighbors)
            {
                v = id_in_new[v];
            }
            sort(s.in_neighbors.begin(), s.in_neighbors.end());
        }
        sort(vertices.begin(), vertices.end());
    }

    /**
     * @brief show how many CC in G
     */
    void connected_component()
    {
        vector<int> vis(n);
        int cc_cnt = 0;
        for (ui i = 0; i < n; i++)
        {
            if (!vis[i])
            {
                ++cc_cnt;
                queue<ui> q;
                q.push(i);
                vis[i] = cc_cnt;
                while (q.size())
                {
                    ui u = q.front();
                    q.pop();
                    for (auto v : vertices[u].out_neighbors)
                    {
                        if (!vis[v])
                        {
                            vis[v] = cc_cnt;
                            q.push(v);
                        }
                    }
                    for (auto v : vertices[u].in_neighbors)
                    {
                        if (!vis[v])
                        {
                            vis[v] = cc_cnt;
                            q.push(v);
                        }
                    }
                }
            }
        }
        printf("#ConnectedComponent= %d\n", cc_cnt);
    }
};

#endif