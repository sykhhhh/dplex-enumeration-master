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
	double T;
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
		}
		// clear the tool-array
		for (ui i = 0; i < n; i++)
		{
			vis[S[i]] = -1;
		}
	}
	/**
	 * @brief induce the subgraph G[S + X]; and ignore the edges between X; note that S C X will be re-mapped
	 * @param vis tool-array: make sure vis[0...n-1]=-1 before and after this function
	 *
	 * @return each vertex in S,C,X will be re-mapped to the new index of subgraph
	 */
	Graph(vector<ui> &S, vector<ui> &C, vector<ui> &X, vector<int> &vis, const Graph &origin_g)
	{
		Timer t;
		vector<ui> V = C;
		V.insert(V.end(), S.begin(), S.end());
		n = V.size() + X.size();
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
		T=t.get_time();
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
				int u = vis[u_origin];
				if (u == -1)
					continue;
				vertex.out_neighbors.push_back(u);
				vertices[u].in_neighbors.push_back(v);
			}
		}
		// out-neighbors of X in S+C
		if (X.size() <= V.size())
		{
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
					int u = vis[u_origin];
					if (u == -1)
						continue;
					if (u >= V.size()) // u in X, then we ignore (v,u)
						continue;
					vertex.out_neighbors.push_back(u);
					vertices[u].in_neighbors.push_back(v);
				}
			}
		}
		else // |S+C| < X, we enumerate in-neighbors of S+C, which is equal to enumerate out-neighbors of X
		{
			for (ui v_origin : X)
			{
				ui v = vis[v_origin];
				assert(v >= V.size());
				auto &vertex = vertices[v];
				auto &origin_vertex = origin_g.vertices[v_origin];
				vertex.id = v;
				vertex.origin_id = origin_vertex.origin_id;
			}
			for (ui v_origin : V)
			{
				ui v = vis[v_origin];
				assert(v < V.size());
				auto &vertex = vertices[v];
				auto &origin_vertex = origin_g.vertices[v_origin];
				int V_size = V.size();
				for (ui u_origin : origin_vertex.in_neighbors)
				{
					int u = vis[u_origin];
					if (u < V_size) // u in V, then we ignore (u,v)
						continue;
					// v in S+C and u in X and u->v
					vertex.in_neighbors.push_back(u);
					vertices[u].out_neighbors.push_back(v);
				}
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
				}
			}
			for (ui u = 0; u < n; u++)
			{
				auto &vertex = vertices[u];
				vertex.id = vertex.origin_id = u;
				vertex.dout = vertex.out_neighbors.size();
				vertex.din = vertex.in_neighbors.size();
				// if(vertex.dout+vertex.din==0) cerr<<"!!!!"<<endl;
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
	 * @brief v_1, v_2, ..., v_n is the out-degeneracy order, which will be dumped to vertex_id[0...n-1]
	 * @return out-degeneracy of G
	 */
	ui compute_degeneracy_order(int k, int l, vector<ui> &vertex_id)
	{
		assert(!vertex_id.size());
		assert(k <= l);
		vector<ui> dout(n); // we only care about out-degree
		for (ui i = 0; i < n; i++)
			dout[i] = vertices[i].dout;
		LinearHeap heap(n, n, dout);
		vertex_id.clear();
		vector<bool> popped(n);
		ui ret = 0;
		while (heap.sz)
		{
			ui u = heap.get_min_node();
			heap.delete_node(u);
			vertex_id.push_back(u);
			popped[u] = 1;
			ret = max(ret, dout[u]);
			for (ui v : vertices[u].in_neighbors)
			{
				if (popped[v])
					continue;
				dout[v]--;
				heap.decrease(dout[v], v);
			}
		}
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
	 * @brief remove the incident edges of u
	 */
	void remove_vertex(ui u)
	{
		auto &out_neighbors = vertices[u].out_neighbors;
		auto &in_neighbors = vertices[u].in_neighbors;
		for (ui v : out_neighbors)
		{
			auto &vertex = vertices[v];
			int din = 0;
			for (ui w : vertex.in_neighbors)
				if (w != u)
					vertex.in_neighbors[din++] = w;
			assert(din + 1 == vertex.din);
			vertex.din = din;
			vertex.in_neighbors.resize(din);
		}
		for (ui v : in_neighbors)
		{
			auto &vertex = vertices[v];
			int dout = 0;
			for (ui w : vertex.out_neighbors)
				if (w != u)
					vertex.out_neighbors[dout++] = w;
			assert(dout + 1 == vertex.dout);
			vertex.dout = dout;
			vertex.out_neighbors.resize(dout);
		}
		m -= out_neighbors.size() + in_neighbors.size();
		out_neighbors.clear();
		in_neighbors.clear();
		vertices[u].din = 0;
		vertices[u].dout = 0;
	}

	/**
	 * @brief remove vertices and edges that can not appear in DPlex larger than 2k-2
	 */
	void strong_reduce(int k, int l,int q)
	{
		ui m_before = m;
		do
		{
			vertex_reduction(k, l,q);
			m_before = m;
			edge_reduction(k, l,q);
		} while (m_before > m);
	}

	/**
	 * @brief remove edges that can not appear in DPlex larger than 2k-2
	 */
	void edge_reduction(int k, int l,int q)
	{
		// only if k==l can we reduce edges
		// if (k != l)
		// 	return;
		// if u->v but not v->u, then we can remove (u,v) if |N^{in}(u) \cap N^{out}(v)| = 0
		vector<bool> vis_in(n),vis_out(n); // cache: vis[v]=1 <=> v->u
		ui previous_m = m;
		for (auto &vertex : vertices)
		{
			ui u = vertex.id;
			// record the in-neighbors of u
			for (ui v : vertex.in_neighbors)
				vis_in[v] = 1;
			for (ui v : vertex.out_neighbors)
				vis_out[v] = 1;
			auto &out_neighbors = vertex.out_neighbors;
			ui new_size = 0; // size of reduced N^{out}(u)
			for (ui i = 0; i < vertex.out_neighbors.size(); i++)
			{
				ui v = vertex.out_neighbors[i]; // u->v
				int lb_in,lb_out,lb_in_out,lb_out_in;
				if (vis_in[v])					 // v->u
				{
					lb_in=q-2*l;
					lb_out=q-2*k;
					lb_in_out=q-l-k;
					lb_out_in=q-l-k;
				}
				else
				{
					lb_in=q-2*l+1;
					lb_out=q-2*k+1;
					lb_in_out=q-l-k+2;
					lb_out_in=q-l-k;
				}
				// check if |N^{in}(u) \cap N^{out}(v)| = 0
				int same_neighbor_in=0,same_neighbor_out=0,same_neighbor_out_in=0,same_neighbor_in_out=0;
				bool removed = 1;
				for (ui w : vertices[v].in_neighbors)
				{
					if (vis_in[w]) same_neighbor_in++;
					if(vis_out[w]) same_neighbor_out_in++;
				}
				for (ui w : vertices[v].out_neighbors)
				{
					if (vis_out[w]) same_neighbor_out++;
					if(vis_in[w]) same_neighbor_in_out++;
				}
				if (same_neighbor_in>=lb_in&&same_neighbor_out>=lb_out&&same_neighbor_out_in>=lb_out_in&&same_neighbor_in_out>=lb_in_out)
					out_neighbors[new_size++] = v;
			}
			m -= out_neighbors.size() - new_size;
			out_neighbors.resize(new_size);
			vertex.dout = new_size;
			// clear the cache
			for (ui v : vertex.in_neighbors)
				vis_in[v] = 0;
			for (ui v : vertex.out_neighbors)
				vis_out[v] = 0;
		}
		// rebuild graph
		if (previous_m > m)
		{
			for (auto &vertex : vertices)
				vertex.in_neighbors.resize(0);
			for (auto &vertex : vertices)
			{
				ui u = vertex.id;
				for (ui v : vertex.out_neighbors)
				{
					vertices[v].in_neighbors.push_back(u);
				}
			}
			for (auto &vertex : vertices)
				vertex.din = vertex.in_neighbors.size();
		}
	}

	/**
	 * @brief remove vertices that can not appear in DPlex larger than 2k-2
	 * i.e., we remove u if pd[u] <= 2k-2
	 */
	void vertex_reduction(int k, int l,int q)
	{
		vector<int> in_queue(n);
		queue<ui> Q;
		for (auto &vertex : vertices)
		{
			ui u = vertex.id;
			ui pseudo_degree = min(vertex.din + l, vertex.dout + k);
			if (pseudo_degree < q)
			{
				in_queue[u] = 1;
				Q.push(u);
			}
		}
		if (!Q.size())
			return;
		while (Q.size())
		{
			ui u = Q.front();
			Q.pop();
			for (ui v : vertices[u].out_neighbors)
			{
				if (in_queue[v])
					continue;
				vertices[v].din--;
				if (vertices[v].din + l < q)
				{
					Q.push(v);
					in_queue[v] = 1;
				}
			}
			for (ui v : vertices[u].in_neighbors)
			{
				if (in_queue[v])
					continue;
				vertices[v].dout--;
				if (vertices[v].dout + k < q)
				{
					Q.push(v);
					in_queue[v] = 1;
				}
			}
		}
		// rebuild graph
		remove_given_vertices(in_queue);
	}

	/**
	 * @brief remove vertices based on diameter or 2-hop
	 * @param u current graph aims to find DPlexes containing u
	 */
	void two_hop_reduction(int k, int l,int q, int u)
	{
		vector<bool> vis_in(n,0),vis_out(n,0);
		vector<int> vertex_removed(n, 0);; // cache: vis[v]=1 <=> v->u
		auto &u_vertex=vertices[u];
		// record the in-neighbors of u
		for (ui v : u_vertex.in_neighbors)
			vis_in[v] = 1;
		for (ui v : u_vertex.out_neighbors)
			vis_out[v] = 1;
		ui new_size = 0; // size of reduced N^{out}(u)
		for (auto &vertex : vertices)
		{
			ui v = vertex.id;
			if (vertex.id == u)
				continue;
			int lb_in,lb_out,lb_in_out,lb_out_in;
			if (vis_in[v]&&vis_out[v])					 // v->u
			{
				lb_in=q-2*l;
				lb_out=q-2*k;
				lb_in_out=q-l-k;
				lb_out_in=q-l-k;
			}
			else if(vis_out[v])
			{
				lb_in=q-2*l+1;
				lb_out=q-2*k+1;
				lb_in_out=q-l-k+2;
				lb_out_in=q-l-k;
			}
			else if(vis_in[v])
			{
				lb_in=q-2*l+1;
				lb_out=q-2*k+1;
				lb_in_out=q-l-k;
				lb_out_in=q-l-k+2;
			}
			else
			{
				lb_in=q-2*l+2;
				lb_out=q-2*k+2;
				lb_in_out=q-l-k+2;
				lb_out_in=q-l-k+2;
			}
			// check if |N^{in}(u) \cap N^{out}(v)| = 0
			int same_neighbor_in=0,same_neighbor_out=0,same_neighbor_out_in=0,same_neighbor_in_out=0;
			bool removed = 1;
			for (ui w : vertices[v].in_neighbors)
			{
				if (vis_in[w]) same_neighbor_in++;
				if(vis_out[w]) same_neighbor_out_in++;
			}
			for (ui w : vertices[v].out_neighbors)
			{
				if (vis_out[w]) same_neighbor_out++;
				if(vis_in[w]) same_neighbor_in_out++;
			}
			if (same_neighbor_in<lb_in||same_neighbor_out<lb_out||same_neighbor_out_in<lb_out_in||same_neighbor_in_out<lb_in_out)
				vertex_removed[v] = 1;
		}
		for (ui v : u_vertex.in_neighbors)
			vis_in[v] = 0;
		for (ui v : u_vertex.out_neighbors)
			vis_out[v] = 0;
		remove_given_vertices(vertex_removed);
	}

	/**
	 * @brief remove vertices and rebuild graph
	 */
	void remove_given_vertices(vector<int> &vertex_removed)
	{
		ui idx = 0;
		auto &vertex_map = vertex_removed;
		for (ui i = 0; i < n; i++)
		{
			if (vertex_removed[i])
				vertex_map[i] = -1;
			else
				vertex_map[i] = idx++;
		}
		vector<Vertex> new_vertices(idx);
		for (ui u_origin = 0; u_origin < n; u_origin++)
		{
			if (vertex_map[u_origin] == -1)
				continue;
			ui u = vertex_map[u_origin];
			auto &vertex = new_vertices[u];
			vertex.id = u;
			auto &origin_vertex = vertices[u_origin];
			vertex.origin_id = origin_vertex.origin_id;
			for (ui v_origin : origin_vertex.out_neighbors)
			{
				int v = vertex_map[v_origin];
				if (v == -1)
					continue;
				vertex.out_neighbors.push_back(v);
				new_vertices[v].in_neighbors.push_back(u);
			}
		}
		swap(new_vertices, vertices);
		n = idx;
		m = 0;
		for (auto &vertex : vertices)
		{
			vertex.in_neighbors.shrink_to_fit();
			vertex.out_neighbors.shrink_to_fit();
			vertex.din = vertex.in_neighbors.size();
			vertex.dout = vertex.out_neighbors.size();
			m += vertex.din;
		}
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