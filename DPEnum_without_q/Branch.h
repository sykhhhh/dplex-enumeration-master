#pragma once
#ifndef BRANCH_H
#define BRANCH_H

#include "Graph.h"
#include "MyBitset.h"
class BranchInstance
{
	using Set = MyBitset;

public:
	ui n;
	AdjacentMatrix *Aout; // Aout[u][v]=1 <=> u->v
	AdjacentMatrix *Ain;  // Ain[u][v]=1 <=> u<-v
	vector<ui> *vertex_map;
	ui k, l;
	Set S, C, X;
	BranchInstance(AdjacentMatrix *out, AdjacentMatrix *in,
				   vector<ui> *_vertex_map,
				   vector<ui> &_S, vector<ui> &_C, vector<ui> &_X,
				   ui _k, ui _l) : Aout(out), Ain(in), vertex_map(_vertex_map), k(_k), l(_l)
	{
		n = _S.size() + _C.size() + _X.size();
		assert(n == out->A.size());
		S = Set(n);
		for (ui v : _S)
			S.set(v);
		C = Set(n);
		for (ui v : _C)
			C.set(v);
		X = Set(n);
		for (ui v : _X)
			X.set(v);
		// reduce C and X
		for (ui v : S)
		{
			ui dout = Aout[0][v].intersect(S);
			if (dout + k == S.size())
			{
				C &= Aout[0][v];
				X &= Aout[0][v];
			}
			ui din = Ain[0][v].intersect(S);
			if (din + l == S.size())
			{
				C &= Ain[0][v];
				X &= Ain[0][v];
			}
		}
		for (ui u : C)
		{
			ui dout = Aout[0][u].intersect(S);
			if (dout + k < S.size() + 1)
			{
				C.reset(u);
				continue;
			}
			ui din = Ain[0][u].intersect(S);
			if (din + l < S.size() + 1)
			{
				C.reset(u);
				continue;
			}
		}
	}

	BranchInstance &operator=(const BranchInstance &other)
	{
		n = other.n;
		Aout = other.Aout;
		Ain = other.Ain;
		vertex_map = other.vertex_map;
		k = other.k;
		l = other.l;
		S = other.S;
		C = other.C;
		X = other.X;
		return *this;
	}

	/**
	 * @brief see paper: DPEnumPivot
	 * @return a pivot from C; if C is empty, return -1
	 */
	ui select_pivot()
	{
		ui min_pd = INF, sel = -1;
		bool removed = 0;
		auto V = S;
		V |= C;
		for (ui u : C)
		{
			ui dout = Aout[0][u].intersect(V);
			ui din = Ain[0][u].intersect(V);
			ui pd = min(dout + k, din + l);
			if (min_pd > pd)
			{
				min_pd = pd;
				sel = u;
			}
		}
		if (sel == -1)
			return -1;
		if (min_pd >= V.size() && !removed) // S+C is DPlex
		{
			bool ok = 1;
			for (ui v : S)
			{
				ui dout = Aout[0][v].intersect(V);
				ui din = Ain[0][v].intersect(V);
				ui pd = min(dout + k, din + l);
				if (pd < V.size())
				{
					ok = 0;
					break;
				}
			}
			if (ok)
			{
				S |= C;
				C.clear();
				for (ui v : S)
				{
					ui dout = Aout[0][v].intersect(S);
					if (dout + k == S.size())
					{
						X &= Aout[0][v];
					}
					ui din = Ain[0][v].intersect(S);
					if (din + l == S.size())
					{
						X &= Ain[0][v];
					}
				}
				return -1;
			}
		}
		return sel;
	}

	void remove_from_C(ui u)
	{
		C.reset(u);
		X.set(u);
	}

	/**
	 * @return false if u+S is not a DPlex
	 */
	bool move_from_C_to_S(ui u)
	{
		C.reset(u);
		ui dout = Aout[0][u].intersect(S);
		ui din = Ain[0][u].intersect(S);
		ui pd = min(dout + k, din + l);
		if (pd < S.size() + 1)
		{
			// assert(0);
			return false;
		}
		S.set(u);
		for (ui v : S)
		{
			if (!Aout[0][u][v])
			{
				ui din = Ain[0][v].intersect(S);
				if (din + l == S.size())
				{
					C &= Ain[0][v];
					X &= Ain[0][v];
				}
			}
			if (!Ain[0][u][v])
			{
				ui dout = Aout[0][v].intersect(S);
				if (dout + k == S.size())
				{
					C &= Aout[0][v];
					X &= Aout[0][v];
				}
			}
		}
		return true;
	}

	/**
	 * @brief check whether S is maximal
	 * @return if no vertex in X can be inserted to S, return true
	 */
	bool is_maximal()
	{
		assert(!C.size());
		for (ui u : X)
		{
			ui dout = Aout[0][u].intersect(S);
			ui din = Ain[0][u].intersect(S);
			ui pd = min(dout + k, din + l);
			if (pd >= S.size() + 1)
				return false;
		}
		return true;
	}

	/**
	 * @brief T(n)=O(n)
	 * @return ub=min_{v in S} { pd[v] }
	 */
	ui get_upper_bound()
	{
		ui ret = S.size() + C.size();
		auto V = S;
		V |= C;
		// vector<ui>pd_list;
		for (ui v : S)
		{
			ui dout = Aout[0][v].intersect(V);
			ui din = Ain[0][v].intersect(V);
			ui pd = min(dout + k, din + l);
			// pd_list.push_back(pd);
			ret = min(ret, pd);
		}
		return ret;
	}

	bool is_plex()
	{
		for (ui v : S)
		{
			ui dout = Aout[0][v].intersect(S);
			ui din = Ain[0][v].intersect(S);
			if (min(dout + k, din + l) < S.size())
				return false;
		}
		return true;
	}
};

/**
 * @brief instance of BK algorithm: partial set S
 * make sure after pushing v, S is still a DPlex
 */
class BK_Instance
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
	BK_Instance(Graph &_g, vector<ui> &v_map, vector<ui> &S,
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
	~BK_Instance() {}
	BK_Instance &operator=(const BK_Instance &other)
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
 * @brief branching instance for DPEnumPivot
 */

/**
 * @brief conduct recursive BK algorithm or DPEnumPivot
 */
class Branch
{
	Graph &g;
	Output &out;
	int k, l;
	int max_plex_size;
	vector<int> int_array_size_N_value_0;	  // default value: 0
	vector<int> int_array_size_N_value_0_another; // default value: 0
	vector<int> int_array_size_N_value_neg1;	  // default value: -1

public:
	ll dfs_cnt=0,cnt_small=0,cnt_large=0;
	double reduce_before_bk_time;

	double DPEnumPivot_time;
	double initial_bitset_time;
	map<int, int> cnt;
	Average subgraph_size;
	Branch(Graph &_g, Output &_o, int paramK, int paramL) : g(_g), out(_o), k(paramK), l(paramL),
															reduce_before_bk_time(0),
															DPEnumPivot_time(0), initial_bitset_time(0),
															subgraph_size(Average("avg-Graph-size"))
	{
		dfs_cnt = 0;
		int_array_size_N_value_0_another.resize(g.n);
		int_array_size_N_value_neg1.resize(g.n, -1);
		int_array_size_N_value_0.resize(g.n, 0);
	}
	~Branch() {}

	void print_results()
	{
		printf("%lld %lld\n",cnt_small,cnt_large);
		print_module_time("reduce-before-bk", reduce_before_bk_time);
		print_module_time("DPEnumPivot", DPEnumPivot_time);
		print_module_time("initial-bitset-time", initial_bitset_time);
		puts("");
		subgraph_size.print_averge_value();
		for (auto &h : cnt)
			cout << h.x << ' ' << h.y << endl;
	}

	void run()
	{
		small_DPlex();
		print_results();
	}

	/*-------------------------------- {enumerate small DPlex} -------------------------------------------------*/

	/**
	 * @brief enumerate maximal DPlexes size <= 2k-2
	 * step-1: IE-framework, include v_i, exclude v_1...v_{i-1}
	 * step-2: enumerate at most k-1 non-neighbors of v_i
	 * step-3 (maybe without this step): enumerate S of size k: v_i, v_i's non-neighbors, v_i's neighbors
	 * step-4: invoke BK
	 */
	void small_DPlex()
	{
		// sort vertices in g according to degeneracy order
		{
			vector<ui> vertices;
			max_plex_size = g.compute_degeneracy_order(k, l, vertices);
			g.re_sort(vertices);
		}

		Timer t;
		ui n = g.n;
		vector<bool> vis(n);	   // cache mask
		vector<int> TT(n);
		for (ui u = 0; u < n; u++) // step-1: IE
		{
			double percentage = 1.0 * u / n;
			print_progress_bar(percentage);
			// mark the out-neighbors of u
			for (ui v : g.vertices[u].out_neighbors)
				vis[v] = 1;
			// {v | (u,v) not in E and v>u}
			vector<ui> non_out_neighbors;
			for (ui v = u + 1; v < n; v++)
			{
				if (!vis[v])
				{
					non_out_neighbors.push_back(v);
				}
			}

			// step-ii: enumerate subsets of non-out-neighbors of u with size <= k-1
			vector<ui> S;
			S.reserve(k);
			S.push_back(u);
			dfs_enumerate_non_neighbor(non_out_neighbors, 0, S, k);
			// erase the out-neighbors of u
			for (ui v : g.vertices[u].out_neighbors)
				vis[v] = 0;
		}
		print_progress_bar(1.0, true);
		printf("#DPlex: %lld , use-time: %.4lf s\n", out.counter, t.get_time_seconds());
	}

	/**
	 * @brief step-2: enumerate subsets (non-out-neighbors) of size at most k-1
	 * @param max_size is equal to k
	 */
	void dfs_enumerate_non_neighbor(vector<ui> &father_set, ui start_of_father_set, vector<ui> &S, ui max_size)
	{
		if (S.size() == max_size)
		{
			// step-4: now |S| = k, we can invoke BK
			auto S_copy = S;
			// Candidate set C={v | (u,v) in E and v>u}
			// eXclusion set X={v | (u,v) in E and v<u}
			vector<ui> C, X;
			ui u = S[0];
			for (ui v : g.vertices[u].out_neighbors)
			{
				if (v < u)
					X.push_back(v);
				else
					C.push_back(v);
			}
			generate_set_for_BK_invocation(S_copy, C, X);
			if(S.size()==max_size) return;
		}
		ui position = S.size();
		S.push_back(0);
		for (ui i = start_of_father_set; i < father_set.size(); i++)
		{
			S[position] = father_set[i];
			dfs_enumerate_non_neighbor(father_set, i + 1, S, max_size);
		}
		S.pop_back();
		// step-3: there are less than k-1 non-neighbors, and we need to enumerate the out-neighbors of u
		ui u = S[0];
		auto &out_neighbors = g.vertices[u].out_neighbors;
		ui start_position = lower_bound(out_neighbors.begin(), out_neighbors.end(), u) - out_neighbors.begin();
		if (start_position < out_neighbors.size())
		{
			assert(out_neighbors[start_position] > u);
			dfs_enumerate_neighbor(out_neighbors, start_position, S, max_size);
		}
	}

	/**
	 * @brief step-3: enumerate subsets (out-neighbors) of size at most k-1
	 * @param max_size is equal to k
	 */
	void dfs_enumerate_neighbor(vector<ui> &father_set, ui start_of_father_set, vector<ui> &S, ui max_size)
	{
		if (S.size() + 1 == max_size)
		{
			ui position = S.size();
			S.push_back(0);
			for (ui i = start_of_father_set; i < father_set.size(); i++)
			{
				S[position] = father_set[i];
				// step-4: now |S| = k, we can invoke BK
				auto S_copy = S;
				// C+X must be a subset of in-neighbors of S
				vector<ui> C; // C = {v | (u,v) in E and v is greater than each vertex in S and v>u}
				vector<ui> X; // X={v | v is an in-neighbor of at least one vertex in S and v not in C}
				const int out_neighbor_of_u = 1;
				const int in_neighbor_of_S = 2;
				const int in_S = 3;
				const int in_C = 4;
				auto &vis = int_array_size_N_value_0_another;
				for (ui j = i + 1; j < father_set.size(); j++) // mark the out-neighbors of u
				{
					ui w = father_set[j];
					vis[w] = out_neighbor_of_u;
				}
				for (ui v : S) // mark the vertices in S (C+X shouldn't have any vertex in S)
					vis[v] = in_S;
				for (ui v : S) // find all in-neighbors of S
				{
					for (ui w : g.vertices[v].in_neighbors)
					{
						if (!vis[w])
						{
							vis[w] = in_neighbor_of_S;
							X.push_back(w);
						}
						else if (vis[w] == out_neighbor_of_u) // not only in-neighbor of S, but also a out-neighbor of u, then w in C
						{
							vis[w] = out_neighbor_of_u | in_neighbor_of_S;
							C.push_back(w);
						}
					}
				}
				// clear vis[]
				for (ui w : C)
					vis[w] = 0;
				for (ui w : X)
					vis[w] = 0;
				for (ui v : S)
					vis[v] = 0;
				generate_set_for_BK_invocation(S_copy, C, X);
			}
			S.pop_back();
			return;
		}
		ui position = S.size();
		S.push_back(0);
		for (ui i = start_of_father_set; i < father_set.size(); i++)
		{
			S[position] = father_set[i];
			dfs_enumerate_neighbor(father_set, i + 1, S, max_size);
		}
		S.pop_back();
		// if |S|<k, then S is not maximal
		return;
	}

	/**
	 * @brief remove u in C\cup X if S+u is not a DPlex; then invoke BK
	 */
	void generate_set_for_BK_invocation(vector<ui> &S, vector<ui> &C, vector<ui> &X)
	{
		assert(S.size() == k);
		if (C.size() || X.size()) // for u in C+X, we remove u if S+u is not a DPlex
		{
			Timer t;
			sort(C.begin(), C.end());
			sort(X.begin(), X.end());
			auto &in_neighbor_cnt = int_array_size_N_value_0_another;
			auto &out_neighbor_cnt = int_array_size_N_value_0;
			// for u in C+X, if u has more than k non-out-neighbors (or l non-in-neighbors) in S, then remove u
			{
				for (ui u : S)
				{
					ui dout = 0, din = 0;
					auto &out_neighbors = g.vertices[u].out_neighbors;
					auto &in_neighbors = g.vertices[u].in_neighbors;
					for (ui v : out_neighbors)
						in_neighbor_cnt[v]++;
					for (ui v : in_neighbors)
						out_neighbor_cnt[v]++;
				}
				ui C_size = 0;
				for (ui i = 0; i < C.size(); i++)
				{
					ui v = C[i];
					if (in_neighbor_cnt[v] + l < S.size() + 1 || out_neighbor_cnt[v] + k < S.size() + 1)
						continue;
					C[C_size++] = v;
				}
				C.resize(C_size);

				ui X_size = 0;
				for (ui i = 0; i < X.size(); i++)
				{
					ui v = X[i];
					if (in_neighbor_cnt[v] + l < S.size() + 1 || out_neighbor_cnt[v] + k < S.size() + 1)
						continue;
					X[X_size++] = v;
				}
				X.resize(X_size);
				// clear tool-array
				for (ui v : C)
					in_neighbor_cnt[v] = out_neighbor_cnt[v] = 0;
				for (ui v : X)
					in_neighbor_cnt[v] = out_neighbor_cnt[v] = 0;
			}
			// for u in C+X, if u is not the in-neighbors of each critical_out vertex, then remove u
			{
				int critical_out_cnt = 0, critical_in_cnt = 0; // if v in S and |S\Nout(v)|=k, then critical_out_cnt++
				for (ui u : S)
				{
					ui dout = out_neighbor_cnt[u], din = in_neighbor_cnt[u];
					auto &out_neighbors = g.vertices[u].out_neighbors;
					auto &in_neighbors = g.vertices[u].in_neighbors;
					assert(min(dout + k, din + l) >= S.size());
					if (dout + k == S.size())
					{
						critical_out_cnt++;
						for (ui v : out_neighbors)
							in_neighbor_cnt[v]++;
					}
					if (din + l == S.size())
					{
						critical_in_cnt++;
						for (ui v : in_neighbors)
							out_neighbor_cnt[v]++;
					}
				}
				ui C_size = 0;
				for (ui i = 0; i < C.size(); i++)
				{
					ui v = C[i];
					if (in_neighbor_cnt[v] < critical_out_cnt || out_neighbor_cnt[v] < critical_in_cnt)
						continue;
					assert(in_neighbor_cnt[v] == critical_out_cnt && out_neighbor_cnt[v] == critical_in_cnt);
					C[C_size++] = v;
				}
				C.resize(C_size);

				ui X_size = 0;
				for (ui i = 0; i < X.size(); i++)
				{
					ui v = X[i];
					if (in_neighbor_cnt[v] < critical_out_cnt || out_neighbor_cnt[v] < critical_in_cnt)
						continue;
					assert(in_neighbor_cnt[v] == critical_out_cnt && out_neighbor_cnt[v] == critical_in_cnt);
					X[X_size++] = v;
				}
				X.resize(X_size);
				// clear tool-array
				for (ui u : S)
				{
					auto &out_neighbors = g.vertices[u].out_neighbors;
					auto &in_neighbors = g.vertices[u].in_neighbors;
					for (ui v : out_neighbors)
						in_neighbor_cnt[v] = 0;
					for (ui v : in_neighbors)
						out_neighbor_cnt[v] = 0;
				}
			}
			reduce_before_bk_time += t.get_time();
		}
		// if (k == 2 && C.size()) // size > 2k-2 = 2, we ignore it
		// 	return;
		if (!X.size() && C.size() <= 1)
		{
			if (C.size() == 1)
			{
				S.push_back(C[0]);
			}
			for (ui &v : S)
			{
				v = g.vertices[v].origin_id;
			}
			out.dump(S);
			return;
		}
		else if (!C.size()) // S is not maximal
			return;
		Timer t_bitset;
		// g = G[S+C+X], and we use adjacent matrix to store g
		ui subg_n = S.size() + C.size() + X.size();
		subgraph_size.add(subg_n);
		// generate adjacent matrix
		AdjacentMatrix Aout(subg_n); // Aout[u][v]=1 <=> u->v
		AdjacentMatrix Ain(subg_n);  // Ain[u][v]=1 <=> u<-v
		// we first consider edges in G[S+C]
		ui idx = 0;
		auto &vis=int_array_size_N_value_neg1;
		for (ui v : S)
			vis[v] = idx++;
		for (ui v : C)
			vis[v] = idx++;
		for (ui u : S)
		{
			ui idx_u = vis[u];
			for (ui v : g.vertices[u].out_neighbors)
			{
				if (vis[v] == -1)
					continue;
				ui idx_v = vis[v];
				Aout.add_edge(idx_u, idx_v);
				Ain.add_edge(idx_v, idx_u);
			}
		}
		for (ui u : C)
		{
			ui idx_u = vis[u];
			for (ui v : g.vertices[u].out_neighbors)
			{
				if (vis[v] == -1)
					continue;
				ui idx_v = vis[v];
				Aout.add_edge(idx_u, idx_v);
				Ain.add_edge(idx_v, idx_u);
			}
		}
		// we then consider edges between S+C and X
		ui S_C_size = idx;
		assert(S_C_size == S.size() + C.size());
		for (ui v : X)
			vis[v] = idx++;
		for (ui u : X)
		{
			ui idx_u = vis[u];
			for (ui v : g.vertices[u].out_neighbors)
			{
				if (vis[v] >= S_C_size || vis[v] == -1)
					continue;
				ui idx_v = vis[v];
				Aout.add_edge(idx_u, idx_v);
				Ain.add_edge(idx_v, idx_u);
			}
			for (ui v : g.vertices[u].in_neighbors)
			{
				if (vis[v] >= S_C_size || vis[v] == -1)
					continue;
				ui idx_v = vis[v];
				Ain.add_edge(idx_u, idx_v);
				Aout.add_edge(idx_v, idx_u);
			}
		}
		vector<ui> vertex_map(subg_n);
		for (ui v : S)
			vertex_map[vis[v]] = g.vertices[v].origin_id;
		for (ui v : C)
			vertex_map[vis[v]] = g.vertices[v].origin_id;
		for (ui v : X)
			vertex_map[vis[v]] = g.vertices[v].origin_id;
		// clear vis[] and re-map S+C+X to id in g
		for (ui &v : S)
		{
			ui id = vis[v];
			vis[v] = -1;
			v = id;
		}
		for (ui &v : C)
		{
			ui id = vis[v];
			vis[v] = -1;
			v = id;
		}
		for (ui &v : X)
		{
			ui id = vis[v];
			vis[v] = -1;
			v = id;
		}
		BranchInstance instance(&Aout, &Ain, &vertex_map, S, C, X, k, l);
		initial_bitset_time+=t_bitset.get_time();
		Timer t;
		DPEnumPivot(instance); // begin recursive search
		DPEnumPivot_time += t.get_time();
	}

	/**
	 * @brief recursive search with better time complexity
	 */
	void DPEnumPivot(BranchInstance &instance)
	{
		dfs_cnt++;
		ui pivot = instance.select_pivot();
		if (pivot == -1) // C is empty
		{
			if (instance.is_maximal())
			{
				vector<ui> res;
				for (ui u : instance.S)
				{
					res.push_back(instance.vertex_map[0][u]);
				}
				out.dump(res);
			}
			return;
		}
		if (instance.C.size() == 1)
		{
			assert(instance.C[pivot]);
			if (instance.move_from_C_to_S(pivot))
				DPEnumPivot(instance);
			else
				DPEnumPivot(instance);
			return;
		}
		auto copy_instance = instance;
		// branch 1: exclude
		instance.remove_from_C(pivot);
		DPEnumPivot(instance);

		// branch 2: include
		if (copy_instance.move_from_C_to_S(pivot))
		{
			DPEnumPivot(copy_instance);
		}
	}
};

#endif