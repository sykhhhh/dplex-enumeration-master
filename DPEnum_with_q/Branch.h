#pragma once
#ifndef BRANCH_H
#define BRANCH_H

#include "Graph.h"
#include "MyBitset.h"

/**
 * @brief branching instance for DPEnumPivot
 */
class BranchInstance
{
	using Set = MyBitset;

public:
	ui n;
	AdjacentMatrix *Aout; // Aout[u][v]=1 <=> u->v
	AdjacentMatrix *Ain;  // Ain[u][v]=1 <=> u<-v
	vector<ui> *vertex_map;
	ui k, l,q;
	Set S, C, X;
	BranchInstance(AdjacentMatrix *out, AdjacentMatrix *in,
				   vector<ui> *_vertex_map,
				   vector<ui> &_S, vector<ui> &_C, vector<ui> &_X,
				   ui _k, ui _l,ui _q) : Aout(out), Ain(in), vertex_map(_vertex_map), k(_k), l(_l),q(_q)
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
			if (pd < q)
			{
				C.reset(u);
				V.reset(u);
				removed = 1;
				continue;
			}
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
		// return C.empty()?-1:*C.begin();
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
		for (ui v : S)
		{
			ui dout = Aout[0][v].intersect(V);
			ui din = Ain[0][v].intersect(V);
			ui pd = min(dout + k, din + l);
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
 * @brief conduct recursive BK algorithm or DPEnumPivot
 */
class Branch
{
	using Set = MyBitset;
	Graph &g;
	Output &out;
	int k, l,q;
	ui max_plex_size;
	vector<int> int_array_size_N_value_0;	  // default value: 0
	vector<int> int_array_size_N_value_0_another; // default value: 0
	vector<int> int_array_size_N_value_neg1;	  // default value: -1

public:
	ll dfs_cnt=0,cnt_small=0,cnt_large=0;
	double bk_time;
	double bk_copy_time;
	double bk_select_time;
	double bk_reduce_time;
	double reduce_before_bk_time;
	double induce_before_bk_time;

	double DPEnumPivot_time;
	double prepare_DPEnumPivot_time;
	double reduce_C_before_DPEnumPivot_time;
	double initial_bitset_time;

	double small_DPlex_time;
	double large_DPlex_time;
	map<int, int> cnt;
	Average subgraph_size;
	Branch(Graph &_g, Output &_o, ui paramK, ui paramL,ui paramQ) : g(_g), out(_o), k(paramK), l(paramL),q(paramQ),
															bk_time(0), bk_select_time(0),bk_copy_time(0), bk_reduce_time(0), reduce_before_bk_time(0),
															induce_before_bk_time(0), DPEnumPivot_time(0), prepare_DPEnumPivot_time(0),
															small_DPlex_time(0), large_DPlex_time(0), reduce_C_before_DPEnumPivot_time(0),initial_bitset_time(0),
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
		print_module_time("DPEnumPivot", DPEnumPivot_time);
		print_module_time("prepare-DPEnumPivot", prepare_DPEnumPivot_time);
		print_module_time("reduce-C-before-DPEnumPivot", reduce_C_before_DPEnumPivot_time);
		print_module_time("initial-bitset", initial_bitset_time);
		puts("");
		subgraph_size.print_averge_value();
		for (auto &h : cnt)
			cout << h.x << ' ' << h.y << endl;
	}

	void run()
	{
		{
			Timer t;
			large_DPlex();
			large_DPlex_time = t.get_time();
		}
		print_results();
	}

	/*-------------------------------- {enumerate large DPlex} -------------------------------------------------*/

	/**
	 * @brief enumerate maximal DPlexes size >= 2k-1
	 * step-1: second order reduction, remove vertices and edges that can not appear in any DPlex larger than 2k-2
	 * step-2: IE-framework, include v_i, exclude v_1...v_{i-1}
	 * step-3: enumerate at most k-1 non-out-neighbors
	 * step-4: recursive procedure
	 */
	void large_DPlex()
	{
		// step-1: second order reduction
		g.strong_reduce(k, l,q);
		printf("now: n= %u m= %u\n", g.n, g.m);
		// sort vertices in g according to degeneracy order
		{
			vector<ui> vertices;
			max_plex_size = g.compute_degeneracy_order(k, l, vertices);
			g.re_sort(vertices);
		}
		// resize the tool-arrays
		int_array_size_N_value_0_another.resize(g.n);
		int_array_size_N_value_0.resize(g.n);
		int_array_size_N_value_neg1.resize(g.n);

		// step-2: IE & 2-hops
		auto &vis = int_array_size_N_value_neg1;
		for (ui u = 0; u < g.n; u++)
		{
			double percentage = u * 1.0 / g.n;
			// if(u>=3100) printf("%d %d\n",u,out.counter);
			print_progress_bar(percentage);
			vector<ui> C;
			for(ui v:g.vertices[u].out_neighbors)
			{
				if(v<u) continue;
				vis[v]=1;
				for(ui x:g.vertices[v].out_neighbors)
				{
					if(x<u) continue;
					vis[x]=1;
				}
			}
			vis[u]=2;
			for(ui v:g.vertices[u].in_neighbors)
			{
				if(v<u) continue;
				if(vis[v]==1) C.push_back(v),vis[v]=2;
				for(ui x:g.vertices[v].in_neighbors)
				{
					if(x<u) continue;
					if(vis[x]==1) C.push_back(x),vis[x]=2;
				}
			}
			for(ui v:g.vertices[u].out_neighbors)
			{
				if(v<u) continue;
				vis[v]=-1;
				for(ui x:g.vertices[v].out_neighbors)
				{
					if(x<u) continue;
					vis[x]=-1;
				}
			}
			vis[u]=-1;
			// generate subgraph and reduce it
			if (1 + C.size() >= q)
			{
				Timer t;
				vector<ui> temp_C = C;
				temp_C.push_back(u);
				unordered_map<ui, ui> inverse_C; // if inv[w]=v in C, then w is the origin id of v
				for (ui v : C)
					inverse_C[g.vertices[v].origin_id] = v;
				Graph subg(temp_C, vis, g);
				ui prev_n = subg.n;
				do // reduce C
				{
					prev_n = subg.n;
					// find the new-index of u in subgraph
					ui origin_id = g.vertices[u].origin_id;
					ui idx = -1;
					for (auto &vertex : subg.vertices)
						if (vertex.origin_id == origin_id)
						{
							idx = vertex.id;
							break;
						}
					if (idx == -1) // this means u is removed from subg
					{
						subg.n = 0;
						break;
					}
					// reduce based on 2hop
					subg.two_hop_reduction(k, l,q, idx);
					subg.vertex_reduction(k, l,q);
					// break; // we find iteratively removing is very time-consuming, thus just 1 round
				} while (prev_n > subg.n && subg.n >= q);
				reduce_C_before_DPEnumPivot_time += t.get_time();
				if (subg.n >= q)
				{
					C.clear();
					ui u_origin_id = g.vertices[u].origin_id;
					for (auto &vertex : subg.vertices)
					{
						if (vertex.origin_id == u_origin_id)
							continue;
						assert(inverse_C.count(vertex.origin_id));
						ui v = inverse_C[vertex.origin_id];
						C.push_back(v);
					}
					// step-3: enumerate at most k-1 non-out-neighbor
					enumerate_non_out_neighbors(u, C);
				}
			}
		}
		print_progress_bar(1.0, true);
	}

	/**
	 * @brief step-3 of large DPlex: enumerate at most k-1 non-neighbors
	 */
	void push_to_d(ui u,ui k,vector<ui> &din_SC,vector<ui> &dout_SC)
	{
		for(ui v:g.vertices[u].out_neighbors)
		{
			din_SC[v]+=k;
		}
		for(ui v:g.vertices[u].in_neighbors)
		{
			dout_SC[v]+=k;
		}
	}
	void S_reduce(ui u,int flag,vector<ui> &S,vector<ui> &C,map<pair<ui,ui>,bool> &mp)
	{
		auto &vis_in=int_array_size_N_value_0_another;
		auto &vis_out=int_array_size_N_value_0;
		auto &vis=int_array_size_N_value_neg1;
		for(ui u:C) vis[u]=1;
		for (ui v : g.vertices[u].in_neighbors)
			vis_in[v] = 1;
		for (ui v : g.vertices[u].out_neighbors)
			vis_out[v] = 1;
		for(ui i=0;i<S.size();i++)
		{
			ui v=S[i];
			int cnt_in=0,cnt_out=0,cnt_in_out=0,cnt_out_in=0;
			if(vis_in[v]&&vis_out[v])
			{
				cnt_in=q-k-max(l-1-flag,0)-max(l-1-flag,0);
				cnt_out=q-k-max(k-1-flag,0)-max(k-1-flag,0);
				cnt_in_out=q-k-max(k-1-flag,0)-max(l-1-flag,0);
				cnt_out_in=q-k-max(k-1-flag,0)-max(l-1-flag,0);
			}
			else if(vis_in[v])
			{
				cnt_in=q-k-max(l-1-flag,0)-max(l-2-flag,0);
				cnt_out=q-k-max(k-1-flag,0)-max(k-2-flag,0);
				cnt_in_out=q-k-max(k-1-flag,0)-max(l-1-flag,0);
				cnt_out_in=q-k-max(k-2-flag,0)-max(l-2-flag,0);
			}
			else if(vis_out[v])
			{
				cnt_in=q-k-max(l-1-flag,0)-max(l-2-flag,0);
				cnt_out=q-k-max(k-1-flag,0)-max(k-2-flag,0);
				cnt_in_out=q-k-max(k-2-flag,0)-max(l-2-flag,0);
				cnt_out_in=q-k-max(k-1-flag,0)-max(l-1-flag,0);
			}
			else
			{
				cnt_in=q-k-max(l-2-flag,0)-max(l-2-flag,0);
				cnt_out=q-k-max(k-2-flag,0)-max(k-2-flag,0);
				cnt_in_out=q-k-max(k-2-flag,0)-max(l-2-flag,0);
				cnt_out_in=q-k-max(k-2-flag,0)-max(l-2-flag,0);
			}
			int same_neighbor_in=0,same_neighbor_out=0,same_neighbor_out_in=0,same_neighbor_in_out=0;
			for (ui w : g.vertices[v].in_neighbors)
			{
				if(vis_in[w]&&vis[w]==1) same_neighbor_in++;
				if(vis_out[w]&&vis[w]==1) same_neighbor_out_in++;
			}
			for (ui w : g.vertices[v].out_neighbors)
			{
				if(vis_out[w]&&vis[w]==1) same_neighbor_out++;
				if(vis_in[w]&&vis[w]==1) same_neighbor_in_out++;
			}
			if (same_neighbor_in<cnt_in||same_neighbor_out<cnt_out||same_neighbor_out_in<cnt_out_in||same_neighbor_in_out<cnt_in_out)
				mp[make_pair(u,v)]=1;
		}
		for (ui v : g.vertices[u].in_neighbors)
			vis_in[v] = 0;
		for (ui v : g.vertices[u].out_neighbors)
			vis_out[v] = 0;
		for(ui u:C) vis[u]=-1;
	}
	void enumerate_non_out_neighbors(ui u, vector<ui> &C)
	{
		vector<ui> neighbor, non_neighbor;
		auto &vis = int_array_size_N_value_neg1; // cache: mark the out-neighbors of u
		for (ui v : g.vertices[u].out_neighbors)
			vis[v] = 1;
		for (ui v : C)
			if (vis[v] != -1)
				neighbor.push_back(v);
			else
				non_neighbor.push_back(v);
		// clear the cache
		for (ui v : g.vertices[u].out_neighbors)
			vis[v] = -1;
		if (neighbor.size() + k < q)
			return;
		vector<ui> S{u};
		S.reserve(k);
		assert(S.size() == 1);
		vector<ui>din_SC(g.n),dout_SC(g.n);
		push_to_d(u,1,din_SC,dout_SC);
		for(ui v:neighbor) push_to_d(v,1,din_SC,dout_SC);
		map<pair<ui,ui>,bool>mp;
		dfs_enumerate_non_neighbor_large(non_neighbor, 0, S, k, neighbor,din_SC,dout_SC,mp);
	}

	/**
	 * @brief step-3: enumerate subsets (non-out-neighbors) of size at most k-1 for large DPlex
	 * @param max_size is equal to k
	 */
	void dfs_enumerate_non_neighbor_large(vector<ui> &father_set, ui start_of_father_set,
										  vector<ui> &S, ui max_size, vector<ui> &out_neighbor,
										  vector<ui> &din_SC,vector<ui> &dout_SC,map<pair<ui,ui>,bool> &mp)
	{
		if (S.size() == max_size)
		{
			// step-4: now |S| = k
			start_DPEnumPivot(S, out_neighbor,din_SC,dout_SC);
			return;
		}
		ui position = S.size();
		S.push_back(0);
		for (ui i = start_of_father_set; i < father_set.size(); i++)
		{
			bool ok=1;
			if(!ok) continue;
			S[position] = father_set[i];
			push_to_d(father_set[i],1,din_SC,dout_SC);
			dfs_enumerate_non_neighbor_large(father_set, i + 1, S, max_size, out_neighbor,din_SC,dout_SC,mp);
			push_to_d(father_set[i],-1,din_SC,dout_SC);
		}
		S.pop_back();
		// |S| < k
		start_DPEnumPivot(S, out_neighbor,din_SC,dout_SC);
	}

	/**
	 * @brief entrance of DPEnumPivot: generate 1) subgraph including X and 2) branching instance
	 * @param S S[0]=u, S[1...] are non-out-neighbors of u
	 * @param out_neighbor out-neighbors of u that serve as candidate set
	 */
	ui get_pd(vector<ui> &S,vector<ui> &C,vector<ui> &din_SC,vector<ui> &dout_SC)
	{
		ui ret = INF;
		for (ui v : S)
		{
			ui pd = min(dout_SC[v] + k, din_SC[v] + l);
			ret = min(ret, pd);
		}
		for (ui v : C)
		{
			ui pd = min(dout_SC[v] + k, din_SC[v] + l);
			ret = min(ret, pd);
		}
		return ret;
	}
	ui find_upperbound(vector<ui> &S,vector<ui> &C,vector<ui> &din_SC,vector<ui> &dout_SC)
	{
		ui ret = INF;
		auto &pd_list=int_array_size_N_value_0;
		ui sz_pd=0;
		auto right=pd_list.begin();
		for (ui v : S)
		{
			ui pd = min(dout_SC[v] + k, din_SC[v] + l);
			pd_list[sz_pd++];right++;
			ret = min(ret, pd);
		}
		for(ui v:C)
		{
			ui pd = min(dout_SC[v] + k, din_SC[v] + l);
			pd_list[sz_pd++];right++;
		}
		sort(pd_list.begin(),right);
		for(ui i=0;i<sz_pd;i++)
		{
			if(pd_list[i]>=sz_pd-i)
			{
				ret=min(ret,(ui)pd_list[i]);
				break;
			}
		}
		for(ui i=0;i<sz_pd;i++) pd_list[i]=0;
		return ret;
	}
	ui is_maximal(vector<ui> &S,vector<ui> &C,vector<ui> &X,vector<ui> &din_SC,vector<ui> &dout_SC)
	{
		auto &cnt_X=int_array_size_N_value_0;
		auto &vis=int_array_size_N_value_neg1;
		ui size_X=0;
		for(ui u:X) vis[u]=1;
		for (ui v : S)
		{
			if(dout_SC[v]+k==S.size()+C.size())
			{
				size_X++;
				for (ui w : g.vertices[v].out_neighbors)
				{
					if(vis[w]==1) cnt_X[w]++;
				}
			}
			if(din_SC[v]+l==S.size()+C.size())
			{
				size_X++;
				for (ui w : g.vertices[v].in_neighbors)
				{
					if(vis[w]==1) cnt_X[w]++;
				}
			}
		}
		for (ui v : C)
		{
			if(dout_SC[v]+k==S.size()+C.size())
			{
				size_X++;
				for (ui w : g.vertices[v].out_neighbors)
				{
					if(vis[w]==1) cnt_X[w]++;
				}
			}
			if(din_SC[v]+l==S.size()+C.size())
			{
				size_X++;
				for (ui w : g.vertices[v].in_neighbors)
				{
					if(vis[w]==1) cnt_X[w]++;
				}
			}
		}
		bool ok=1;
		for(ui i:X)
		{
			if(cnt_X[i]!=size_X) continue;
			if(min(dout_SC[i]+k,din_SC[i]+l)>=S.size()+C.size()+1)
			{
				ok=0;
				break;
			}
		}
		for(ui u:X) vis[u]=-1,cnt_X[u]=0;
		return ok;
	}
	void start_DPEnumPivot(vector<ui> &_S, vector<ui> &out_neighbor,vector<ui> &din_SC,vector<ui> &dout_SC)
	{
		Timer t_prepare;
		auto S = _S;
		auto C = out_neighbor;
		if (find_upperbound(S,C,din_SC,dout_SC) < q)
			return;
		vector<ui> X; // X={ v | v in V\(S+C)  and v+S+C may contain a DPlex larger than 2k-2  }
		auto &vis = int_array_size_N_value_neg1;
		for (ui v : S)
			vis[v] = 1;
		for (ui v : C)
			vis[v] = 1;
		/*
		 for v in X:
			if v is out-neighbor or in-neighbor of u
			then we have: v<u
			else we have: v is 2hop of u, i.e., v has a out-neighbor in C
		*/
		ui u = S[0];
		for (ui v : g.vertices[u].out_neighbors)
		{
			if (v >= u) // note that for those v: v>u and u->v and v is not in C, then v can not appear in DPlex larger than 2k-2
				break;
			X.push_back(v);
			vis[v] = 1;
		}
		if (S.size() < k) // if |S|<k, X contains non-out-neighbors of u
		{
			for (ui v : g.vertices[u].in_neighbors)
			{
				if (v >= u)
					break;
				if (vis[v] != -1)
					continue;
				X.push_back(v);
				vis[v] = 1;
			}
			for (ui v : C)
			{
				for (ui w : g.vertices[v].in_neighbors)
				{
					if (vis[w] != -1)
						continue;
					X.push_back(w);
					vis[w] = 1;
				}
			}
		}
		// clear vis[]
		for (ui v : S)
			vis[v] = -1;
		for (ui v : C)
			vis[v] = -1;
		for (ui v : X)
			vis[v] = -1;
		// further reduce X based on degree
		{
			// record S+C
			for (ui v : S)
				vis[v] = 1;
			for (ui v : C)
				vis[v] = 1;
			ui X_sz = 0;
			for (ui i = 0; i < X.size(); i++)
			{
				ui v = X[i];
				ui dout = 0;
				for (ui w : g.vertices[v].out_neighbors)
					if (vis[w] == 1)
						dout++;
				if (dout + k < q)
					continue;
				ui din = 0;
				for (ui w : g.vertices[v].in_neighbors)
					if (vis[w] == 1)
						din++;
				if (din + l < q)
					continue;
				X[X_sz++] = v;
			}
			X.resize(X_sz);
			// clear vis[]
			for (ui v : S)
				vis[v] = -1;
			for (ui v : C)
				vis[v] = -1;
		}
		prepare_DPEnumPivot_time += t_prepare.get_time();
		ui pd=get_pd(S,C,din_SC,dout_SC);
		if(pd>=S.size()+C.size())
		{
			if(S.size()+C.size()>=q&&is_maximal(S,C,X,din_SC,dout_SC)) out.dump(S);
			return;
		}
		Timer t_bitset;
		// g = G[S+C+X], and we use adjacent matrix to store g
		ui subg_n = S.size() + C.size() + X.size();
		subgraph_size.add(subg_n);
		// generate adjacent matrix
		AdjacentMatrix Aout(subg_n); // Aout[u][v]=1 <=> u->v
		AdjacentMatrix Ain(subg_n);  // Ain[u][v]=1 <=> u<-v
		// we first consider edges in G[S+C]
		ui idx = 0;
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
		BranchInstance instance(&Aout, &Ain, &vertex_map, S, C, X, k, l,q);
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
		if (instance.get_upper_bound() < q)
			return;
		ui pivot = instance.select_pivot();
		if (pivot == -1) // C is empty
		{
			if (instance.S.size() >= q && instance.is_maximal())
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