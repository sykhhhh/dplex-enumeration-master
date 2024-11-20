#ifndef PREPRO_H
#define PREPRO_H

#include "Utility.h"
#include "Graph.h"

/**
 * @brief preprocessing: some vertices have no out-neighbors (or in-neighbors),
 * so we can first output all DPlexes related to those vertices, then we can remove them safely
 * and needn't consider whether those vertices can be inserted to the later DPlexes
 */
class Prepro
{
	Graph &g;
	Output &out;
	int k, l;

public:
	Prepro(Graph &_g, Output &_o, int paramK, int paramL) : g(_g), out(_o), k(paramK), l(paramL)
	{
	}
	~Prepro() {}

	/**
	 * @brief entrance of preprocessing stage
	 */
	void begin_preprocessing()
	{
#ifndef NO_PREPRO
	ui prev_n;
		do
		{
			prev_n=g.n;
			no_out_or_no_in();
			if (k == l && k == 2)
			{
				special_judge_when_k_l_2();
			}
			else if (k == l && k == 3)
			{
				special_judge_when_k_l_3();
			}
			simple_reduce();
			// break;
		} while(prev_n > g.n);
#endif
	}
	/**
	 * @brief preprocessing when (k,l) = (2,2): if the 2hop-in-neighbors of u and out-neighbors of u are not intersected, then UB(u)=2
	 */
	int dfs_out(int x,int t,vector<bool> &vis,vector<int> &in_queue)
	{
		if(t<0) return 0;
		vis[x]=1;
		for(auto y:g.vertices[x].out_neighbors)
		{
			if(in_queue[y]) continue;
			if(vis[y]||dfs_out(y,t-1,vis,in_queue)) return vis[x]=0,1;
		}
		vis[x]=0;
		return 0;
	}
	int dfs_in(int x,int t,vector<bool> &vis,vector<int> &in_queue)
	{
		if(t<0) return 0;
		vis[x]=1;
		for(auto y:g.vertices[x].in_neighbors)
		{
			if(in_queue[y]) continue;
			if(vis[y]||dfs_in(y,t-1,vis,in_queue)) return vis[x]=0,1;
		}
		vis[x]=0;
		return 0;
	}
	void simple_reduce()
	{
		ui n=g.n;
		queue<int>q;
		vector<bool>vis(n);
		vector<int>in_queue(n);
		vector<bool> no_directed(n, true); //|Nout(u) \cap Nin(u)| = 0
		for (ui u = 0; u < n; u++)		 // find the vertices without bi-directed edges
		{
			if (!no_directed[u])
				continue;
			for (ui v : g.vertices[u].out_neighbors)
				vis[v] = 1;
			no_directed[u] = 1;
			for (ui v : g.vertices[u].in_neighbors)
				if (vis[v]) // u->v and v->u
				{
					no_directed[u] = 0;
					no_directed[v] = 0;
					break;
				}
			for (ui v : g.vertices[u].out_neighbors)
				vis[v] = 0;
		}
		int is_reduced=1;
		while(is_reduced)
		{
			is_reduced=0;
			for(int i=0;i<n;i++)
			{
				if(in_queue[i]) continue;
				int flag=1;
				if(k<=l) flag&=dfs_out(i,k,vis,in_queue);
				// if(k>=l) flag&=dfs_in(i,l,vis,in_queue);
				if(!flag)
				{
					q.push(i);
					in_queue[i]=1;
					is_reduced=1;
				}
			}
		}
		vector<ui> candidate(n), index(n);
		for (ui i = 0; i < n; i++)
			candidate[i] = index[i] = i;
		while (q.size() && n >= k)
		{
			ui u = q.front();
			q.pop();
			// remove u from candidate
			{
				ui pos = index[u];
				assert(pos < n);
				assert(candidate[pos] == u);
				index[candidate[n - 1]] = pos;
				index[u] = n;
				swap(candidate[pos], candidate[n - 1]);
				candidate.pop_back();
				n--;
			}
			// now generate DPlex containing u using k-1 vertices in candidate[0...n-1]
			// enumerate all the subsets of size k-1 using depth-first-search (dfs)
			vector<ui> plex(k);
			plex[0] = u;
			dfs(candidate, 0, plex, 1);
		}
		g.remove_given_vertices(in_queue);
	}
	void special_judge_when_k_l_2()
	{
		ui n = g.n;
		queue<ui> q;
		vector<bool> vis(n); // cache & mark
		vector<int> in_queue(n);
		for (ui u = 0; u < n; u++)
		{
			auto &out_neighbors = g.vertices[u].out_neighbors;
			for (ui v : out_neighbors) // mark the out-neighbors
				vis[v] = 1;
			bool ok = 1; // ok <=> | Nout(u) \cap (Nin2hop(u)) | = 0
			for (ui v : g.vertices[u].in_neighbors)
				if (vis[v])
				{
					ok = 0;
					break;
				}
				else
				{
					for (ui w : g.vertices[v].in_neighbors)
						if (vis[w])
						{
							ok = 0;
							break;
						}
					if (!ok)
						break;
				}
			for (ui v : out_neighbors) // clear the cache
				vis[v] = 0;
			if (ok)
			{
				q.push(u);
				in_queue[u] = 1;
			}
		}
		if (!q.size())
			return;
		// if i in V: candidate[index[i]]=i; if i not in V, then index[i] > n
		// i.e., candidate is the vertex set of the not removed vertices
		vector<ui> candidate(n), index(n);
		for (ui i = 0; i < n; i++)
			candidate[i] = index[i] = i;
		while (q.size() && n >= k)
		{
			ui a = q.front();
			q.pop();
			// remove a from candidate
			{
				ui pos = index[a];
				assert(pos < n);
				assert(candidate[pos] == a);
				index[candidate[n - 1]] = pos;
				index[a] = n;
				swap(candidate[pos], candidate[n - 1]);
				candidate.pop_back();
				n--;
			}
			auto &vertex = g.vertices[a];
			ui u_origin = vertex.origin_id;
			// enumerate DPlexes containing u: their sizes are 2
			vector<ui> S(2);
			S[0] = u_origin;
			for (ui v : candidate)
			{
				S[1] = g.vertices[v].origin_id;
				out.dump(S, 2);
			}
			for (ui u : vertex.out_neighbors)
			{
				if (in_queue[u])
					continue;
				auto &out_neighbors = g.vertices[u].out_neighbors;
				for (ui v : out_neighbors)
					if (!in_queue[v])
						vis[v] = 1;
				bool ok = 1; // ok <=> | Nout(u) \cap (Nin2hop(u)) | = 0
				for (ui v : g.vertices[u].in_neighbors)
					if (vis[v])
					{
						ok = 0;
						break;
					}
					else if (!in_queue[v])
					{
						for (ui w : g.vertices[v].in_neighbors)
							if (vis[w])
							{
								ok = 0;
								break;
							}
						if (!ok)
							break;
					}
				for (ui v : out_neighbors)
					vis[v] = 0;
				if (ok)
				{
					q.push(u);
					in_queue[u] = 1;
				}
			}
			for (ui u : vertex.in_neighbors)
			{
				if (in_queue[u])
					continue;
				auto &out_neighbors = g.vertices[u].out_neighbors;
				for (ui v : out_neighbors)
					if (!in_queue[v])
						vis[v] = 1;
				bool ok = 1; // ok <=> | Nout(u) \cap (Nin2hop(u)) | = 0
				for (ui v : g.vertices[u].in_neighbors)
					if (vis[v])
					{
						ok = 0;
						break;
					}
					else if (!in_queue[v])
					{
						for (ui w : g.vertices[v].in_neighbors)
							if (vis[w])
							{
								ok = 0;
								break;
							}
						if (!ok)
							break;
					}
				for (ui v : out_neighbors)
					vis[v] = 0;
				if (ok)
				{
					q.push(u);
					in_queue[u] = 1;
				}
			}
		}
		g.remove_given_vertices(in_queue);
	}
	/**
	 * @brief preprocessing when (k,l) = (3,3): if the 2hop-in-neighbors of u and 2hop-out-neighbors of u are not intersected, then UB(u)=3
	 */
	void special_judge_when_k_l_3()
	{
		ui n = g.n;
		queue<ui> q;
		vector<bool> vis(n); // cache & mark
		vector<int> in_queue(n);
		for (ui u = 0; u < n; u++)
		{
			auto &out_neighbors = g.vertices[u].out_neighbors;
			for (ui v : out_neighbors) // mark the 2hop-out-neighbors
			{
				vis[v] = 1;
				for (ui w : g.vertices[v].out_neighbors)
					vis[w] = 1;
			}
			bool ok = 1; // ok <=> | Nout2hop(u) \cap (Nin2hop(u)) | = 0
			for (ui v : g.vertices[u].in_neighbors)
				if (vis[v])
				{
					ok = 0;
					break;
				}
				else
				{
					for (ui w : g.vertices[v].in_neighbors)
						if (vis[w])
						{
							ok = 0;
							break;
						}
					if (!ok)
						break;
				}
			for (ui v : out_neighbors) // clear the cache
			{
				vis[v] = 0;
				for (ui w : g.vertices[v].out_neighbors)
					vis[w] = 0;
			}
			if (ok)
			{
				q.push(u);
				in_queue[u] = 1;
			}
		}
		if (!q.size())
			return;
		// if i in V: candidate[index[i]]=i; if i not in V, then index[i] > n
		// i.e., candidate is the vertex set of the not removed vertices
		vector<ui> candidate(n), index(n);
		for (ui i = 0; i < n; i++)
			candidate[i] = index[i] = i;
		while (q.size() && n >= k)
		{
			ui a = q.front();
			q.pop();
			// remove a from candidate
			{
				ui pos = index[a];
				assert(pos < n);
				assert(candidate[pos] == a);
				index[candidate[n - 1]] = pos;
				index[a] = n;
				swap(candidate[pos], candidate[n - 1]);
				candidate.pop_back();
				n--;
			}
			auto &vertex = g.vertices[a];
			ui u_origin = vertex.origin_id;
			// enumerate DPlexes containing u: their sizes are 3
			vector<ui> S(3);
			S[0] = u_origin;
			for (ui v : candidate)
			{
				S[1] = g.vertices[v].origin_id;
				for (ui w : candidate)
				{
					if (w == v)
						break;
					S[3] = g.vertices[w].origin_id;
					out.dump(S, 3);
				}
			}
			for (ui u : vertex.out_neighbors)
			{
				if (in_queue[u])
					continue;
				auto &out_neighbors = g.vertices[u].out_neighbors;
				for (ui v : out_neighbors) // mark the 2hop-out-neighbors
				{
					if (in_queue[v])
						continue;
					vis[v] = 1;
					for (ui w : g.vertices[v].out_neighbors)
						if (!in_queue[w])
							vis[w] = 1;
				}
				bool ok = 1; // ok <=> | Nout2hop(u) \cap (Nin2hop(u)) | = 0
				for (ui v : g.vertices[u].in_neighbors)
					if (vis[v])
					{
						ok = 0;
						break;
					}
					else if (!in_queue[v])
					{
						for (ui w : g.vertices[v].in_neighbors)
							if (vis[w])
							{
								ok = 0;
								break;
							}
						if (!ok)
							break;
					}
				for (ui v : out_neighbors) // clear the cache
				{
					vis[v] = 0;
					for (ui w : g.vertices[v].out_neighbors)
						vis[w] = 0;
				}
				if (ok)
				{
					q.push(u);
					in_queue[u] = 1;
				}
			}
			for (ui u : vertex.in_neighbors)
			{
				if (in_queue[u])
					continue;
				auto &out_neighbors = g.vertices[u].out_neighbors;
				for (ui v : out_neighbors) // mark the 2hop-out-neighbors
				{
					if (in_queue[v])
						continue;
					vis[v] = 1;
					for (ui w : g.vertices[v].out_neighbors)
						if (!in_queue[w])
							vis[w] = 1;
				}
				bool ok = 1; // ok <=> | Nout2hop(u) \cap (Nin2hop(u)) | = 0
				for (ui v : g.vertices[u].in_neighbors)
					if (vis[v])
					{
						ok = 0;
						break;
					}
					else if (!in_queue[v])
					{
						for (ui w : g.vertices[v].in_neighbors)
							if (vis[w])
							{
								ok = 0;
								break;
							}
						if (!ok)
							break;
					}
				for (ui v : out_neighbors) // clear the cache
				{
					vis[v] = 0;
					for (ui w : g.vertices[v].out_neighbors)
						vis[w] = 0;
				}
				if (ok)
				{
					q.push(u);
					in_queue[u] = 1;
				}
			}
		}
		g.remove_given_vertices(in_queue);
	}
	/**
	 * enumerate the subsets of candidate of size k-1 by depth-first-searching
	 */
	void dfs(vector<ui> &candidate, ui start_of_candidate, vector<ui> &plex, ui current_plex_size)
	{
		if (current_plex_size + 1 == plex.size())
		{
			for (ui i = start_of_candidate; i < candidate.size(); i++)
			{
				plex[current_plex_size] = candidate[i];
				out.dump(plex);
			}
			return;
		}
		for (ui i = start_of_candidate; i < candidate.size(); i++)
		{
			plex[current_plex_size] = candidate[i];
			dfs(candidate, i + 1, plex, current_plex_size + 1);
		}
	}
	/**
	 * @brief for u in V, if din[u]=0 or dout[u]=0, then we can directly dump all the maximal DPlexes containing u
	 */
	void no_out_or_no_in()
	{
		if (k != l)
		{
			no_out_neighbors();
			// we ignore din[u]=0 since (k<L) we can not guarantee any subset of size L is a DPlex
			return;
		}
		queue<ui> q;
		ui n = g.n;
		vector<bool> in_que(n);
		for (auto &vertex : g.vertices)
		{
			if (!vertex.dout || !vertex.din)
			{
				q.push(vertex.id);
				in_que[vertex.id] = 1;
			}
		}
		if (!q.size())
			return;
		// if i in V: candidate[index[i]]=i; if i not in V, then index[i] > n
		// i.e., candidate is the vertex set of the not removed vertices
		vector<ui> candidate(n), index(n);
		for (ui i = 0; i < n; i++)
		{
			candidate[i] = index[i] = i;
		}
		while (q.size() && k <= n)
		{
			ui u = q.front();
			q.pop();
			// remove u from candidate
			{
				ui pos = index[u];
				assert(pos < n);
				assert(candidate[pos] == u);
				index[candidate[n - 1]] = pos;
				index[u] = n;
				swap(candidate[pos], candidate[n - 1]);
				candidate.pop_back();
				n--;
			}
			// now generate DPlex containing u using k-1 vertices in candidate[0...n-1]
			// enumerate all the subsets of size k-1 using depth-first-search (dfs)
			vector<ui> plex(k);
			plex[0] = u;
			dfs(candidate, 0, plex, 1);
			// remove u from g
			for (ui v : g.vertices[u].out_neighbors)
			{
				if (!in_que[v])
				{
					if (--g.vertices[v].din == 0)
					{
						q.push(v);
						in_que[v] = 1;
					}
				}
			}
			for (ui v : g.vertices[u].in_neighbors)
			{
				if (!in_que[v])
				{
					if (--g.vertices[v].dout == 0)
					{
						q.push(v);
						in_que[v] = 1;
					}
				}
			}
		}
		if (k > n)
		{
			// there is no maximal DPlex in g now
			g.n = g.m = 0;
			g.vertices.clear();
		}
		else
		{
			// rebuild the graph using the rest vertices
			if (g.n > n)
			{
				assert(candidate.size() == n);
				vector<int> vis(g.n, -1);
				vector<ui> &S = candidate;
				g = Graph(S, vis, g);
			}
		}
	}

	/**
	 * @brief for u in V, if dout[u]=0, then we can directly dump all the maximal DPlexes containing u
	 */
	void no_out_neighbors()
	{
		queue<ui> q;
		ui n = g.n;
		vector<bool> in_que(n);
		for (auto &vertex : g.vertices)
		{
			if (!vertex.dout)
			{
				q.push(vertex.id);
				in_que[vertex.id] = 1;
			}
		}
		if (!q.size())
			return;
		// if i in V: candidate[index[i]]=i; if i not in V, then index[i] > n
		// i.e., candidate is the vertex set of the not removed vertices
		vector<ui> candidate(n), index(n);
		for (ui i = 0; i < n; i++)
		{
			candidate[i] = index[i] = i;
		}
		while (q.size() && k <= n)
		{
			ui u = q.front();
			q.pop();
			// remove u from candidate
			{
				ui pos = index[u];
				assert(pos < n);
				assert(candidate[pos] == u);
				index[candidate[n - 1]] = pos;
				index[u] = n;
				swap(candidate[pos], candidate[n - 1]);
				candidate.pop_back();
			}
			n--;
			// now generate DPlex containing u using k-1 vertices in candidate[0...n-1]
			// enumerate all the subsets of size k-1 using depth-first-search (dfs)
			vector<ui> plex(k);
			plex[0] = u;
			dfs(candidate, 0, plex, 1);
			// remove u from g, and we just consider updating dout[]
			for (ui v : g.vertices[u].in_neighbors)
			{
				if (!in_que[v])
				{
					if (--g.vertices[v].dout == 0)
					{
						q.push(v);
						in_que[v] = 1;
					}
				}
			}
		}
		if (k > n)
		{
			// there is no maximal DPlex in g now
			g.n = g.m = 0;
			g.vertices.clear();
		}
		else
		{
			// rebuild the graph using the rest vertices
			if (g.n > n)
			{
				assert(candidate.size() == n);
				vector<int> vis(g.n, -1);
				vector<ui> &S = candidate;
				g = Graph(S, vis, g);
			}
		}
	}
};

#endif