#include "graph.h"

VI operator+(VI &v1, VI &v2)
{
	VI res(v1.size() + v2.size());
	int sz = set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
	res.resize(sz);
	return res;
};

VI operator-(VI &v1, VI &v2)
{
	VI res(v1.size());
	int sz = set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
	res.resize(sz);
	return res;
};

VI operator*(VI &v1, VI &v2)
{
	VI res(v1.size());
	int sz = set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
	res.resize(sz);
	return res;
};

void Graph::dataInput(string graphname)
{
	ifstream infile;
	infile.open(graphname);
	infile >> n >> m >> T;

	Node.reserve(n);
	U.reserve(n);
	Selected.resize(n, 0);
	VI emp(T, 0);
	TSelected.resize(n, emp);

	for (I i = 0; i < n; i++)
	{
		U.push_back(i);

		Node1 n1;
		n1.Deg = 0;
		n1.Sup = 0;
		n1.exist = 1;

		VI emp = {};
		for (auto j = 0; j < T; j++)
		{
			n1.T_adj.push_back(emp);
			n1.T_adjedge.push_back(emp);
			n1.T_Deg.push_back(0);
			n1.T_exist.push_back(0);
		}

		Node.push_back(n1);
	}

	I u, v, t;
	I pu = -1, pv = -1;
	k=2*K;
	for (I j = 0; j < m; j++)
	{
		infile >> u >> v >> t;

		Edge1 e;
		e.u = u;
		e.v = v;
		e.t = t;
		e.exist = 1;
		Edge.push_back(e);

		Node[u].Deg++;
		Node[v].Deg++;

		Node[u].T_adj[t].push_back(v);
		Node[u].T_adjedge[t].push_back(j);
		Node[u].T_Deg[t]++;
		Node[u].T_exist[t] = 1;

		Node[v].T_adj[t].push_back(u);
		Node[v].T_adjedge[t].push_back(j);
		Node[v].T_Deg[t]++;
		Node[v].T_exist[t] = 1;

		if (u != pu || v != pv)
		{
			Node[u].adj.push_back(v);
			Node[v].adj.push_back(u);
			Node[u].nadj.insert(v);
			Node[v].nadj.insert(u);
		}
		pu = u;
		pv = v;
	}

	infile.close();
}

int Graph::sest(I u)
{
	unordered_set<int> cl;
    cl.insert(u);

    vector<int> candidates(Node[u].nadj.begin(), Node[u].nadj.end());

    while (!candidates.empty()) {
        int best = -1;
        int maxCommon = -1;

        for (int v : candidates) {
            int common = 0;
            for (int c : cl) {
                if (Node[v].nadj.count(c)) ++common;
            }
            if (common == cl.size() && common > maxCommon) {
                maxCommon = common;
                best = v;
            }
        }

        if (best == -1) break;
        cl.insert(best);

        vector<int> newCandidates;
        for (int v : candidates) {
            if (v == best) continue;
            bool valid = true;
            for (int c : cl) {
                if (!Node[v].nadj.count(c)) {
                    valid = false;
                    break;
                }
            }
            if (valid) newCandidates.push_back(v);
        }
        candidates = std::move(newCandidates);
    }

    return cl.size();
}

void Graph::est()
{
	I maxD = 0, umax = -1;
	for (auto u : U)
	{
		if (Node[u].Deg > maxD)
		{
			maxD = Node[u].Deg;
			umax = u;
		}
	}
	threshold = max(sest(umax), 2*k-1);
}

void Graph::GP_Core()
{
	for (I t = 0; t < T; t++)
	{
		for (auto u : U)
		{
			if (Node[u].T_Deg[t] > 0)
				Node[u].Sup++;
			else
				Node[u].T_exist[t] = 0;
		}
	}

	for (I t = 0; t < T; t++)
	{
		for (auto u : U)
		{
			if ((Node[u].T_Deg[t] > 0 && Node[u].T_Deg[t] < threshold - K) || (Node[u].Sup > 0 && Node[u].Sup < L))
				sub_GP_Core(u, t);
		}
	}

	for (auto u : U)
	{
		if (Node[u].Sup == 0)
			Node[u].exist = 0;
	}

	I new_m = 0;

	for (I t = 0; t < T; t++)
	{
		int _m = 0;
		for (auto u : U)
		{
			if (Node[u].T_exist[t] == 0)
			{
				del_tn++;
				Node[u].T_adj[t].clear();
				Node[u].T_adjedge[t].clear();
				continue;
			}
			I p = 0;

			for (I i = 0; i < Node[u].T_adj[t].size(); i++)
			{
				vI v = Node[u].T_adj[t][i];
				if (Node[v].T_exist[t])
					Node[u].T_adj[t][p++] = v;
			}
			Node[u].T_adj[t].resize(p);
			Node[u].T_Deg[t] = p;
			Node[u].Deg += p;
			_m += p;
		}
		new_m += _m;
		assert(_m % 2 == 0);
	}
	new_m = new_m/2;
	del_te = m - new_m;

	VI del_U = {};
	del_U.reserve(n);
	for (auto u : U)
	{
		if (Node[u].exist == 0)
		{
			del_U.push_back(u);
			del_n++;
			Node[u].adj.clear();
			continue;
		}

		I p = 0;
		for (I i = 0; i < Node[u].adj.size(); i++)
		{
			vI v = Node[u].adj[i];
			if (Node[v].exist)
				Node[u].adj[p++] = v;
		}
		Node[u].adj.resize(p);
	}
	
	U = U - del_U;
}

void Graph::sub_GP_Core(uI u, I t)
{
	Node[u].T_Deg[t] = 0;
	Node[u].T_exist[t] = 0;

	for (auto v : Node[u].T_adj[t])
	{
		if (Node[v].T_Deg[t] > 0)
		{
			Node[v].T_Deg[t]--;
			if (Node[v].T_Deg[t] < threshold - K)
				sub_GP_Core(v, t);
		}
	}

	if (Node[u].Sup > 0)
	{
		Node[u].Sup--;
		if (Node[u].Sup < L)
		{
			Node[u].Sup = 0;
			Node[u].exist = 0;
			for (I t = 0; t < T; t++)
			{
				if (Node[u].T_Deg[t] > 0)
					sub_GP_Core(u, t);
			}
		}
	}
}

void Graph::MUM()
{
	I U_size = U.size();

	ms = 2 * K - 1;

	vector<PII> order_U;
	for (I i = 0; i < U_size; i++)
		order_U.emplace_back(U[i], Node[U[i]].Deg);

	sort(order_U.begin(), order_U.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
		return a.second < b.second;
	});

	TDeg_inAns.resize(n);
	for (I i = 0; i < n; i++)
		TDeg_inAns[i].resize(T, 0);

	Cand.resize(U.size() + 2);
	Cand_t.resize(U.size() + 2);

	Ans.resize(U.size());

	inAns.resize(n, false);

	for (I i = 0; i < U_size; i++)
	{
		uI u = order_U[i].first;
		if (Node[u].exist == 0)
			continue;

		int_sel++;
		Selected[u] = int_sel;
		for (auto v : Node[u].adj)
		{
			if (Node[v].exist) Selected[v] = int_sel;
			for (auto w : Node[v].adj)
				if (Node[w].exist) Selected[w] = int_sel;
		}

		for (auto v : U)
		{
			Node[v].Sup = 0;
			if (Node[v].exist == 0 || Selected[v] != int_sel) continue;
			for (I t = 0; t < T; t++)
			{
				Node[v].T_Deg[t] = 0;
				if (Node[v].T_exist[t] == 0) continue;
				TSelected[v][t] = int_sel;
				for (auto w : Node[v].T_adj[t])
					if (Node[w].T_exist[t] && Selected[w] == int_sel)
						Node[v].T_Deg[t]++;
			}
		}

		in_GP_Core(ms);
		if (Selected[u] != int_sel)
		{
			Node[u].exist = 0;
			for (I t = 0; t < T; t++) Node[u].T_exist[t] = 0;
			continue;
		}

		Ans[0] = u;
		inAns[u] = true;
		Cand[1].clear();
		for (auto v : U)
			if (Selected[v] == int_sel && v != u) Cand[1].push_back(v);
		Cand_t[1].clear();
		for (I t = 0; t < T; t++)
			if (TSelected[u][t] == int_sel) Cand_t[1].push_back(t);
		for (auto t : Cand_t[1])
			for (auto v : Node[u].T_adj[t])
				if (TSelected[v][t] == int_sel) TDeg_inAns[v][t]++;

		sub_MUM(1);

		for (auto t : Cand_t[1])
			for (auto v : Node[u].T_adj[t])
				if (TSelected[v][t] == int_sel) TDeg_inAns[v][t]--;
		inAns[u] = false;

		Node[u].exist = 0;
		for (I t = 0; t < T; t++) Node[u].T_exist[t] = 0;
	}
}


void Graph::sub_MUM(I p)
{
	if ((double)clock() / CLOCKS_PER_SEC - throwtime > TIMEOVER)
		throw string("TIMEOVER");
	
	for (auto u : Cand[p])
	{
		Ans[p] = u;
		inAns[u] = true;

		DVB_Pruning(p);

		if (p + 1 + Cand[p+1].size() > ms && Cand_t[p+1].size() >= L)
		{
			sub_MUM(p + 1);

			if (p + 1 > ms)
			{
				ms = p + 1;
			}
		}

		for (auto t : Cand_t[p])
		{
			if (TDeg_inAns[u][t] < p + 1 - K)
				continue;

			for (auto v : Node[u].T_adj[t])
				if (TSelected[v][t] == int_sel)
					TDeg_inAns[v][t]--;
		}

		inAns[u] = false;
	}
}

void Graph::in_GP_Core(I lb)
{
	for (I t = 0; t < T; t++)
	{
		for (auto u : U)
		{
			if (Node[u].exist == 0 || Selected[u] != int_sel) continue;
			if (Node[u].T_Deg[t] > 0)
				Node[u].Sup++;
			else
				TSelected[u][t]--;
		}
	}

	for (I t = 0; t < T; t++)
	{
		for (auto u : U)
		{
			if (Selected[u] != int_sel || TSelected[u][t] != int_sel) continue;
			if ((Node[u].T_Deg[t] > 0 && Node[u].T_Deg[t] < lb - K) || (Node[u].Sup > 0 && Node[u].Sup < L))
				in_sub_GP_Core(u, t, lb);
		}
	}

	for (auto u : U)
	{
		if (Node[u].Sup == 0)
			Selected[u]--;
	}
}

void Graph::in_sub_GP_Core(uI u, I t, I lb)
{
	Node[u].T_Deg[t] = 0;
	TSelected[u][t]--;

	for (auto v : Node[u].T_adj[t])
	{
		if (Selected[v] == int_sel && TSelected[v][t] == int_sel && Node[v].T_Deg[t] > 0)
		{
			Node[v].T_Deg[t]--;
			if (Node[v].T_Deg[t] < lb - K)
				in_sub_GP_Core(v, t, lb);
		}
	}

	if (Node[u].Sup > 0)
	{
		Node[u].Sup--;
		if (Node[u].Sup < L)
		{
			Node[u].Sup = 0;
			Selected[u]--;
			for (I t = 0; t < T; t++)
			{
				if (TSelected[u][t] == int_sel && Node[u].T_Deg[t] > 0)
					in_sub_GP_Core(u, t, lb);
			}
		}
	}
}

void Graph::DVB_Pruning(I p)
{
	Cand[p+1].clear();
	Cand_t[p+1].clear();

	uI u = Ans[p];
	
	for (auto t : Cand_t[p])
	{
		if (TDeg_inAns[u][t] < p + 1 - K)
			continue;

		for (auto v : Node[u].T_adj[t])
			if (TSelected[v][t] == int_sel)
				TDeg_inAns[v][t]++;

		I flag = 0;
		for (I i = 0; i < p; i++)
		{
			vI v = Ans[i];
			if (TDeg_inAns[v][t] < p + 1 - K)
			{
				flag = 1;
				break;
			}
		}

		if (flag)
			continue;

		Cand_t[p+1].push_back(t);
	}

	if (Cand_t[p+1].size() < L)
		return;

	for (auto v : Cand[p])
	{
		if (v <= u)
			continue;

		I time_count = 0;
		for (auto t : Cand_t[p+1])
		{
			if (TDeg_inAns[v][t] >= p + 1 - K)
				time_count++;
		}
		if (time_count >= L)
			Cand[p+1].push_back(v);
	}
}