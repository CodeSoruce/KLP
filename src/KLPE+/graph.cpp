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
	ifstream infile(graphname);
	infile >> n >> m >> T;

	Node.resize(n);
	U.resize(n);
	for (I i = 0; i < n; i++) {
		U[i] = i;
		auto &node = Node[i];
		node.Sup = 0;
		node.exist = 1;
		node.T_adj.resize(T);
		node.T_adjedge.resize(T);
		node.T_Deg.assign(T, 0);
		node.T_exist.assign(T, 0);
	}

	Edge.resize(m);
	I u, v, t;
	for (I j = 0; j < m; j++) {
		infile >> u >> v >> t;
		Edge[j] = {u, v, t, 1};

		auto &nu = Node[u];
		auto &nv = Node[v];

		nu.T_adj[t].push_back(v);
		nv.T_adj[t].push_back(u);
		nu.T_adjedge[t].push_back(j);
		nv.T_adjedge[t].push_back(j);
		nu.T_Deg[t]++;
		nv.T_Deg[t]++;
		nu.T_exist[t] = 1;
		nv.T_exist[t] = 1;

		if (nu.adj.empty() || nu.adj.back() != v)
			nu.adj.push_back(v);
		if (nv.adj.empty() || nv.adj.back() != u)
			nv.adj.push_back(u);
	}
	infile.close();
}

void Graph::GP_Core()
{
    for (I t = 0; t < T; t++) {
        for (auto u : U) {
            auto &nu = Node[u];
            if (nu.T_Deg[t] > 0)
                nu.Sup++;
            else
                nu.T_exist[t] = 0;
        }
    }

    for (I t = 0; t < T; t++) {
        for (auto u : U) {
            auto &nu = Node[u];
            if ((nu.T_Deg[t] > 0 && nu.T_Deg[t] < threshold - K) || (nu.Sup > 0 && nu.Sup < L)) {
                sub_GP_Core(u, t);
            }
        }
    }

    for (auto u : U) {
        if (Node[u].Sup == 0)
            Node[u].exist = 0;
    }

    I new_m = 0;
    for (I t = 0; t < T; t++) {
        I edge_count = 0;
        for (auto u : U) {
            auto &nu = Node[u];
            if (!nu.T_exist[t]) {
                del_tn++;
                nu.T_adj[t].clear();
                nu.T_adjedge[t].clear();
                continue;
            }
            I p = 0;
            for (I i = 0; i < nu.T_adj[t].size(); i++) {
                int v = nu.T_adj[t][i];
                if (Node[v].T_exist[t]) {
                    nu.T_adj[t][p++] = v;
                }
            }
            nu.T_adj[t].resize(p);
            edge_count += p;
        }
        new_m += edge_count;
        assert(edge_count % 2 == 0);
    }
    new_m /= 2;
    del_te = m - new_m;

    I write_ptr = 0;
    for (I i = 0; i < U.size(); ++i) {
        int u = U[i];
        auto &nu = Node[u];
        if (!nu.exist) {
            del_n++;
            nu.adj.clear();
            continue;
        }

        U[write_ptr++] = u;

        I p = 0;
        for (I j = 0; j < nu.adj.size(); ++j) {
            int v = nu.adj[j];
            if (Node[v].exist)
                nu.adj[p++] = v;
        }
        nu.adj.resize(p);
    }
    U.resize(write_ptr);
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

void Graph::OPT()
{
	TDeg_inAns.resize(n);
	for (I i = 0; i < n; i++)
		TDeg_inAns[i].resize(T, 0);

	Hop_adj.resize(n);

	Cand.resize(U.size()+2);
	Xand.resize(U.size()+2);
	Cand_t.resize(U.size()+2);

	Ans.resize(U.size());

	inAns.resize(n, false);

	Cand[0] = U;

	for (auto u : U)
	{
		Ans[0] = u;
		inAns[u] = true;

		Cand[1].clear();
		Cand[1].push_back(u);
		for (auto v : Node[u].adj)
		{
			if (v < u)
				continue;
			Cand[1].push_back(v);
			for (auto w : Node[v].adj)
			{
				if (w < u)
					continue;
				Cand[1].push_back(w);
			}
		}
		sort(Cand[1].begin(), Cand[1].end());
		I p = 0;
		I pv = -1;
		for (I i = 0; i < Cand[1].size(); i++)
		{
			uI v = Cand[1][i];
			if (pv != v)
			{
				Cand[1][p++] = v;
				pv = v;
			}
		}
		Cand[1].resize(p);

		Cand_t[1].clear();
		for (I t = 0; t < T; t++)
		{
			if (Node[u].T_exist[t])
				Cand_t[1].push_back(t);
		}

		unit[0] = u;
		Cand[1] = Cand[1] - unit;
		Hop_adj[u] = Cand[1];

		for (auto t : Cand_t[1])
		{
			for (auto v : Node[u].T_adj[t])
			{
				TDeg_inAns[v][t]++;
			}
		}
		

		sub_OPT(1);

		for (auto t : Cand_t[1])
		{
			for (auto v : Node[u].T_adj[t])
			{
				TDeg_inAns[v][t]--;
			}
		}
		inAns[u] = false;
	}
}





I Graph::sub_OPT(I p)
{
	I flag = 0;
	for (auto u : Cand[p])
	{
		Ans[p] = u;
		inAns[u] = true;

		Processing(p);


		if (p + 1 + Cand[p+1].size() >= threshold && Cand_t[p+1].size() >= L)
		{
			I in_flag = sub_OPT(p+1);

			flag = 1;
			if (p + 1 >= threshold)
			{
				if (in_flag == 0 && Xand[p+1].size() == 0)
				{
					ms++;	
				}
			}
		}

		for (auto t : Cand_t[p])
		{
			if (TDeg_inAns[u][t] < p + 1 - K)
				continue;

			for (auto v : Node[u].T_adj[t])
			{
				TDeg_inAns[v][t]--;
			}
		}

		inAns[u] = false;
	}

	return flag;
}



void Graph::Processing(I p)
{
	Cand[p+1].clear();
	Xand[p+1].clear();
	Cand_t[p+1].clear();

	uI u = Ans[p];
	
	for (auto t : Cand_t[p])
	{
		if (TDeg_inAns[u][t] < p + 1 - K)
			continue;

		for (auto v : Node[u].T_adj[t])
		{
			TDeg_inAns[v][t]++;
		}

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

	for (auto v : Cand[1])
	{
		if (v >= u)
			break;

		if (inAns[v])
			continue;

		I time_count = 0;
		for (auto t : Cand_t[p+1])
		{
			if (TDeg_inAns[v][t] >= p + 1 - K)
				time_count++;
		}
		if (time_count >= L)
			Xand[p+1].push_back(v);
	}
}

