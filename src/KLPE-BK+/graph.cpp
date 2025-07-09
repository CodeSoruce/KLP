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

	for (I i = 0; i < n; i++)
	{
		U.push_back(i);

		Node1 n1;
		n1.Sup = 0;
		n1.exist = 1;

		VI emp = {};
		for (auto j = 0; j < T; j++)
			n1.T_adj.push_back(emp);

		for (auto j = 0; j < T; j++)
			n1.T_adjedge.push_back(emp);

		for (auto j = 0; j < T; j++)
			n1.T_Deg.push_back(0);

		for (auto j = 0; j < T; j++)
			n1.T_exist.push_back(0);

		Node.push_back(n1);
	}

	I u, v, t;
	for (I j = 0; j < m; j++)
	{
		infile >> u >> v >> t;

		Edge1 e;
		e.u = u;
		e.v = v;
		e.t = t;
		e.exist = 1;
		Edge.push_back(e);

		Node[u].T_adj[t].push_back(v);
		Node[u].T_adjedge[t].push_back(j);
		Node[u].T_Deg[t]++;
		Node[u].T_exist[t] = 1;

		Node[v].T_adj[t].push_back(u);
		Node[v].T_adjedge[t].push_back(j);
		Node[v].T_Deg[t]++;
		Node[v].T_exist[t] = 1;
	}

	infile.close();
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
		}
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

void Graph::Base(VI &U, VI &C, VI &X)
{

	if ((double)clock() / CLOCKS_PER_SEC - throwtime > TIMEOVER)
		throw string("TIMEOVER");


	if (U.size() + C.size() < threshold)
	{
		return;
	}

	VI doing(C);
	VI X1(X);

	I flag = 0;
	I Csize = C.size();

	for (auto u : doing)
	{
		unit[0] = u;
		C = C - unit;
		VI U_1 = U + unit;

		if (CheckFP(U_1))
		{
			flag = 1;
			VI C_1(C);
			Base(U_1, C_1, X);
			unit[0] = u;
			X = X + unit;
		}
	}

	if ((flag == 0 || Csize == 0) && U.size() >= threshold)
	{
		I in_flag = 1;
		for (auto v : X1)
		{
			unit[0] = v;
			VI U_1 = U + unit;
			if (CheckFP(U_1))
			{
				in_flag=0;
			}
		}
		if (in_flag)
			ms++;
	}
}

bool Graph::CheckFP(VI &U)
{
	I Usize = U.size();
	I con_t = 0;

	for (I t = 0; t < T; t++)
	{
		I flag = 1;
		for (auto u : U)
		{
			if (Usize - (Node[u].T_adj[t] * U).size() > K)
			{
				flag = 0;
			}
		}

		if (flag)
		{
			con_t++;
			if (con_t >= L)
				return true;
		}
	}

	return false;
}


