#include "graph.h"
#include <ctime>

int main(int argc, char *argv[])
{
	string graphname, graphfile;
	I K, L, threshold;
	graphname=argv[1];
	K=atoi(argv[2]);
	L=atoi(argv[3]);
	threshold=atoi(argv[4]);

	cout << graphname << "  K=" << K << "  L=" << L << "  Theta=" << threshold << endl;
	graphfile="../../sample_dataset/" + graphname + ".txt";

	Graph g;

	g.K = K;
	g.L = L;
	g.threshold = threshold;

	g.dataInput(graphfile);

	auto start = chrono::steady_clock::now();
	g.GP_Core();

	string s = "-";
	try{
		g.OPT();
	}
	catch (string s1){
		s = s1;
	}
	auto end = chrono::steady_clock::now();
	double RuningTime = chrono::duration_cast<chrono::microseconds>(end - start).count();

	cout << "precent_delE = " << float(g.del_te)/float(g.m) << endl;

	cout << "RuningTime : " << RuningTime/1000000 << endl;

	return 0;
}