// https://www.optil.io/optilion/problem/3028
// tar zcvf AntColony.tgz AntColony.cpp CMakeLists.txt

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <functional>
#include <set>
#include <queue>
#include <signal.h>
#include <cstring>
#include <iostream>
#include <cmath>

volatile sig_atomic_t tle = 0;

#define INF 2000000007
//#define repNodes(i, n) for(int i = 1; i <= n && !tle; ++i)
//#define repEdges(i, n) for(int i = 0; i < n && !tle; ++i)

using namespace std;

typedef pair<int, int> P;

typedef unsigned long long ull;

/*****     paramater    ******/
int numEpisode;
double epsilon;
int resetPheromone = 0;
unsigned long long Q = 1;// = 100000;
int a = 5; // pheromone rate
int b = 3; // cost rate
		   //double gamma = 0.3;

class Ant;
class Graph;
class Vertex;
class Edge;

unsigned long long bestCost;
vector<Edge*> bestRoute;
vector<double> pheList;

//#define repV(it, V) for(vector<int>::iterator it = V.begin(), end = V.end(); it != end && !tle; ++it)
//#define repVecE(it, Route) for(vector<Edge*>::iterator it = Route.begin(), end = Route.end(); it != end && !tle; ++it)

/*************           input variable           ***************/
int numNodes = INF, numEdges = INF, numTerminals = INF;


void file_input();


// random
random_device seed_gen;
mt19937 engine(1);

/*************			  ant colony               ***************/
class Graph {
public:
	vector<Vertex> v;
	vector<Edge> e;
	vector<int> t;
};
class Vertex {
public:
	int id;
	vector<Edge*> e;
	//double pheromone;
	/*Vertex() {
		pheromone = 1000;
	}//*/
};
class Edge {
public:
	int id;
	Vertex *v1, *v2; // connected g.v
	double pheromone;
	int cost;
	Edge() {
		pheromone = 1;
	}
};
class Ant {
public:
	int id;
	vector<Edge*> tabu_list;  // visited node by this ant.
	Vertex *pos;             // now position of this ant.
	int index;  //tabu_list index
};

Graph g;

void handler(int signum){
	tle = 1;
}
void term(int signum) {
	if(signum == SIGTERM || signum == SIGINT){
		printf("VALUE %llu\n", bestCost);
		#ifdef DEBUG
		printf("Edges %lu\n",bestRoute.size());
		#else

		for (vector<Edge*>::iterator it = bestRoute.begin(),endl=bestRoute.end(); it != endl; ++it) {
			printf("%d %d\n", (*it)->v1->id, (*it)->v2->id);
		}
		#endif
		exit(0);
	}
}

/*************           input function           ***************/
void file_input() {

	char s1[10], s2[10];
	scanf("%s %s", s1, s2);
	//printf("%s %s\n", s1, s2);   //SECTION

	// Number of g.v
	scanf("%s %d", s1, &numNodes);

	pheList.resize(numEdges + 1);
	g.v.resize(numNodes + 1);
	for(int i = 1; i <= numNodes; ++i) g.v[i].id = i;
	
	// Number of Edges
	scanf("%s %d", s1, &numEdges);
	g.e.resize(numEdges + 1);
	
	// Edges
	int start, end, weight;
	for(int i = 0; i < numEdges; ++i) {
		scanf("%s %d %d %d", s1, &start, &end, &weight);
		g.e[i].id = i;
		g.e[i].cost = weight;
		//Q += weight;
		g.e[i].v1 = &g.v[start];                // connected vertexs of edge
		g.e[i].v2 = &g.v[end];
		g.v[start].e.push_back(&g.e[i]);   // connected edges of vertex
		g.v[end].e.push_back(&g.e[i]);
	}

	scanf("%s", s1);   // END
	scanf("%s %s", s1, s2);   // SECTION Terminals

	scanf("%s %d", s1, &numTerminals);   // Number of Terminals
										 //Q *= ((double)numTerminals / numEdges)*100;
										 // Terminals
	int terminal;
	for(int i = 1; i <= numTerminals; ++i) {
		scanf("%s %d", s1, &terminal);
		g.t.push_back(terminal);
	}
}

/*************               union find           ***************/
int find(vector<int> &T, int x) {
	//if (x < 0) return x;
	if (T[x] < 0) return x;
	return T[x] = find(T, T[x]);
}

void unite(vector<int> &T, int x, int y) {
	x = find(T, x);
	y = find(T, y);
	if (x == y) return;
	T[y] = x;
	T[x]--;
}

/*************             ant function           ***************/
bool antMove(Ant *ant, vector<int> &passedAnt, vector<int> &T) {
	Vertex *v = ant->pos;
	vector<Edge*> e = v->e;    // possible edge of this ant
	int moveEdges = e.size();

	// already passed edge was deleted
	
	for(int i = 0; i < moveEdges; ++i) {
		if (passedAnt[e[i]->v1->id] != -1 && passedAnt[e[i]->v2->id] != -1) {
			int v1Ant = find(T, passedAnt[e[i]->v1->id]), v2Ant = find(T, passedAnt[e[i]->v2->id]);  // find parent
			if (v1Ant == ant->id && v2Ant == ant->id) {   // already passed edge by myself
				e[i] = e.back();
				e.pop_back();
				moveEdges--;
				i--;
			}
		}
	}
	

	if (moveEdges == 0) {
		if (ant->index < 0) ant->index = ant->tabu_list.size() - 1;
		Edge *beforeEdge = ant->tabu_list[ant->index--];
		//beforeEdge->pheromone = 1e-3;
		if (ant->pos == beforeEdge->v1) ant->pos = beforeEdge->v2;
		else ant->pos = beforeEdge->v1;

		return true;
	}
	
	// calculate probability
	double sumPheromone = 0, sumCost = 0;
	vector<double> addPheromone(moveEdges);
	fill(addPheromone.begin(), addPheromone.end(), 0);
	double maxPheromone = e[0]->pheromone;
	for(int i = 1; i < moveEdges; ++i) {
		//double pherom;
		//if(e[i]->v1 == ant->pos) pherom = e[i]->v2->pheromone;
		//else pherom = e[i]->v1->pheromone;

		//if (pherom > maxPheromone) maxPheromone = e[i]->pheromone;
		maxPheromone = max(maxPheromone,e[i]->pheromone);
	}

	for(int i = 0; i < moveEdges; ++i) {
		if (passedAnt[e[i]->v1->id] != -1 && passedAnt[e[i]->v2->id] != -1) { // already ant passed
//			double pherom;
			//if(e[i]->v1 == ant->pos) pherom = e[i]->v2->pheromone;
			//else pherom = e[i]->v1->pheromone;
			addPheromone[i] = maxPheromone;//e[i]->pheromone;
		}

		//double pherom;
		//if(e[i]->v1 == ant->pos) pherom = e[i]->v2->pheromone;
		//else pherom = e[i]->v1->pheromone;

		//if (i==0 || pherom > maxPheromone) maxPheromone = pherom;
		if (e[i]->pheromone < 0.1) {
			e[i]->pheromone = 1;// 0.1;
		}
		if (e[i]->pheromone > 1e10) {
			e[i]->pheromone = 1;// 1e10;
		}

		sumPheromone += pow(e[i]->pheromone + addPheromone[i], a);    // sum of pheromons
		sumCost += pow(e[i]->cost + 0.0001 / moveEdges, b);  // sum of cost
	}
	
	// probability
	vector<double> p(moveEdges);
	long double sumP = 0;

	for(int i = 0; i < moveEdges; ++i) {
		//double pherom;
		//if(e[i]->v1 == ant->pos) pherom = e[i]->v2->pheromone;
		//else pherom = e[i]->v1->pheromone;

		p[i] = pow(e[i]->pheromone + addPheromone[i], a) / sumPheromone;
		p[i] *= sumCost / pow(e[i]->cost + 0.0001, b);
		sumP += p[i];
	}

	// uniform distribution
	uniform_real_distribution<> dist1(0, sumP);
	double r = dist1(engine);
	// move
	for(int i = 0; i < moveEdges; ++i) {
		r -= p[i];
		if (r <= 0) {
			// update pos
			if (ant->pos == e[i]->v1) ant->pos = e[i]->v2;
			else if (ant->pos == e[i]->v2) ant->pos = e[i]->v1;
			// update tabu-list
			ant->tabu_list.push_back(e[i]);
			ant->index = ant->tabu_list.size() - 1;  //index must be initilized

			break;
		}
	}
	

	return false;
}

typedef pair<bool, ull> BP;

BP branchCut(vector<Edge*> &Route) {
	vector<int> onlyVertex(g.v.size());
	fill(onlyVertex.begin(), onlyVertex.end(), -1);
	int cnt = 0;
	ull res = 0;
	for(vector<Edge*>::iterator it = Route.begin(), end = Route.end(); it != end; ++it) {
		int v1 = (*it)->v1->id;
		int v2 = (*it)->v2->id;

		// v1
		if (onlyVertex[v1] == -1)     onlyVertex[v1] = cnt;  // anyone not passed
		else if (onlyVertex[v1] >= 0) onlyVertex[v1] = -2; // already passed

		// v2
		if (onlyVertex[v2] == -1)     onlyVertex[v2] = cnt;
		else if (onlyVertex[v2] >= 0) onlyVertex[v2] = -2;
		cnt++;
		res += (*it)->cost;
	}

	// rest terminal
	for(vector<int>::iterator it = g.t.begin(), end = g.t.end(); it != end; ++it) onlyVertex[*it] = -1;

	sort(onlyVertex.begin(), onlyVertex.end(), greater<int>());  // because delete back of list

	bool boo = false;
	for(vector<int>::iterator it = onlyVertex.begin(), end = onlyVertex.end(); it != end; ++it) {
		if (*it < 0) break;
		boo = true;
		res -= Route[*it]->cost;
		Route[*it] = Route.back();
		Route.pop_back();
	}
	BP ans;
	ans = BP(boo, res);
	return ans;
}
BP branchCut_Tmp(vector<Edge*> &Route){
	vector<vector<int> >	nextList(g.v.size()+1);		// 点を含有する辺を保持
	vector<int>				degree(g.v.size()+1,0);			// 点の次数
	int allcost = 0;
	vector<int> 			deleteEdges;

	// init
	//fill(degree.begin(),degree.end(),0);

	for(int i=0,endlnum = Route.size();i<endlnum;++i){
	//for (vector<Edge*>::iterator it = Route.begin(),endl=Route.end(); it != endl; ++it) {
		int v1 = Route[i]->v1->id;
		int v2 = Route[i]->v2->id;

		++degree[v1],++degree[v2];		// 次数を加算する
		// 含有する辺番号を追加
		nextList[v1].push_back(i);
		nextList[v2].push_back(i);
		// コストを加算
		allcost += Route[i]->cost;
	}

	// ターミナル点を区別させる
	for (vector<int>::iterator it = g.t.begin(),endl=g.t.end(); it != endl; ++it) degree[*it] = -1;

	for(int i=1,endlnum=g.v.size();i<=endlnum;++i){
		if(degree[i]==1){		// 次数が１==取り除ける点
			int p = i;			// 注目している点
			int next;			// 点p の次の点
	//		cout << p << ":"<<nextList[p].size()<<"=";
			while (degree[p]==1){
				next = Route[nextList[p][0]]->v1->id == p ? Route[nextList[p][0]]->v2->id : Route[nextList[p][0]]->v1->id;
	//			cout << next << "[" << degree[next] << "]" << nextList[p].size() << ":";
				allcost -= Route[nextList[p][0]]->cost;		// Cost 
				--degree[p],--degree[next];					// Degree
				// Delete 
				for(int j=0;j<nextList[next].size();++j){
					if(Route[nextList[next][j]]->v1->id == p || Route[nextList[next][j]]->v2->id == p){
	//					cout << Route[nextList[next][j]]->v1->id <<":"<<Route[nextList[next][j]]->v2->id<<" ";
						nextList[next][j] = nextList[next].back();
						nextList[next].pop_back();
						break;
					}
				}
				deleteEdges.push_back(nextList[p][0]);
				p = next;
			}
		}
	}
	
	sort(deleteEdges.begin(), deleteEdges.end(), greater<int>());		// sort
	for (vector<int>::iterator it = deleteEdges.begin(),endl=deleteEdges.end(); it != endl; ++it){
		Route[*it] = Route.back();
		Route.pop_back();
	}

	return BP(false,allcost);
}
typedef pair<int, Edge*> PE;
typedef pair<int, PE> PPE;
typedef pair<ull, vector<Edge*> > PVE;

PVE prim(vector<Edge*> Route) {

	set<int> st;  // Duplication
	for(vector<Edge*>::iterator it = Route.begin(), end = Route.end(); it != end; ++it){
		int v1 = (*it)->v1->id, v2 = (*it)->v2->id;   // add vertex
		st.insert(v1);
		st.insert(v2);
	}

	vector<bool> used(numNodes + 1);       // already used
	vector<int> d(numNodes + 1);   // most minimum distance of v
	fill(d.begin(), d.end(), INF);
	fill(used.begin(), used.end(), false);

	priority_queue<PPE, vector<PPE>, greater<PPE> > Q;
	set<int>::iterator it = st.begin(); // begin vertex
	d[*it] = 0;
	Q.push(PPE(0, PE(*it, nullptr)));  // (cost, Vertex Number) start 

	vector<Edge*> primRoute;
	while (!Q.empty()) {
		PPE pp = Q.top(); Q.pop();
		int v = pp.second.first;
		if (used[v]) continue;   // already used Vertex
		used[v] = true;
		primRoute.push_back(pp.second.second);

		for(vector<Edge*>::iterator it = g.v[v].e.begin(), end = g.v[v].e.end(); it != end; ++it) {// edges of v
			Edge *e = *it;
			set<int>::iterator first = st.find(e->v1->id);
			if (first != st.end() && d[e->v1->id] > e->cost) {
				d[e->v1->id] = e->cost;
				Q.push(PPE(e->cost, PE(e->v1->id, e)));
			}
			first = st.find(e->v2->id);
			if (first != st.end() && d[e->v2->id] > e->cost) {
				d[e->v2->id] = e->cost;
				Q.push(PPE(e->cost, PE(e->v2->id, e)));
			}
		}
	} // end prim
	primRoute[0] = primRoute.back();  // delete nullptr
	primRoute.pop_back();

	// branch cut
	ull ansSum = INF;
	vector<Edge*> ansRoute = primRoute, beforeRoute = primRoute;
	//vector<Edge*> beforeRoute = Route;
	//bool boo = true;
	//while (boo) {
	BP bcut = branchCut_Tmp(beforeRoute);
		//boo = bcut.first;
	if (bcut.second < ansSum) {
		ansRoute = beforeRoute;
		ansSum = bcut.second;
	}

	return PVE(ansSum, ansRoute);
}

void updatePheromone(vector<Edge*> allRoute, int episode) {

	/******           prim             *****/
	vector<Edge*> resRoute;

	PVE pve = prim(allRoute);

	ull sum = pve.first;
	resRoute = pve.second;
	//cout << " sum = " << sum;
	if (episode == 0) {
		bestCost = sum;
		bestRoute = resRoute;//allRoute;
		Q = bestCost;
		//for(int i = 1; i <= numEdges; ++i) pheList[i] = g.e[i].pheromone;
	}
	else {
		if (sum < bestCost) { // save the best route and cost
			bestCost = sum;
			bestRoute = resRoute; //allRoute;
		}

		// update pheromone
		double tau = Q / (double)sum+1;
		double rho = 0.3;

		set<int> st;  // Duplication
		for (vector<Edge*>::iterator it = resRoute.begin(); it != resRoute.end(); it++) {
			st.insert((*it)->v1->id);
			st.insert((*it)->v2->id);
		}

		for(set<int>::iterator it = st.begin(); it != st.end(); it++)
			pheList[*it] = (1 - rho) * pheList[*it] + tau;
		

		if(episode % 50 == 0)
			for(int i = 1; i <= numEdges; ++i) g.e[i].pheromone = pheList[i];
	#ifdef DEBUG
		printf("\tCost : %d\n",sum);
	#endif
	}

}

void deathGame(int episode) {
	/***    Ant variable   ***/
	vector<Ant> ants(numTerminals);  // initial ants
	int numAnts = numTerminals;      // number of ants
	vector<int> passedAnt(g.v.size());    // insert already passed ant's ID
	vector<Ant*> aliveAnt(numTerminals);  // not died ant ID
	vector<int> T(numTerminals);          // union find tree

	fill(T.begin(), T.end(), -1);
	fill(passedAnt.begin(), passedAnt.end(), -1);

	// Ant initial position
	for (int i = 0; i < numTerminals; i++) {
		ants[i].id = i;
		ants[i].pos = &g.v[g.t[i]];
		passedAnt[g.t[i]] = i;
		aliveAnt[i] = &ants[i];
	}

	while (numAnts > 1) {
		if(tle) return;
		//for (int i = 0; i < numAnts; i++) {
			int k = rand() % numAnts;
			Ant *nowAnt = aliveAnt[k];

			/****    move   ****/
			uniform_real_distribution<> dist2(0, 1);   // uniform distribution
			double r = dist2(engine);
			bool backEdge = false;
			backEdge = antMove(nowAnt, passedAnt, T);


			// ant[nowAntId] die.
			if (passedAnt[nowAnt->pos->id] != -1 && !backEdge) {

				int parentID = find(T, passedAnt[nowAnt->pos->id]);  // alive ant
				if (nowAnt->id != parentID) {  // not my child
					Ant *parentAnt = &ants[parentID];
					unite(T, parentID, nowAnt->id);

					// merge
					parentAnt->tabu_list.reserve(parentAnt->tabu_list.size() + nowAnt->tabu_list.size());
					copy(nowAnt->tabu_list.begin(), nowAnt->tabu_list.end(), back_inserter(parentAnt->tabu_list)); // merge

					aliveAnt[k] = aliveAnt.back();
					aliveAnt.pop_back();
			//		i--;
					numAnts--;
				}
			}
			// no die
			else if (!backEdge) {
				passedAnt[nowAnt->pos->id] = nowAnt->id;
			}
			if (numAnts == 1) break;
		//}
	}
	updatePheromone(aliveAnt[0]->tabu_list, episode);    // update pheromone
}

int main() {
	struct sigaction action;
	memset(&action, 0, sizeof(struct sigaction));
	action.sa_handler = handler;
	sigaction(SIGTERM, &action, NULL);
	signal(SIGINT,handler);

	file_input();
	
	for (int episode = 0; !tle; ++episode) {
		deathGame(episode); 		  // deathGame
		#ifdef DEBUG
			printf("episode = %d ", episode);
			printf(" bestCost = %llu\n", bestCost);
		#endif
	}

	printf("VALUE %llu\n", bestCost);
	#ifdef DEBUG
	printf("SIZE %lu\n",bestRoute.size());
	#else
	for (vector<Edge*>::iterator it = bestRoute.begin(); it != bestRoute.end(); ++it) {
		printf("%d %d\n", (*it)->v1->id, (*it)->v2->id);
	}
	#endif

	return 0;
}

