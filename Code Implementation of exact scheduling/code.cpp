//######################################//
// Author: Sid Mahapatra                //
// Last Modified: 03 September 2020		//					
// Description: Code for an exact method// 
// for scheduling of the alternative    //
// technologies in R&D projects (paper) //
//######################################//

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <bits/stdc++.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <string>
#include <stack>



#define DBG 1   // set DBG 1 for debugging code and 0 for normal run
#define R 0.05
#define PTC 8

using namespace std;

struct alternative {
	int duration;
	int es, ef, ls, lf, st;  // es : earliest stasrt time , ef : earliest finish time
							 // ls : latest start time ,  lf : latest finish time
							 // st : slack time 
	int nos;
	int cost;
	double p;
};

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}


// returns vector of n numbers for input
std::vector<int> ReadNumbers()
{
    std::vector<int> numbers ;

    do
    {
        int input ;
        if ( std::cin >> input )
            numbers.push_back(input) ;
    }while ( std::cin && std::cin.peek() != '\n' );

    return numbers ;
}

// utility for topological sorting of alternative graph
void topologicalSortUtil(int v, vector<bool> &visited,  stack<int> &Stack, vector< vector<int> > &adj)
{
    visited[v] = true;
 
    vector<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            topologicalSortUtil(*i, visited, Stack, adj);
 
    //cout<<v;
    Stack.push(v);
}

double payoffc (vector<int> s, alternative nodes[] )
{
	int TotalCost = 0;
	double payoff = 0;
	for(int i = 0 ; i < s.size(); i++) {
		if (s[i] != -1){
		TotalCost = nodes[i].cost + TotalCost;
		}
	}
	//cout<<"totalcost"<<TotalCost<<endl;
	payoff = TotalCost * PTC;
	return payoff;

}

vector<int> NonUniqueFinish(vector<int> s, alternative nodes[])
{
	vector<int> nufv;
	int f = 0, j = 0, count = 0, k = 0, flag = 0;
	for (int i = 0; i<s.size(); i++)
	{
		if(s[i] != -1)
		{
			f = s[i] + nodes[i].duration;
			count = 0;
			for (j = 0; j<s.size(); j++)
			{
				if (s[j] != -1)
				{
					if (f == (s[j] + nodes[j].duration))
					{
						count++; 
					}
				
				}
			}
			if (count > 1)
			{   
				flag = 0;
				count = 0;
				for (k = 0; k<nufv.size(); k++)
				{
					if (f == nufv[k])
					{flag = 1;}
				}
				if (flag == 0)
				{nufv.push_back(f);}
			}
		}
	}
	
	return nufv;
	
}

vector<int> FinishBeforeAt(vector<int> s, int time, alternative nodes[])
{
	vector<int> fbav;
	int f;
	for (int i = 0; i<s.size(); i++)
	{
		if (s[i] != -1)
		{
			f = s[i] + nodes[i].duration;
			if (f <= time)
			{
				fbav.push_back(i);
			}
		}
	}
	
	return fbav;
	
}

vector<int> FinishBefore(vector<int> s, int time, alternative nodes[])
{
	vector<int> fbv;
	int f;
	for (int i = 0; i<s.size(); i++)
	{
		if (s[i] != -1)
		{
			f = s[i] + nodes[i].duration;
			if (f < time)
			{
				fbv.push_back(i);
			}
		}
	}
	
	return fbv;
}

vector<int> FinishAt(vector<int> s, int time, alternative nodes[])
{
	vector<int> fav;
	int f = 0;
	for (int i = 0; i<s.size(); i++)
	{
		if (s[i] != -1)
		{
			f = s[i] + nodes[i].duration;
			if (f == time)
			{
				fav.push_back(i);
			}
		}
	}
	
	return fav;
}

vector<int> UniqueFinish(vector<int> s, alternative nodes[])
{
	vector<int> uf;
	int f = 0, j = 0, count = 0;
	for (int i = 0; i<s.size(); i++)
	{
		if(s[i] != -1)
		{
			f = s[i] + nodes[i].duration;
			count = 0;
			//uf.push_back(i);
			for (j = 0; j<s.size(); j++)
			{
				if (s[j] != -1)
				{
					if (f == (s[j] + nodes[j].duration))
					{
						count++; 
					}
				
				}
			}
			if (count < 2)
			{   
				uf.push_back(i);
			}
		}
	}
	
	return uf;
} 

double enpv(vector<int> s, alternative nodes[])
{
	double npv = 0;
	double prob = 1, first1 = 0, first = 0;
	double payoff = payoffc(s,nodes);
	int j = 0, i = 0;
	vector<int> fba;
	for (i = 0; i<s.size(); i++) //for first term
	{
		if (s[i] != -1)
		{
			first1 = exp(-R*s[i]) * nodes[i].cost;	
			//cout<<first1;
			fba = FinishBeforeAt(s, s[i], nodes);
			prob = 1;
			for (j = 0; j<fba.size(); j++)
			{
				//cout<<fba[j];
				prob = prob * (1-nodes[fba[j]].p);
			}
			first1 = first1 * prob;
			first = first + first1;
			//cout<<first;
		}
	}
	//cout<<"first"<<first<<endl;
	
	double prob1 = 1, prob2 = 1, second=0, second1=0, second2=0, second3=0;
	vector<int> nuf = NonUniqueFinish(s, nodes);
	int k = 0;
	vector<int> fa, fb;
	for (i = 0; i<nuf.size(); i++) //for second term | t = nuf(i)
	{	
		second1 = exp(-R*nuf[i]);
		fa = FinishAt(s, nuf[i], nodes);
		prob1 = 1;
		for (k = 0; k<fa.size(); k++)
		{
			prob1 = prob1 * (1-nodes[fa[k]].p);
		}
		second2 = 1 - prob1;
		fb = FinishBefore(s, nuf[i], nodes);
		prob2 = 1;
		for (j = 0; j<fb.size(); j++)
		{
			prob2 = prob2 * (1-nodes[fb[j]].p);
		}
		second3 = prob2;
		second  = second + (second1 * second2 * second3);
	}
	second = second * payoff;
	//cout<<"second"<<second<<endl;
	
	double third = 0, third1 = 0, third2 = 0; prob = 1;  
	vector<int> ufa = UniqueFinish(s, nodes); 
	for(i = 0; i<ufa.size(); i++) //for third term
	{
		third1 = exp(-R*(s[ufa[i]] + nodes[ufa[i]].duration)) * nodes[ufa[i]].p;
		fb = FinishBefore(s, (s[ufa[i]] + nodes[ufa[i]].duration), nodes);
		prob = 1;
		for (j = 0; j<fb.size(); j++)
		{
			prob = prob * (1-nodes[fb[j]].p);
		}
		third2 = prob;
		third = third + (third1 * third2);
	}
	third = third * payoff;
	//cout<<"third"<<third<<endl;
	
	npv = first + second + third;
	//cout<<"final"<<npv<<endl;
	return npv;
}

void topologicalSort(int n, vector<bool> &visitem, vector<int> &update, vector< vector<int> > &adj)
{
    visitem[n] = true;
 
    vector<int>::iterator i;
    for (i = adj[n].begin(); i != adj[n].end(); ++i)
        if (!visitem[*i])
            topologicalSort(*i, visitem, update, adj);
 
    //cout<<n;
    update.push_back(n);
}

vector<int> SA(vector<int> s)
{
	int i;
	vector<int> sav;
	for(i=0; i<s.size();i++)
	{
		if(s[i] != -1)
		{sav.push_back(i);}
	}
	return sav;
}

vector<int> RA(vector<int> s)
{
	int i;
	vector<int> rav;
	for(i=0; i<s.size();i++)
	{
		if(s[i] == -1)
		{rav.push_back(i);}
	}
	return rav;
}

vector<int> intersect(vector<int> v1, vector<int> v2)
{
	vector<int> intersection;
	int i, j;
	for (i = 0; i < v1.size(); i++)
	{
		for (j=0; j < v2.size(); j++)
		{
			if(v1[i] == v2[j])
			{
				intersection.push_back(v1[i]);
			}
		}
	}
	return intersection;
}

double UB(vector<int> s, alternative nodes[], vector< vector<int> > &adj, vector< vector<int> > &pred)
{
		// updating earliest start and finish times for each task

	int top, i, j, k, t;
	int n_tasks = s.size();
	vector<int> update;
	vector<bool> visits(n_tasks, false);
	topologicalSort(0, visits, update, adj);
	nodes[0].es = 0;
	nodes[0].ef = 0;
	update.pop_back();

	while(!update.empty()) {
		top = update.back();
		//cout<<top;
		int max_f = -1;
		for(i = 0; i < pred[top].size(); i++) {
			if(max_f < nodes[pred[top][i]].ef) {
				max_f = nodes[pred[top][i]].ef;
			} 
		}
		if (s[top] != -1)
		{
			nodes[top].es = s[top];
			nodes[top].ef = s[top] + nodes[top].duration;
		}
		else
		{
		nodes[top].es = max_f;
		nodes[top].ef = max_f + nodes[top].duration;
		}
		update.pop_back();
	}

	/*  //FOR DEBUG
		cout<<"Es and Ef : \n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" "<<i+1<<" "<<nodes[i].es<<" "<<nodes[i].ef<<endl;
		}
	*/
	
	// updating latest start and finish time for each task

	vector<int> update2;
	vector<bool> visits2(n_tasks, false);

	topologicalSort(n_tasks-1, visits2, update2, pred);

	nodes[n_tasks-1].ls = nodes[n_tasks-1].es;
	nodes[n_tasks-1].lf = nodes[n_tasks-1].ef;
	update2.pop_back();
	while(!update2.empty()) {
		top = update2.back();
		int min_s = 99999;
		for(i = 0; i < adj[top].size(); i++) {
			if(min_s > nodes[adj[top][i]].ls) {
				min_s = nodes[adj[top][i]].ls;
			} 
		}
		nodes[top].lf = min_s;
		nodes[top].ls = min_s - nodes[top].duration;
		update2.pop_back();
	}
	/*	// FOR DEBUG
		cout<<"Ls and Lf : \n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" "<<i+1<<" "<<nodes[i].ls<<" "<<nodes[i].lf<<endl;
		}
	*/

	//--------start UB-----------//
	double bound;
	vector<int> salt = SA(s);
	vector<int> sav;
	vector<int> fba;
	vector<int> sa_fba;
	vector<int> rav;
	double first1 = 0, first2 = 0, first3 = 0, first = 0, prob1 = 1, prob2 = 1;
	for(i=0;i<salt.size();i++)
	{
		first1 = exp(-R*s[salt[i]]) * nodes[salt[i]].cost;
		sav = SA(s);
		fba = FinishBeforeAt(s, s[salt[i]], nodes);
		sa_fba = intersect(sav,fba);
		prob1 = 1;
		for(j=0; j<sa_fba.size(); j++)
		{
			prob1 = prob1 * (1-nodes[sa_fba[j]].p);
		}
		first2 = prob1;
		rav = RA(s);
		prob2 = 1;
		for(k=0; k<rav.size(); k++)
		{
			if(nodes[rav[k]].ef <= s[salt[i]])
			{
				prob2 = prob2 * (1-nodes[rav[k]].p);
			}
		}
		first3 = prob2;
		first = first + (first1 * first2 * first3);
	}
	
	vector<int> ralt = RA(s);
	double second = 0, second1 = 0, second2 = 0, second3 = 0;
	for(i=0;i<ralt.size();i++)
	{
		second1 = exp(-R*nodes[ralt[i]].ls) * nodes[ralt[i]].cost;
		sav = SA(s);
		prob1 = 1;
		for(j=0; j<sav.size(); j++)
		{
			if((s[sav[j]] + nodes[sav[j]].duration) <= nodes[ralt[i]].ls)
			{
				prob1 = prob1 * (1-nodes[sav[j]].p);
			}
		}
		second2 = prob1;
		rav = RA(s);
		prob2 = 1;
		for(k=0; k<rav.size(); k++)
		{
			if(nodes[rav[k]].ef <= nodes[ralt[i]].ls)
			{
				prob2 = prob2 * (1-nodes[rav[k]].p);
			}
		}
		second3 = prob2;
		second = second + (second1 * second2 * second3);
	}
	
	double payoff = payoffc(s,nodes);
	int prob3 = 1;
	vector<int> fa;
	vector<int> sa_fa;
	vector<int> sa_fb;
	vector<int> fb;
	double third = 0, third1 = 0, third2 = 0, third3 = 0, third4 = 0;
	vector<int> nuf = NonUniqueFinish(s, nodes);
	for(t = 0; t<nuf.size(); t++)
	{
		third1 = exp(-R*nuf[t]);
		sav = SA(s);
		fa = FinishAt(s, nuf[t], nodes);
		sa_fa = intersect(sav,fa);
		prob1 = 1;
		for(k=0; k<sa_fa.size(); k++)
		{
			prob1 = prob1 * (1 - nodes[sa_fa[k]].p);
		}
		third2 = 1 - prob1;
		
		sav = SA(s);
		fb = FinishBefore(s, nuf[t], nodes);
		sa_fb = intersect(sav,fb);
		prob2 = 1;
		for(j=0; j<sa_fb.size(); j++)
		{
			prob2 = prob2 * (1 - nodes[sa_fb[j]].p);
		}
		third3 = prob2;
		
		rav = RA(s);
		int l;
		prob3 = 1;
		for(l = 0; l<rav.size(); l++)
		{
			if(nodes[rav[l]].lf < nuf[t])
			{
				prob3 = prob3 * (1 - nodes[rav[l]].p);
			}
		}
		third4 = prob3;
		third = third + (third1 * third2 * third3 * third4);
	}
	third = third * payoff;
	
	double fourth = 0, fourth1 = 0, fourth2 = 0, fourth3 = 0;
	int finish, f, flag = 0;
	salt = SA(s);
	nuf = NonUniqueFinish(s, nodes);
	for(i=0; i<salt.size(); i++)
	{
		finish = s[salt[i]] + nodes[salt[i]].duration;
		flag = 0;
		for (f = 0; f<nuf.size(); f++)
		{
			if(finish == nuf[f])
			{ flag = 1; }
		}
		if(flag == 0)
		{
			fourth1 = exp(-R*finish) * nodes[salt[i]].p;
			sav = SA(s);
			fb = FinishBefore(s, finish, nodes);
			sa_fb = intersect(sav,fb);
			prob1 = 1;
			for(j=0; j<sa_fb.size(); j++)
			{
				prob1 = prob1 * (1 - nodes[sa_fb[j]].p);
			}
			fourth2 = prob1;
			
			rav = RA(s);
			prob2=1;
			for(k=0; k<rav.size(); k++)
			{
				if(nodes[rav[k]].lf < finish)
				{
					prob2 = prob2 * (1 - nodes[rav[k]].p);
				}
			}
			fourth3 = prob2;
			fourth = fourth + (fourth1 * fourth2 * fourth3);
		}
	}
	fourth = fourth * payoff;
	double fifth=0, fifth1 = 0, fifth2 = 0, fifth3 = 0;
	ralt = RA(s);
	for(i=0; i<ralt.size(); i++)
	{
		fifth1 = exp(-R*(nodes[ralt[i]].es + nodes[ralt[i]].duration)) * nodes[ralt[i]].p;
		sav = SA(s);
		prob1 = 1;
		for(j=0; j<sav.size(); j++)
		{
			if ((s[sav[j]] + nodes[sav[j]].duration) < nodes[ralt[i]].ef)
			{
				prob1 = prob1 * (1 - nodes[sav[j]].p);
			}
		}
		fifth2 = prob1;
		rav = RA(s);
		prob2 = 1;
		for(k=0; k<rav.size(); k++)
		{
			if(nodes[rav[k]].lf < nodes[ralt[i]].ef)
			{
				prob2 = prob2 * (1 - nodes[rav[k]].p);
			}
		}
		fifth3 = prob2;
		fifth = fifth + (fifth1 * fifth2 * fifth3);
	}
	fifth = fifth * payoff;
	bound = first + second + third + fourth + fifth;
	return bound;
}

vector<int> EL(vector<int> s, vector< vector<int> > &pred)
{
	vector<int> sav = SA(s);
	vector<int> rav = RA(s);
	vector<int> eligible;
	int count = 0;
	int check;
	for(int i = 0; i<rav.size(); i++)
	{
		count = 0;
		for(int j = 0 ; j < pred[rav[i]].size(); j++) 
		{
			check = pred[rav[i]][j];
			for(int k = 0; k<sav.size(); k++)
			{
				if(check == sav[k])
				{
					count++;
					break;
				}
			}	
		}
		if(count == pred[rav[i]].size())
		{ eligible.push_back(rav[i]); }
	}
	sort(eligible.begin(), eligible.end());
	return eligible;
}

double EST(int elem, vector<int> s,  vector< vector<int> > &pred, alternative nodes[]) //max finish time amongst predecessors
{
	double start = 0;
	for(int i = 0; i<pred[elem].size(); i++)
	{
		if (start < (s[pred[elem][i]] + nodes[pred[elem][i]].duration))
		{
			start = s[pred[elem][i]] + nodes[pred[elem][i]].duration;
		}
	}
	return start;
}

vector<int> Possible_Start_Times(int elem, vector<int> s, int est, int eft, alternative nodes[])
{
	vector<int> pst;
	pst.push_back(est);
	vector<int> sav = SA(s);
	int flag = 0;
	for(int i = 0; i<sav.size(); i++)
	{
		if(s[sav[i]] > est)
		{
			flag = 0;
			for(int j=0; j<pst.size(); j++)
			{
				if(s[sav[i]] == pst[j])
				{flag = 1;}
			}
			if(flag == 0)
			{pst.push_back(s[sav[i]]);}
		}
	}
	
	for(int i = 0; i<sav.size(); i++)
	{
		if((s[sav[i]] + nodes[sav[i]].duration) > est)
		{
			flag = 0;
			for(int j=0; j<pst.size(); j++)
			{
				if((s[sav[i]] + nodes[sav[i]].duration) == pst[j])
				{flag = 1;}
			}
			if(flag == 0)
			{pst.push_back(s[sav[i]] + nodes[sav[i]].duration);}
		}
	}
	
	for(int i = 0; i<sav.size(); i++)
	{
		if(s[sav[i]] > eft)
		{
			flag = 0;
			for(int j=0; j<pst.size(); j++)
			{
				if((s[sav[i]] - nodes[elem].duration) == pst[j])
				{flag = 1;}
			}
			if(flag == 0)
			{pst.push_back(s[sav[i]] - nodes[elem].duration);}
		}
	}

	for(int i = 0; i<sav.size(); i++)
	{
		if((s[sav[i]] + nodes[sav[i]].duration) > eft)
		{
			flag = 0;
			for(int j=0; j<pst.size(); j++)
			{
				if(((s[sav[i]] + nodes[sav[i]].duration) - nodes[elem].duration) == pst[j])
				{flag = 1;}
			}
			if(flag == 0)
			{pst.push_back((s[sav[i]] + nodes[sav[i]].duration) - nodes[elem].duration);}
		}
	}
	
	sort(pst.begin(), pst.end());
	return pst;
}

int main() { 

	
	int i,n_tasks,top,j;

	cout<<"############## Exact Method for Alternative Scheduling ################\n\n";
	cout<<"Enter the number of tasks : ";
	cin>>n_tasks; // n_tasks is the number of tasks
	


	struct alternative nodes[n_tasks]; // number of activities here 0th alternative is the start
									  // and the (n+1)th alternative refers finish 
									  // both having duration 0

	nodes[0].duration = 0;
	nodes[0].cost = 0;
	nodes[0].p = 1.000000;
	nodes[n_tasks-1].duration = 0;
	nodes[n_tasks-1].cost = 0;
	nodes[n_tasks-1].p = 0.000000;
	 
/*
	// FOR DEBUG
	int dur[] = {5,7,4,3,5,10};
	int cos[] = { 24  ,    69  ,    50  ,    41  ,    15  ,    64     };
	double pr[] = {      0.583117   ,   0.831523  ,    0.725394  ,    0.676061    ,  0.528520    ,  0.803842 };
	for(i = 1 ; i < n_tasks-1; i++) {
		nodes[i].cost = cos[i-1];
		nodes[i].duration = dur[i-1];
		nodes[i].p = pr[i-1];
	}
*/		
	// FOR INPUT
int dur[] = {7,2,8,3,6,9,1,6,1,3,7,5};
int cos[] = { 30    ,  87 ,     30   ,   32    ,  58   ,   98  ,    14   ,   17    ,  57 ,     48    ,  18  ,    33  };
double pr[] = {  0.614246   ,   0.931791     , 0.614719    ,  0.624775     , 0.771203   ,   0.992416  ,    0.526902 ,     0.540712    ,  0.762337 ,     0.713401 ,     0.547334    ,  0.629398  };
	
	for(i = 1 ; i < n_tasks-1; i++) {
		nodes[i].cost = cos[i-1];
		nodes[i].duration = dur[i-1];
		nodes[i].p = pr[i-1];
	}

/*	// input of all the tasks manually
	for(i = 1 ; i <= n_tasks-2; i++) {
		cout<<"Enter duration for "<<i+1<<" : ";
		cin>>nodes[i].duration;
	}

	for(i = 1 ; i <= n_tasks-2; i++) {
		cout<<"Enter cost for "<<i+1<<" : ";
		cin>>nodes[i].cost;
	}
	
		for(i = 1 ; i <= n_tasks-2; i++) {
		cout<<"Enter PTS for "<<i+1<<" : ";
		cin>>nodes[i].p;
	}
*/

	cout<<"\n\n\t\tTasks and Durations :\n";
	for(i = 0 ; i <= n_tasks-1; i++) {
		cout<<"\t\t"<<i+1<<". "<<" "<<nodes[i].duration<<endl;
	}
	
	vector< vector<int> > adj;  // adj represents sucessor list
	vector< vector<int> > pred; // pred reperesents predecessor list


	// initialization of both lists with empty vectors
	for(i = 0 ; i <= n_tasks-2; i++) {
		vector<int> temp;
		adj.push_back(temp);
		pred.push_back(temp);
	}

	// initialization of successor list based on terminal input
	// NOTE : input all the tasks with no predecessors as the successor of 0
	cout<<"\n\nNOTE : Enter all the tasks with no predecessors as the successor of 0";
	for(i = 0 ; i <= n_tasks-2; i++) {
		cout<<"\n\nEnter successors for task "<<i+1<<" : ";
		vector<int> temp = ReadNumbers();
		if(temp.size()==0){
			adj[i].push_back(n_tasks-2);
			pred[n_tasks-2].push_back(i);

		}
		else {
			for(int j=0; j<temp.size(); j++)
				adj[i].push_back(temp[j]-1);
			for(int j=0;j < temp.size(); j++)
				pred[temp[j]-1].push_back(i);

		}

	}
	if(DBG) {
		//debugging
		cout<<"\nSuccessor matrix :\n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" : ";
			for(j = 0 ; j < adj[i].size(); j++) {
				cout<<adj[i][j]+1<<", ";
			}
			cout<<endl;
		}

		cout<<"Predecessor matrix :\n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" : ";
			for(j = 0 ; j < pred[i].size(); j++) {
				cout<<pred[i][j]+1<<", ";
			}
			cout<<endl;
		}
	}

	// ----------- Begin CPM for ef,lf,es,ls-----------------------------//

	// calculating earliest start and finish times for each task
	// topological sort of task is required here

	stack<int> Stack;
	vector<bool> visit(n_tasks, false);
	topologicalSortUtil(0,visit, Stack, adj);
	
	nodes[0].es = 0;
	nodes[0].ef = 0;
	Stack.pop();

	while(!Stack.empty()) {
		top = Stack.top();
		//cout<<top;
		int max_f = -1;
		for(i = 0; i < pred[top].size(); i++) {
			if(max_f < nodes[pred[top][i]].ef) {
				max_f = nodes[pred[top][i]].ef;
			} 
		}
		nodes[top].es = max_f;
		nodes[top].ef = max_f + nodes[top].duration;
		Stack.pop();
	}

/*
	if(DBG) {
		cout<<"Es and Ef : \n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" "<<i+1<<" "<<nodes[i].es<<" "<<nodes[i].ef<<endl;
		}
	}
*/
	// calculating latest start and finish time for each task

	stack<int> Stack2;
	vector<bool> visit2(n_tasks, false);

	topologicalSortUtil(n_tasks-1, visit2, Stack2, pred);

	nodes[n_tasks-1].ls = nodes[n_tasks-1].es;
	nodes[n_tasks-1].lf = nodes[n_tasks-1].ef;
	Stack2.pop();
	while(!Stack2.empty()) {
		top = Stack2.top();
		int min_s = 99999;
		for(i = 0; i < adj[top].size(); i++) {
			if(min_s > nodes[adj[top][i]].ls) {
				min_s = nodes[adj[top][i]].ls;
			} 
		}
		nodes[top].lf = min_s;
		nodes[top].ls = min_s - nodes[top].duration;
		Stack2.pop();
	}
	//cout<<"\n\n";

	 //queue<int> q2;
	 //vector<int> visited2(n_tasks,0);
	 //q2.push(n_tasks-2);

/*
	if(DBG) {
		cout<<"Ls and Lf : \n";
		for(i = 0 ; i < n_tasks; i++) {
			cout<<i+1<<" "<<i+1<<" "<<nodes[i].ls<<" "<<nodes[i].lf<<endl;
		}
	}
*/

// ----------- END CPM for ef,lf,es,ls-----------------------------//
	
	// UB DEBUG
	// ENPV DEBUG
	
	int pin;
	cout<<endl<<"Schedule entry (1) or schedule compute (0) ?: ";
	cin>>pin;
	cout<<endl;
	if(pin)
	{
		vector<int> test;
		int test_val;
		for (i = 0; i<n_tasks; i++)
		{
			cout<<"enter schedule index "<<i+1<<": ";
			cin>>test_val;
			test.push_back(test_val);
		}
		
		double test_enpv = enpv(test, nodes);
		cout<<endl<<"ENPV(TEST): "<<test_enpv;

		double test2 = UB(test,nodes,adj,pred);
		cout<<endl<<"UB: "<<test2<<endl;
	}

//--------------start Branch and Bound (F-B&B)---------------------//

	if(!pin)
	{
		vector< vector<int> > frontier;
		vector< vector<int> > explored;
		vector<int> front;
		int flag = 0;
		// ---Initialize frontier vector and LB---//
		vector<int> tmp;
		//vector<int> tmp2;
		tmp.push_back(0);
		for(i = 1; i<n_tasks; i++)
			{tmp.push_back(-1);}
		frontier.push_back(tmp);
		double lb = 0;
		double bound;
		int goal_flag = 0;
		int c = 0;
		double upper = 0;
		
		// ----------algo start------------//
		while(!frontier.empty())
		{
			c++; //debug
			front = frontier.back();
			frontier.pop_back();
			flag = 0;
			//check if front is in explored
			for(i = 0; i<explored.size(); i++)
			{
				if (front == explored[i])
					{flag = 1; break;}
			}
			if (flag == 1)
			{
				continue;
			}
			//for (int x=0; x<front.size(); x++)
			//cout<<front[x];
			bound = UB(front, nodes, adj, pred);
			//cout<< bound;
			if (bound < lb)
			{
				explored.push_back(front);
				continue;
			}
			
			//lb = bound;
			if(bound > upper)
			upper = bound;
			double lb_local = enpv(front, nodes);
			if (lb_local > lb)
			lb = lb_local;
			 
			explored.push_back(front);
			vector<int> goal = RA(front);
			if(goal.empty())
			{
				goal_flag = 1;
				break;
			}
			vector<int> eligible = EL(front, pred);
			vector<int> pst;
			vector<int> vec;
			
			while(!eligible.empty())
			{
				vec = front;
				int elem = eligible.back();
				eligible.pop_back();
				int est = EST(elem, front, pred, nodes);
				int eft = est + nodes[elem].duration;
				pst = Possible_Start_Times(elem, front, est, eft, nodes);
				while(!pst.empty())
				{
					int add = pst.back();
					pst.pop_back();
					vec[elem] = add;
					frontier.push_back(vec);
				}
			}
		}
		if(goal_flag == 1)
		{
			cout<<endl<<"Solution found: ";
			for (int n=0; n<front.size()-1; n++)
			cout<<front[n]<<" , ";
			cout<<front[front.size()-1];
			double fub = UB(front, nodes, adj, pred);
			double fnpv = enpv(front, nodes);
			cout<<endl<<"NPV: "<<fnpv;
			cout<<endl<<"UB: "<<fub;
			//cout<<endl<<"upper: "<<upper;
		}
		else{
			cout<<"error!";
			cout<<c;
		}
	}
	return 0;

}
