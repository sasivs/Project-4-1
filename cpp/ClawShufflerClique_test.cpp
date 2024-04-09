#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

using namespace std;

string EdgeFile;
int NodeNum;
double Eps, Eps_l, Eps_1, Eps_2, c;
string Eps_s, Eps_l_s, Eps_1_s, Eps_2_s, c_s, t_s, delta_s;

int ItrNum;
string Alg;
int t;

// Initialization of statslib
stats::rand_engine_t engine(1776);

//opens file and returns a file pointer to it
FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

int compare_double(double *x, double *y) {
	if (*x > *y)       return(1);   /* return positive integer */
	else if (*x == *y) return(0);   /* return zero     integer */
	else              return(-1);  /* return negative integer */
}

bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num) {
	int rnd;
	int *ordperm;
	int i, j;

	// 0, 1, 2, ..., size-1 --> ordperm
	ordperm = (int *)malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		ordperm[i] = i;
	}

	for (i = 0; i < num; i++) {
		rnd = genrand_int32() % (size - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < size - i; j++) {
			ordperm[j - 1] = ordperm[j];
		}
	}

	free(ordperm);
}

// Read edges from the edge file
void ReadEdges(map<int, int> *a_mat, int *node_order){
	int node1, node2;
	int i;
	char s[1025];
	char *tok;
	FILE *fp;

	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<3;i++) fgets(s, 1024, fp);
	while(fgets(s, 1024, fp) != NULL){
		// 1st node --> node1
		tok = strtok(s, ",");
		node1 = atoi(tok);
		// 2nd node --> node2
		tok = strtok(NULL, ",");
		node2 = atoi(tok);
		if(node1 == node2) continue;
		// If both nodes exist, add the edge
		if(node_order[node1] < NodeNum && node_order[node2] < NodeNum){
			a_mat[node_order[node1]][node_order[node2]] = 1;
			a_mat[node_order[node2]][node_order[node1]] = 1;
		}
	}
	fclose(fp);
}

void swap(int* randperm, int i, int j){
	int temp;
	temp = randperm[i];
	randperm[i] = randperm[j];
	randperm[j] = temp;
}


void CalIFCliqCountingShuffler(map<int, int> *a_mat, string outfile, double &cliq_4_num1, double &cliq_4_num2, double &cliq_4_num3, \
double &clip_cliq_4_num1, double &clip_cliq_4_num2, double &clip_cliq_4_num3){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<pair<int,int>, int> *wedge_ns;  // noisy wedge matrix
	int i,j,k,itr,itr1;
	int zij, zjk, zki, zji, zkj, zik, partial_sum_emp1, partial_sum_emp2, partial_sum_emp3;
    int  wi, wj, wk;
	double ql, q, q_1_3, ql_1, q_1_2;
	double rnd;
	double star_3;
	int* randperm;
	double cliq_4_num_emp1, cliq_4_num_clip1, cliq_4_num_emp2, cliq_4_num_clip2, cliq_4_num_emp3, cliq_4_num_clip3;
	int *deg_ns;
	double deg_avg, deg_th;

	map<int, int>:: iterator aitr;

	a_mat_ns = new map<int, int>[NodeNum];

	wedge_ns = new map<pair<int,int>, int>[NodeNum];

	malloc1D(&randperm, NodeNum);

	malloc1D(&deg_ns, NodeNum);
	
	MakeRndPerm(randperm, NodeNum, NodeNum);

	ql = 1 / (exp(Eps_l) + 1);

	q = 1 / (exp(Eps_2) + 1);

	q_1_3 = 1 / pow(1-2*q, 3);
	
	q_1_2 = 1 / pow(1-2*q, 2);

	ql_1 = 1 / (1-2*ql);

	cout<<ql<<" "<<q<<" "<<q_1_3<<" "<<q_1_2<<" "<<ql_1<<endl;

	deg_avg = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg_ns[i] += 1;
		deg_ns[i] += stats::rlaplace(0.0, 1/Eps_1, engine);
		deg_avg += deg_ns[i];
	}
	deg_avg = deg_avg / (double) NodeNum;

	deg_th = c * deg_avg;

	cliq_4_num_emp1 = 0; cliq_4_num_clip1 = 0; cliq_4_num_emp2 = 0; cliq_4_num_clip2 = 0; cliq_4_num_emp3 = 0; cliq_4_num_clip3 = 0;
	for(itr=0; itr<t; itr+=1){
		i = randperm[itr];
		j = randperm[itr+1];
		k = randperm[itr+2];
		star_3 = 0;
		for(itr1=0; itr1<NodeNum; itr1++){
			if(itr1 == i || itr1 == j || itr1 == k) continue;
			rnd = genrand_real2();
			star_3 -= ql;
			if(a_mat[itr1].count(i) == 1 && a_mat[itr1].count(j) == 1 && a_mat[itr1].count(k) == 1){
				if(rnd >= ql) star_3++;
			}
			else if (rnd < ql) star_3++; 
		}

		if (a_mat_ns[i].count(j) > 0) {
			zij = a_mat_ns[i][j];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[i].count(j) == 0) || (rnd >= q && a_mat[i].count(j) == 1)) zij = 1;
			else zij = 0;
		}
		if (a_mat_ns[j].count(k) > 0) {
			zjk = a_mat_ns[j][k];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[j].count(k) == 0) || (rnd >= q && a_mat[j].count(k) == 1)) zjk = 1;
			else zjk = 0;
		}
		if (a_mat_ns[k].count(i) > 0){
			zki = a_mat_ns[k][i];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[k].count(i) == 0) || (rnd >= q && a_mat[k].count(i) == 1)) zki = 1;
			else zki = 0;
		}
		if (a_mat_ns[j].count(i) > 0){
			zji = a_mat_ns[j][i];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[j].count(i) == 0) || (rnd >= q && a_mat[j].count(i) == 1)) zji = 1;
			else zji = 0;
		}
		if (a_mat_ns[k].count(j) > 0){
			zkj = a_mat_ns[k][j];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[k].count(j) == 0) || (rnd >= q && a_mat[k].count(j) == 1)) zkj = 1;
			else zkj = 0;
		}
		if (a_mat_ns[i].count(k) > 0){
			zik = a_mat_ns[i][k];
		}
		else{
			rnd = genrand_real2();
			if ((rnd < q && a_mat[i].count(k) == 0) || (rnd >= q && a_mat[i].count(k) == 1)) zik = 1;
			else zik = 0;
		}

		partial_sum_emp1 = (zij-q) * (zjk-q) * (zki-q) + (zij-q) * (zjk-q) * (zik-q) + 
					(zij-q) * (zkj-q) * (zki-q) + (zij-q) * (zkj-q) * (zik-q) + 
					(zji-q) * (zjk-q) * (zki-q) + (zji-q) * (zjk-q) * (zik-q) + 
					(zji-q) * (zkj-q) * (zki-q) + (zji-q) * (zkj-q) * (zik-q);

		pair<int, int> key_pair_i = {j,k};
		if (wedge_ns[i].count(key_pair_i) > 0){
			wi = wedge_ns[i][key_pair_i];
		}
		else{
			rnd = genrand_real2();
			if (a_mat[i].count(j) > 0 && a_mat[i].count(k) > 0){
				if (rnd >= q) wi=1;
			} 
			else if (rnd < q) wi=1;
			else wi=0;
		}

		pair<int,int> key_pair_j = {i,k};
		if (wedge_ns[j].count(key_pair_j) > 0){
			wj = wedge_ns[j][key_pair_j];
		}
		else{
			rnd = genrand_real2();
			if (a_mat[j].count(i) > 0 && a_mat[j].count(k) > 0){
				if (rnd >= q) wj=1;
			} 
			else if (rnd < q) wj=1;
			else wj=0;
		}

		pair<int,int> key_pair_k = {i,j};
		if (wedge_ns[k].count(key_pair_k) > 0){
			wk = wedge_ns[k][key_pair_k];
		}
		else{
			rnd = genrand_real2();
			if (a_mat[k].count(i) > 0 && a_mat[k].count(j) > 0){
				if (rnd >= q) wk=1;
			} 
			else if (rnd < q) wk=1;
			else wk=0;
		}
		
		partial_sum_emp2 = (wi-q) * (zjk-q) + (wi-q) * (zkj-q) +
					(wj-q) * (zik-q) + (wj-q) * (zki-q) +
					(wk-q) * (zij-q) + (wk-q) * (zji-q);

		partial_sum_emp3 = (wi-q) * (wj-q) + (wi-q) * (wk-q) + (wj-q) * (wk-q);

		cliq_4_num_emp1 += partial_sum_emp1 * star_3 * ql_1 * q_1_3 / 8.0;

		cliq_4_num_emp2 += partial_sum_emp2 * star_3 * ql_1 * q_1_2 / 6.0;

		cliq_4_num_emp3 += partial_sum_emp3 * star_3 * ql_1 * q_1_2 / 3.0;

		if(min(min(deg_ns[i], deg_ns[j]), deg_ns[k]) > deg_th){
			cliq_4_num_clip1 += partial_sum_emp1 * star_3 * ql_1 * q_1_3 / 8.0;
			cliq_4_num_clip2 += partial_sum_emp2 * star_3 * ql_1 * q_1_2 / 6.0;
			cliq_4_num_clip3 += partial_sum_emp3 * star_3 * ql_1 * q_1_3 / 3.0;
		}
	} 

	cliq_4_num_emp1 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_emp1 / (24*t);
	cliq_4_num_emp2 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_emp2 / (24*t);
	cliq_4_num_emp3 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_emp3 / (24*t);

	cliq_4_num_clip1 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_clip1 / (24*t);
	cliq_4_num_clip2 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_clip2 / (24*t);
	cliq_4_num_clip3 = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_clip3 / (24*t);

	cliq_4_num1 = cliq_4_num_emp1;
	cliq_4_num2 = cliq_4_num_emp2;
	cliq_4_num3 = cliq_4_num_emp3;

	clip_cliq_4_num1 = cliq_4_num_clip1;
	clip_cliq_4_num2 = cliq_4_num_clip2;
	clip_cliq_4_num3 = cliq_4_num_clip3;

	free1D(deg_ns);
	free1D(randperm);
	delete[] a_mat_ns;
	delete[] wedge_ns;
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int **node_order;
	int *deg;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<int, int>::iterator aitr3;

	long long cliq_num;

	double cliq_4_num1, cliq_4_num2, cliq_4_num3, clip_cliq_num1, clip_cliq_num2, clip_cliq_num3, cliq_4_num_re1, cliq_4_num_re2, cliq_4_num_re3, 
	clip_cliq_num_re1, clip_cliq_num_re2, clip_cliq_num_re3, cliq_4_num_avg_re1, cliq_4_num_avg_re2, cliq_4_num_avg_re3, 
	clip_cliq_num_avg_re1, clip_cliq_num_avg_re2, clip_cliq_num_avg_re3; 
	
	int itr;
	int i, j, k, x, l;
	string outdir;
	string outfile;
	char s[1025], *str;
	char str_1[] = "1";
	char *tok;
	FILE *fp;

	int fix_perm;

	// Initialization of Mersennne Twister
	unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
	init_by_array(init, length);
	//	init_genrand(1);

	if (argc < 2) {
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon] [epsilon-1] [epsilon-2] [epsilon-l] [delta] [c] [t (default: Nodenum/2)] [#itr(-1) (default: 1)] \n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon]: Parameters epsilon whole privacy budget and (ep=ep1+ep2)\n");
		printf("[epsilon-1]: Parameters epsilon-1 (Lap Noise Parameter)\n");
		printf("[epsilon-2]: Parameters epsilon-2 (Wedge shuffle alg parameter to calculate ep-l)\n");
		printf("[epsilon-l]: Parameters epsilon-l (Wedge shuffle amplification parameter)\n");
		printf("[delta]: Parameters delta (Wedge shuffle alg parameter to calculate ep-l)\n");
		printf("[c]: Parameters c (To calculate deg-threshold, bias and accuracy tradeoff)\n");
		printf("[t]: Parameters t (To generate random user pairs, time and accuracy tradeoff)\n");
		printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
		return -1;
	}

	EdgeFile = argv[1]; //name of the edge file

	NodeNum = -1;
	if (argc >= 3) NodeNum = atoi(argv[2]); //#of nodes

	cout<<"#Nodes: "<<NodeNum<<endl;

	Eps = 1.0; 			//epsilon value
	Eps_s = "1";		//epsilon string value

	if (argc >= 4){
		Eps = atof(argv[3]);   	//converting string to float
		Eps_s = argv[3];
	}

	cout<<"Epsilon: "<<Eps<<endl;

	if (argc >= 5){
		Eps_1 = atof(argv[4]);   	//converting string to float
		Eps_1_s = argv[4];
	}
	cout<<"Laplace Epsilon: "<<Eps_1<<endl;
	if (argc >= 6){
		Eps_2 = atof(argv[5]);   	//converting string to float
		Eps_2_s = argv[5];
	}
	cout<<"WSLE Epsilon: "<<Eps_2<<endl;

    if (argc >= 7){
		Eps_l = atof(argv[6]);   	//converting string to float
		Eps_l_s = argv[6];
	}

    cout<<"Shuffle Epsilon: "<<Eps_l<<endl;

    if (argc >= 8){
        delta_s = argv[7];
    }
    cout<<"Delta(10^(-)): "<<delta_s<<endl;

	if (argc >= 9){
		c = atof(argv[8]);   	//converting string to float
		c_s = argv[8];
	}
	cout<<"c(deg_th): "<<c<<endl;

	if (argc >= 10){
		t = atoi(argv[9]);   	//converting string to float
		t_s = argv[9];
	}
	if(t==-1 && NodeNum != -1){
		t = (int)floor(NodeNum/3.0);
	}
	t = NodeNum-2;
	cout<<"t: "<<t<<endl;

	ItrNum = 1;			//Number of iterations (set #itr-1 to fix the permutation of nodes)
	fix_perm = 0;
	if (argc >= 11){
		tok  = strtok(argv[10], "-");
		cout<<tok<<endl;
		ItrNum = atoi(tok);
		if((tok  = strtok(NULL, "-")) != NULL){
			if (strcmp(tok, "1") != 0){
				printf("Error: incorrect [#itr(-1)]\n");
				exit(-1);
			}
			else fix_perm = 1;
		}
	}

	cout<<"#Itr: "<<ItrNum<<endl;

	cout<<"Fix Perm: "<<fix_perm<<endl;

	Alg = "ClawShuffle";	

	cout<<"Algorithm: "<<Alg<<endl;

	//Total number of nodes --> all_node_num
	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<2;i++) fgets(s, 1024, fp);
	all_node_num = atoi(s);
	fclose(fp);

	// malloc
	malloc2D(&node_order, ItrNum, all_node_num);

	// Use all nodes
	if (NodeNum == -1){
		NodeNum = all_node_num;
        if(t==-1){
            t = (int)floor(NodeNum/3.0);
			t = NodeNum-2;
            cout<<"t: "<<t<<endl;
        } 
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
	}
	// Randomly generate the order of nodes --> node_order
	else{
		i = EdgeFile.find_last_of("\\");
		outdir = EdgeFile.substr(0, i+1);
		outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
		if(checkFileExistence(outfile)){
			fp = FileOpen(outfile, "r");
			for(j=0;j<all_node_num;j++){
				fgets(s, 1024, fp);
				strtok(s, ",");
				for(i=0;i<ItrNum;i++){
					node_order[i][j] = atoi(strtok(NULL, ","));
				}
			}
			fclose(fp);
		}
		else{
			for(i=0;i<ItrNum;i++){
				MakeRndPerm(node_order[i], all_node_num, all_node_num);
			}
			fp = FileOpen(outfile, "w");
			for(j=0;j<all_node_num;j++){
				fprintf(fp, "%d,", j);
				for(i=0;i<ItrNum;i++) fprintf(fp, "%d,", node_order[i][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}

		// Use only the first permutation
		if (fix_perm){
			for(j=0;j<all_node_num;j++){
				for(i=1;i<ItrNum;i++) node_order[i][j] = node_order[0][j];
			}
		}
	}

	cout<<"NodeNum: "<<NodeNum<<endl;

	malloc1D(&deg, NodeNum);

	// Output the header
	i = EdgeFile.find_last_of("\\");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_eps_1_" + Eps_1_s + "_eps_2_" + Eps_2_s + "_eps_l_" + Eps_l_s + "_delta_" + delta_s + "_c_" + c_s + "_t_" + t_s + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_eps_1_" + Eps_1_s + "_eps_2_" + Eps_2_s + "_eps_l_" + Eps_l_s + "_delta_" + delta_s + "_c_" + c_s + "_t_" + t_s + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "4-Clique(true),4-Clique(emp1),4-Clique(emp2),4-Clique(emp3),4-Clique(clip1),4-Clique(clip2),4-Clique(clip3),4-Clique(rel-err(emp1)),4-Clique(rel-err(emp2)), 4-Clique(rel-err(emp3)),4-Clique(rel-err(clip1)), 4-Clique(rel-err(clip2)),4-Clique(rel-err(clip3))\n");
		fclose(fp);
	}

	cliq_4_num_avg_re1 = 0; cliq_4_num_avg_re2 = 0; cliq_4_num_avg_re3 = 0; clip_cliq_num_avg_re1 = 0; clip_cliq_num_avg_re2 = 0; clip_cliq_num_avg_re3 = 0;
	// For each iteration
	for(itr=0;itr<ItrNum;itr++){
		// Read edges for each iteration when NodeNum < all_node_num
		if(NodeNum < all_node_num || itr == 0){
			// Initialization
			a_mat = new map<int, int>[NodeNum];

			// Read edges from the edge file --> a_mat
			ReadEdges(a_mat, node_order[itr]);

            long long total_edges = 0;

            for(i=0;i<NodeNum;i++){
                for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					total_edges += 1;
					deg[i] += 1;
				}
            }

			total_edges /= 2;
            cout<<"Total edges: "<<total_edges<<endl;

			cliq_num = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					if (i >= j) continue;
					for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
						k = aitr2->first;
						if (j >= k) continue;
						if(a_mat[i].count(k) > 0){
							cliq_num++;
						} 
					}
				}
			}
        }
		cout<<"4-Clique: "<<cliq_num<<endl;

		/************************ Calculate sub-graph counts ************************/

        CalIFCliqCountingShuffler(a_mat, outfile, cliq_4_num1, cliq_4_num2, cliq_4_num3, clip_cliq_num1, clip_cliq_num2, clip_cliq_num3);
        // cout<<"4-Clique counts unbiased: "<<cliq_num_ns<<endl;
        // cout<<"4-Clique counts clipping: "<<cliq_num_clip<<endl;

        clip_cliq_num_re1 = fabs(fabs(clip_cliq_num1) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        clip_cliq_num_re2 = fabs(fabs(clip_cliq_num2) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        clip_cliq_num_re3 = fabs(fabs(clip_cliq_num3) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        cliq_4_num_re1 = fabs(fabs(cliq_4_num1) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        cliq_4_num_re2 = fabs(fabs(cliq_4_num2) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        cliq_4_num_re3 = fabs(fabs(cliq_4_num3) - fabs((double)cliq_num)) / max((double)cliq_num, 0.001 * NodeNum);
        
        /**************************** Output the results ****************************/
        fp = FileOpen(outfile, "a");
        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", 
        (double)cliq_num, cliq_4_num1, cliq_4_num2, cliq_4_num3, clip_cliq_num1, clip_cliq_num2, clip_cliq_num3, cliq_4_num_re1, 
		cliq_4_num_re2, cliq_4_num_re3, clip_cliq_num_re1, clip_cliq_num_re2, clip_cliq_num_re3);
        fclose(fp);

		cliq_4_num_avg_re1 += cliq_4_num_re1;
		cliq_4_num_avg_re2 += cliq_4_num_re2;
		cliq_4_num_avg_re3 += cliq_4_num_re3;
		clip_cliq_num_avg_re1 += clip_cliq_num_re1;
		clip_cliq_num_avg_re2 += clip_cliq_num_re2;
		clip_cliq_num_avg_re3 += clip_cliq_num_re3;

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	cliq_4_num_avg_re1 /= (double)ItrNum;
	cliq_4_num_avg_re2 /= (double)ItrNum;
	cliq_4_num_avg_re3 /= (double)ItrNum;
	clip_cliq_num_avg_re1 /= (double)ItrNum;
	clip_cliq_num_avg_re2 /= (double)ItrNum;
	clip_cliq_num_avg_re3 /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err(emp1)),AVG(rel-error(emp2)),AVG(rel-err(emp3)),AVG(rel-error(clip1)),AVG(rel-err(clip2)),AVG(rel-error(clip3))\n");
	fprintf(fp, "4-Cliques,%e,%e,%e,%e,%e,%e\n", cliq_4_num_avg_re1, cliq_4_num_avg_re2, cliq_4_num_avg_re3, clip_cliq_num_avg_re1, clip_cliq_num_avg_re2, clip_cliq_num_avg_re3);
	fclose(fp);

	// free
	free2D(node_order, ItrNum);
	free1D(deg);

	return 0;
}
