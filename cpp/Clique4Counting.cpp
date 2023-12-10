#include <iostream>
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
double Eps, p2, p1, Eps_l, Eps_1;
string Eps_s, p2_s, Eps_l_s, Eps_1_s, c_s, t_s;

int ItrNum;
int Alg;
double Mu;

int t;
int c;

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

// Calculate non interactive local four clique counts
// Problem is with the nooisy count calculation where we loop through all the vertices
void CalNILocFCliq(map<int, int> *a_mat, string outfile, double &cliq_4_num_ns, double &cliq_4_num_ns_emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	long long cliq_4_num;
	double cliq_4_num_emp;
	double q;
	double rho;
	double rnd;
	int i, j, k, l;

	a_mat_ns = new map<int, int>[NodeNum];

    q = exp(Eps) / (exp(Eps) + 1.0) * p2;

	rho = 1 / exp(Eps); 

	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < q && a_mat[i].count(j) == 0){ //a_mat[i].count(j) gives whether there is an edge between i,j
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd >= q && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	cliq_4_num = 0; cliq_4_num_ns = 0;
	for (i=0; i<NodeNum; i++){
		for (j=i+1; j<NodeNum; j++){
			for (k=j+1; k<NodeNum; k++){
				for (l=k+1; l<NodeNum; l++){
					if(a_mat_ns[i].count(j) > 0 && a_mat_ns[j].count(k) > 0 && a_mat_ns[k].count(l) > 0 && 
					a_mat_ns[l].count(i) > 0 && a_mat_ns[i].count(k) > 0 && a_mat_ns[j].count(l) > 0){
						cliq_4_num_ns++;
					}
					cliq_4_num += (((a_mat_ns[i].count(j)/q)-rho) * ((a_mat_ns[j].count(k)/q)-rho) * ((a_mat_ns[k].count(l)/q)-rho) *
									((a_mat_ns[l].count(i)/q)-rho) * ((a_mat_ns[i].count(k)/q)-rho) * ((a_mat_ns[j].count(l)/q)-rho)) *
									pow((1-rho), -6);
				}
			}
		}
	}
	cliq_4_num_ns = (double)cliq_4_num;
	cliq_4_num_ns_emp = (double)cliq_4_num_emp;
	cout<<cliq_4_num_ns<<" "<<cliq_4_num_ns_emp<<endl;
}

void CalNIFCliqCountingARR(map<int, int> *a_mat, string outfile, double &cliq_4_num_ns, double &cliq_4_num_ns_emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr1;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long cliq_4_num;
	double cliq_4_num_emp;
	long long edge_5, edge_4, edge_3, edge_2, edge_1, edge_0, tri_num, st2_num, total_edges_ns; 
	double cliq_4_num_bs, edge_5_bs, edge_4_bs, edge_3_bs, edge_2_bs, edge_1_bs, edge_0_bs; 
	double rho, murho;
	double rnd;
	int i, j, k, l;
	double q, q2;
	double alp, alp_1_6, q_inv_11, q_inv_21, q_inv_31, q_inv_41, q_inv_51, q_inv_61, q_inv_71;

	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);

	murho = Mu / exp(Eps);

	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < murho && a_mat[i].count(j) == 0){ //a_mat[i].count(j) gives whether there is an edge between i,j
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd >= Mu && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Degree --> deg_ns
	total_edges_ns = 0;
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = 0;
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
		total_edges_ns += deg_ns[i];
	}

	st2_num = 0;
	for(i=0;i<NodeNum;i++){
		st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
	}

	tri_num = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i >= j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num++;
			}
		}
	}

	// cliq_4_num = 0;
	// for (i=0; i<NodeNum; i++){
	// 	for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++){
	// 		j = aitr->first;
	// 		if (i >= j) continue;
	// 		for (aitr1 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++){
	// 			k = aitr1->first;
	// 			if (j >= k) continue;
	// 			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++){
	// 				l = aitr2->first;
	// 				if (k >= l) continue;
	// 				if(a_mat_ns[j].count(k) > 0 && a_mat_ns[k].count(l) > 0 && a_mat_ns[j].count(l) > 0){
	// 					cliq_4_num++;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	cliq_4_num = 0; edge_5 = 0; edge_4 = 0; edge_3 = 0;
	for (i=0; i<NodeNum; i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++){
			j = aitr->first;
			if (i >= j) continue;
			for (aitr1 = a_mat_ns[j].begin(); aitr2 != a_mat_ns[j].end(); aitr2++){
				k = aitr1->first;
				if (j >= k) continue;
				//But actual triangles get counted twice
				edge_3 += deg_ns[i] + deg_ns[j] + deg_ns[k];
				for (aitr2 = a_mat_ns[k].begin(); aitr2 != a_mat_ns[k].end(); aitr2++){
					l = aitr2->first;
					if (k >= l) continue;
					if(a_mat_ns[i].count(k) > 0 && a_mat_ns[j].count(l) > 0 && a_mat_ns[i].count(l) > 0){
						cliq_4_num++;
					}
					else if((a_mat_ns[i].count(k) > 0 && a_mat_ns[j].count(l) > 0  && a_mat_ns[i].count(l) == 0) || 
						(a_mat_ns[i].count(k) > 0 && a_mat_ns[i].count(l) > 0  && a_mat_ns[j].count(l) == 0) || 
						(a_mat_ns[i].count(l) > 0 && a_mat_ns[j].count(l) > 0  && a_mat_ns[i].count(k) == 0)){
							edge_5++;
					}
					else if((a_mat_ns[i].count(k) > 0 && a_mat_ns[j].count(l) == 0  && a_mat_ns[i].count(l) == 0) || 
						(a_mat_ns[i].count(k) == 0 && a_mat_ns[i].count(l) > 0  && a_mat_ns[j].count(l) == 0) || 
						(a_mat_ns[i].count(l) == 0 && a_mat_ns[j].count(l) > 0  && a_mat_ns[i].count(k) == 0)){
							edge_4++;
					}
				}
			}
		}
	}

	//Subtractng the triangles that are counted twice
	edge_3 -= tri_num;

	//Counting 2-edges with respect to 4-clique
	//Selecting any two edges gives 2-edge
	edge_2 = (total_edges_ns*(total_edges_ns-1)-st2_num) + st2_num*(NodeNum-3) - 3*edge_3 - 6*edge_4 - 10*edge_5 - 15*cliq_4_num;

	//No of 4-vertex tuples occuring as single edge
	edge_1 = total_edges_ns*(NodeNum-2)*(NodeNum-3) - 2*edge_2 - 3*edge_3 - 4*edge_4 - 5*edge_5 - 6*cliq_4_num;

	q2 = 1.0 - p2;

	alp = exp(Eps) / (exp(Eps) + 1);
	alp_1_6 = pow(((alp*alp)-1.0), 6);

	cliq_4_num_bs = (double)cliq_4_num / pow(p2, 6);
	edge_5_bs = (double)edge_5 / pow(p2, 5) - 6 * q2 * cliq_4_num_bs;
	edge_4_bs = (double)edge_4 / pow(p2, 4) - 5 * q2 * edge_5_bs - 15 * q2 * q2 * cliq_4_num_bs;
	edge_3_bs = (double)edge_3 / pow(p2, 3) - 4 * q2 * edge_4_bs - 10 * q2 * q2 * edge_5_bs - 20 * q2 * q2 * q2 * cliq_4_num_bs;
	edge_2_bs = (double)edge_2 / pow(p2, 2) - 3 * q2 * edge_3_bs - 6 * q2 * q2 * edge_4_bs - 10 * pow(q2,3) * edge_5_bs - 15 * pow(q2,4) * cliq_4_num_bs;
	edge_1_bs = (double)edge_1 / p2 - 2 * q2 * edge_2_bs - 3 * pow(q2, 2) * edge_3_bs - 4 * pow(q2, 3) * edge_4_bs - 5 * pow(q2,4) * edge_5_bs - 6 * pow(q2,5) * cliq_4_num_bs;
	edge_0_bs = (double)NodeNum*(NodeNum-1)*(NodeNum-2)*(NodeNum-3)/24 - cliq_4_num_bs - edge_5_bs - edge_4_bs - edge_3_bs - edge_2_bs - edge_1_bs;

	q_inv_11 = pow(alp, 6) / alp_1_6;
	q_inv_21 = pow(alp, 5) / alp_1_6;
	q_inv_31 = pow(alp, 4) / alp_1_6;
	q_inv_41 = pow(alp, 3) / alp_1_6;
	q_inv_51 = pow(alp, 2) / alp_1_6;
	q_inv_61 = pow(alp, 1) / alp_1_6;
	q_inv_71 = 1.0 / alp_1_6;

	cliq_4_num_emp = cliq_4_num_bs*q_inv_11 - edge_5_bs*q_inv_21 + edge_4_bs*q_inv_31 - edge_3_bs*q_inv_41 +
						edge_2_bs*q_inv_51 - edge_1_bs*q_inv_61 + edge_0_bs*q_inv_71; 

	cliq_4_num_ns = (double)cliq_4_num;
	cliq_4_num_ns_emp = (double)cliq_4_num_emp;
	cout<<cliq_4_num_ns<<" "<<cliq_4_num_ns_emp<<endl;
}

void CalIFCliqCountingShuffler(map<int, int> *a_mat, string outfile, double &cliq_4_num_ns_unb, double &cliq_4_num_ns_clip){
	int i,j,k,itr,itr1;
	int zij, zjk, zki, zji, zkj, zik, partial_sum;
	double ql, q, q_1_3, ql_1;
	double rnd;
	double star_3;
	int* randperm;
	double cliq_4_num, cliq_4_num_clip;
	int *deg_ns;
	double deg_avg, deg_th;

	map<int, int>:: iterator aitr;

	malloc1D(&randperm, NodeNum);

	malloc1D(&deg_ns, NodeNum);
	
	cout<<NodeNum<<endl;
	MakeRndPerm(randperm, NodeNum, 3*t);

	ql = 1 / (exp(Eps_l) + 1);

	q = 1 / (exp(Eps) + 1);

	q_1_3 = 1 / pow(1-2*ql, 3);

	ql_1 = 1 / (1-2*ql);

	cout<<"Epsilon_L: "<<Eps_l<<endl;
	cout<<"q_L: "<<ql<<endl;
	cout<<"Epsilon: "<<Eps<<endl;
	cout<<"q: "<<q<<endl;
	cout<<"q_1: "<<ql_1<<endl;
	cout<<"ql_1_3: "<<q_1_3<<endl;

	deg_avg = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg_ns[i] += 1;
		deg_ns[i] += stats::rlaplace(0.0, 1/Eps_1, engine);
		deg_avg += deg_ns[i];
	}
	deg_avg = deg_avg / (double) NodeNum;
	cout<<"Degree Avg: "<<deg_avg<<endl;

	deg_th = c * deg_avg;

	cout<<"Degree Threshold: "<<deg_th<<endl;

	cliq_4_num = 0; cliq_4_num_clip = 0;
	for(itr=0; itr<3*t; itr+=3){
		i = randperm[itr];
		j = randperm[itr+1];
		k = randperm[itr+2];
		star_3 = 0;
		for(itr1=0; itr1<NodeNum; itr1++){
			if(itr1 == i || itr1 == j || itr1 == k) continue;
			rnd = genrand_real2();
			star_3 -= ql;
			if(a_mat[itr1].count(i) == 1 && a_mat[itr1].count(j) == 1 && a_mat[itr1].count(k) == 1){
				if(rnd < (1-ql)) star_3++;
			}
			else if (rnd < ql) star_3++; 
		}
		rnd = genrand_real2();
		if ((rnd < q && a_mat[i].count(j) == 0) || (rnd < (1-q) && a_mat[i].count(j) == 1)) zij = 1;
		else zij = 0;
		if ((rnd < q && a_mat[j].count(k) == 0) || (rnd < (1-q) && a_mat[j].count(k) == 1)) zjk = 1;
		else zjk = 0;
		if ((rnd < q && a_mat[k].count(i) == 0) || (rnd < (1-q) && a_mat[k].count(i) == 1)) zki = 1;
		else zki = 0;

		if ((rnd < q && a_mat[j].count(i) == 0) || (rnd < (1-q) && a_mat[j].count(i) == 1)) zji = 1;
		else zji = 0;
		if ((rnd < q && a_mat[k].count(j) == 0) || (rnd < (1-q) && a_mat[k].count(j) == 1)) zkj = 1;
		else zkj = 0;
		if ((rnd < q && a_mat[i].count(k) == 0) || (rnd < (1-q) && a_mat[i].count(k) == 1)) zik = 1;
		else zik = 0;

		partial_sum = (zij-q) * (zjk-q) * (zki-q) + (zij-q) * (zjk-q) * (zik-q) + 
					(zij-q) * (zkj-q) * (zki-q) + (zij-q) * (zkj-q) * (zik-q) + 
					(zji-q) * (zjk-q) * (zki-q) + (zji-q) * (zjk-q) * (zik-q) + 
					(zji-q) * (zkj-q) * (zki-q) + (zji-q) * (zkj-q) * (zik-q);

		cliq_4_num += partial_sum * star_3 * ql_1 * q_1_3 / 8.0;
		if(min(min(deg_ns[i], deg_ns[j]), deg_ns[k]) > deg_th) cliq_4_num_clip += partial_sum * star_3 * ql_1 * q_1_3 / 8.0;
	} 

	cliq_4_num = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num / (24*t);

	cliq_4_num_clip = NodeNum*(NodeNum-1)*(NodeNum-2)*cliq_4_num_clip / (24*t);

	cout<<"4-Clique counts unbiased: "<<cliq_4_num<<endl;
	cout<<"4-Clique counts clipping: "<<cliq_4_num_clip<<endl;

	cliq_4_num_ns_unb = cliq_4_num;

	cliq_4_num_ns_clip = cliq_4_num_clip;

	free1D(deg_ns);
	free1D(randperm);
}

// KL divergence function
double CalcKL(double p1, double p2){
	return p1 * log(p1 / p2) + (1.0 - p1) * log((1.0 - p1) / (1.0 - p2));
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int **node_order;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<int, int>::iterator aitr3;

	long long cliq_4_num;

	double cliq_4_num_ns, cliq_4_num_ns_emp, cliq_4_num_re, cliq_4_num_avg_re;
	double cliq_4_num_ns_unb, cliq_4_num_ns_clip;
	
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
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon-mu/mu2/mu3 (default: 1-1)] [NSType (default: -1)] [tclip-eclip (default: -1)] [#itr(-1) (default: 1)] [alg (default: 0)] [Balloc (default: 1-1)])\n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon-mu/mu2/mu3]: Parameters epsilon and mu/mu2/mu3 (I/II/III) (mu/mu2/mu3 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1)\n");
		printf("[NSType]: Noise type (-1: no noise, 0: Lap (max degree), 1: Lap (true degree + clip), 2: Lap (double clip: noisy degree + clip))\n");
		printf("[tclip(-eclip)]: Triangle and edge clipping parameters (-1: no clipping) (alg=2-4; set eclip when NSType=2)\n");
		printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
		printf("[alg]: Algorithm (1: interactive local, 2: efficient interactive local I, 3: efficient interactive local II, 4: efficient interactive local III, 5: non-interactive local (RR w/ emp), 6: non-interactive local (RR w/o emp), 7: [Ye+, T-KDE (mean)], 8: [Ye+, T-KDE (median)], 9: [Ye+, T-KDE (most frequent degree)], 10: non-interactive local (ARR w/ emp))\n");
		printf("[Balloc]: Privacy budget allocation (alg=1-3): Eps1st-Eps2ndTrSt\n");
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
		Eps_l = atof(argv[4]);   	//converting string to float
		Eps_l_s = argv[4];
	}
	cout<<"Permute Epsilon: "<<Eps_l<<endl;
	if (argc >= 6){
		Eps_1 = atof(argv[5]);   	//converting string to float
		Eps_1_s = argv[5];
	}
	cout<<"Laplace Epsilon: "<<Eps_1<<endl;
	if (argc >= 7){
		c = atof(argv[6]);   	//converting string to float
		c_s = argv[6];
	}
	cout<<"C(deg_th): "<<c<<endl;
	if (argc >= 8){
		t = atoi(argv[7]);   	//converting string to float
		t_s = argv[7];
	}
	if(t==-1){
		t = (int)floor(NodeNum/3.0);
	}
	cout<<"T: "<<t<<endl;

    //Sampling Prob
    if (argc >= 9){
		p2 = pow(10,-atof(argv[8]));   	//converting string to float
		p2_s = argv[8];
	}

    cout<<"Sampling Prob: "<<p2<<endl;

	ItrNum = 1;			//Number of iterations (set #itr-1 to fix the permutation of nodes)
	fix_perm = 0;
	if (argc >= 10){
		tok  = strtok(argv[9], "-");
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

	// We require 1,2,3,4 algorithms
	Alg = 0;	//Algorithm (1: interactive local, 2: efficient interactive local I, 3: efficient interactive local II,
				//4: efficient interactive local III, 5: non-interactive local (RR w/ emp), 6: non-interactive local (RR w/o emp), 
				//7: [Ye+, T-KDE (mean)], 8: [Ye+, T-KDE (median)], 9: [Ye+, T-KDE (most frequent degree)], 10: non-interactive local (ARR w/ emp))
	if (argc >= 11) Alg = atoi(argv[10]);
	// if (Alg <= 0 || Alg > 10){
	// 	printf("Error: incorrect [Alg]\n");
	// 	exit(-1);
	// }

	cout<<"Algorithm: "<<Alg<<endl;


    //Setting RR
    p1 = exp(Eps) / (exp(Eps) + 1.0);
    
    cout<<"RR Prob: "<<p1<<endl;

    //Setting Mu
    Mu = p1*p2;

    cout<<"MU: "<<Mu<<endl;
	// Total number of nodes --> all_node_num
	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<2;i++) fgets(s, 1024, fp);
	all_node_num = atoi(s);
	fclose(fp);

	// malloc
	malloc2D(&node_order, ItrNum, all_node_num);

	// Use all nodes
	if (NodeNum == -1){
		NodeNum = all_node_num;
		t = (int)floor(NodeNum/3.0);
		cout<<"T_all: "<<t<<endl;
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

	cout<<NodeNum<<endl;

	//Output the header
	i = EdgeFile.find_last_of("\\");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "#4-cliq(true),#4-cliq(est),#4-cliq(emp-est),#4-cliq(rel-err),#4-cliq(l2-loss)\n");
		fclose(fp);
	}

	// For each iteration
	for(itr=0;itr<ItrNum;itr++){
		// Read edges for each iteration when NodeNum < all_node_num
		if(NodeNum < all_node_num || itr == 0){
			// Initialization
			a_mat = new map<int, int>[NodeNum];

			// Read edges from the edge file --> a_mat
			ReadEdges(a_mat, node_order[itr]);

			cliq_4_num = 0;
			for (i=0; i<NodeNum; i++){
				for (aitr = a_mat[i].begin(); aitr!=a_mat[i].end(); aitr++){
					j = aitr->first;
					if (i>=j) continue;
					for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++){
						k = aitr2->first;
						if (j>=k) continue;
						for (aitr3 = a_mat[k].begin(); aitr3 != a_mat[k].end(); aitr3++){
							l = aitr3->first;
							if (k>=l) continue;
							if(a_mat[l].count(i) > 0 && a_mat[l].count(j) > 0 && a_mat[i].count(k) > 0){
								cliq_4_num++;
							}
						}
					}
				}
			}

		}
		cout<<"Clique-4: "<<cliq_4_num<<endl;

		/************************ Calculate sub-graph counts ************************/

		if (Alg == 1){
			if (NodeNum <= 10000) CalNILocFCliq(a_mat, outfile, cliq_4_num_ns, cliq_4_num_ns_emp);
			cout<<"Main 4-cliques Noisy: "<<cliq_4_num_ns<<endl;
			cout<<"Unbiased 4-cliques Noisy: "<<cliq_4_num_ns_emp<<endl;

			cliq_4_num_re = fabs(cliq_4_num_ns_emp - (double)cliq_4_num) / max((double)cliq_4_num, 0.001 * NodeNum);

			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			fprintf(fp, "%lld,%e,%e, %e\n", 
			cliq_4_num, cliq_4_num_ns, cliq_4_num_ns_emp, cliq_4_num_re);
			fclose(fp);
		}
		else if(Alg == 2){
			cout<<"Calling function"<<endl;
			CalIFCliqCountingShuffler(a_mat, outfile, cliq_4_num_ns_unb, cliq_4_num_ns_clip);
			cout<<"4-cliques counts unbiased: "<<cliq_4_num_ns_unb<<endl;
			cout<<"4-cliques counts clipping: "<<cliq_4_num_ns_clip<<endl;

			cliq_4_num_re = fabs(cliq_4_num_ns_clip - (double)cliq_4_num) / max((double)cliq_4_num, 0.001 * NodeNum);
			
			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			fprintf(fp, "%lld,%e,%e,%e\n", 
			cliq_4_num, cliq_4_num_ns_unb, cliq_4_num_ns_clip, cliq_4_num_re);
			fclose(fp);
		}
		else if(Alg == 3){
			CalNIFCliqCountingARR(a_mat, outfile, cliq_4_num_ns, cliq_4_num_ns_emp);
			cout<<"Main 4-cliques Noisy: "<<cliq_4_num_ns<<endl;
			cout<<"Unbiased 4-cliques Noisy: "<<cliq_4_num_ns_emp<<endl;

			cliq_4_num_re = fabs(cliq_4_num_ns_emp - (double)cliq_4_num) / max((double)cliq_4_num, 0.001 * NodeNum);

			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			fprintf(fp, "%lld,%e,%e, %e\n", 
			cliq_4_num, cliq_4_num_ns, cliq_4_num_ns_emp, cliq_4_num_re);
			fclose(fp);
		}

		cliq_4_num_avg_re += cliq_4_num_re;

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	cliq_4_num_avg_re /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err)\n");
	fprintf(fp, "4-Cliques,%e\n", cliq_4_num_avg_re);
	fclose(fp);

	cout<<"Averaged Relative Error: "<<cliq_4_num_avg_re<<endl;

	// free
	free2D(node_order, ItrNum);

	return 0;
}
