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
int t;
double Mu;

double c;

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

//Calculate the triangle counts with applying the privacy amplification by shuffling
//While building noisy adj matrix, you do random permutation for each user and then apply 
//ARR to that bit to get noisy bit
void CalcNLocTriARRPerm(map<int, int> *a_mat, string outfile, double &tri_num_ns_emp, double &tri_num_ns){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long tot_edge_num_ns;
	double tri_num, st2_num, ed2_num, ed1_num, non_num;
	double tri_num_bs, ed2_num_bs, ed1_num_bs, non_num_bs;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int rnd_ind, curr_ind;
	int i, j, k;
	double murho, q2;
	int* randperm;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&randperm, NodeNum);

	MakeRndPerm(randperm, NodeNum, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps_l);
    cout<<"Murhho: "<<murho<<endl;
	
    cout<<"Mu: "<<Mu<<endl;

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	for(i=0;i<NodeNum;i++){
		if(i%1000==0){
			cout<<"Heart Beat: "<<i<<endl;
		}
		curr_ind = 0;
		for(j=i+1;j<NodeNum;j++){
			rnd_ind = curr_ind + (int)(genrand_int32() % (NodeNum-curr_ind));
			swap(randperm, curr_ind, rnd_ind);
			rnd = genrand_real2();
			if(rnd < murho && a_mat[i].count(randperm[curr_ind]) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd < Mu && a_mat[i].count(curr_ind) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			curr_ind++;
		}
	}

	// Degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
	// Here we consider only upper triangle matrix of a_mat to compute a_mat_ns
	// So I think we should not divide it by 2 for total edges
	tot_edge_num_ns /= 2;
    cout<<"Noisy Edges: "<<tot_edge_num_ns<<endl;
	// #triangles --> tri_num
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

    cout<<"Triangle Noisy: "<<tri_num<<endl;

	// With empirical estimation
	// #2-stars --> st2_num
	st2_num = 0;
	for(i=0;i<NodeNum;i++){
		st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
	}

	// #2-edges --> ed2_num
	ed2_num = st2_num - 3*tri_num;
	// #1-edge --> ed1_num
	ed1_num = (long long)tot_edge_num_ns*(NodeNum-2) - 2*ed2_num - 3*tri_num;

	// Calculate #triangles, #2-edges, #1-edge before sampling --> tri_num_bs, ed2_num_bs, ed1_num_bs
	q2 = 1.0 - p2;
	tri_num_bs = (double)tri_num / (p2 * p2 * p2);
	ed2_num_bs = (double)ed2_num / (p2 * p2) - 3.0 * q2 * tri_num_bs;
	ed1_num_bs = (double)ed1_num / p2 - 2.0 * q2 * ed2_num_bs - 3.0 * q2 * q2 * tri_num_bs;

	// #none --> non_num_bs
	non_num_bs = (double)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num_bs - ed2_num_bs - ed1_num_bs;

	alp = exp(Eps);
	alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
	q_inv_11 = (alp*alp*alp) / alp_1_3;
	q_inv_21 = - alp*alp / alp_1_3;
	q_inv_31 = alp / alp_1_3;
	q_inv_41 = - 1.0 / alp_1_3;

	tri_num_ns_emp = tri_num_bs * q_inv_11 + ed2_num_bs * q_inv_21 + ed1_num_bs * q_inv_31 + non_num_bs * q_inv_41;

	// Without empirical estimation
	tri_num_ns = (double)tri_num;

	cout<<"unbiased Noisy emp est: "<<tri_num_ns_emp<<endl;
	cout<<"Noisy emp est: "<<tri_num_ns<<endl;

	delete[] a_mat_ns;
	free(randperm);
	free1D(deg_ns);
}

void CalITriangleCountingShuffler(map<int, int> *a_mat, string outfile, double &tri_num_ns_unb, double &tri_num_ns_clip){
	int i,j,itr,itr1;
	int* randperm;	
	double ql, q, ql_1, q_1;
	int zi, zj;
	double star_2;
	double tri_num, tri_num_clip;
	double rnd;
	int *deg_ns;
	double deg_avg, deg_th;
	
	map<int, int>:: iterator aitr;

	malloc1D(&randperm, NodeNum);

	malloc1D(&deg_ns, NodeNum);
	
	MakeRndPerm(randperm, NodeNum, 2*t);


	ql = 1 / (exp(Eps_l) + 1);

	q = 1 / (exp(Eps) + 1);

	ql_1 = 1 / (1-2*ql);

	q_1 = 1 / (1-2*q);

	cout<<"Epsilon_L: "<<Eps_l<<endl;
	cout<<"q_L: "<<ql<<endl;
	cout<<"Epsilon_1: "<<Eps_1<<endl;
	cout<<"q_1: "<<q_1<<endl;
	cout<<"Epsilon: "<<Eps<<endl;
	cout<<"q: "<<q<<endl;
	
	deg_avg = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg_ns[i] += 1;
		deg_ns[i] += stats::rlaplace(0.0, 1/Eps_1, engine);
		deg_avg += deg_ns[i];
	}
	deg_avg = deg_avg / (double) NodeNum;

	deg_th = c * deg_avg;

	cout<<"Degree Threshold: "<<deg_th<<endl;

	tri_num = 0; tri_num_clip = 0;
	for(itr=0; itr<2*t; itr+=2){
		i = randperm[itr];
		j = randperm[itr+1];
		star_2 = 0;
		for(itr1=0; itr1<NodeNum; itr1++){
			if(itr1 == i || itr1 == j) continue;
			rnd = genrand_real2();
			star_2 -= ql;
			if(a_mat[itr1].count(i) == 1 && a_mat[itr1].count(j) == 1){
				if(rnd < (1-ql)) star_2++;
			}
			else if(rnd < ql) star_2++;
		}
		// cout<<"i,j,star_2: "<<i<<","<<j<<","<<star_2<<endl;
		rnd = genrand_real2();
		if ((rnd < q && a_mat[i].count(j) == 0) || (rnd < (1-q) && a_mat[i].count(j) == 1)) zi = 1;
		else zi = 0;
		if ((rnd < q && a_mat[j].count(i) == 0) || (rnd < (1-q) && a_mat[j].count(i) == 1)) zj = 1;
		else zj = 0;
		
		tri_num += (zi + zj - 2*q) * star_2 * q_1 * ql_1/ 2.0;
		if(min(deg_ns[i], deg_ns[j]) > deg_th) tri_num_clip += (zi + zj - 2*q) * star_2 * q_1 * ql_1 / 2.0;
	}

	tri_num = NodeNum*(NodeNum-1)*tri_num / (6*t);

	tri_num_clip = NodeNum*(NodeNum-1)*tri_num_clip / (6*t);

	cout<<"Triangle counts unbiased: "<<tri_num<<endl;
	cout<<"Triangle counts clipping: "<<tri_num_clip<<endl;

	tri_num_ns_unb = tri_num;

	tri_num_ns_clip = tri_num_clip;
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

	long long tri_num;

	double tri_num_ns, tri_num_ns_emp, tri_num_ns_unb, tri_num_ns_clip, tri_num_re, tri_num_avg_re;
	
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
		t = (int)ceil(NodeNum/2.0);
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
    p1 = exp(Eps_l) / (exp(Eps_l) + 1.0);
    
    cout<<"RR Prob: "<<p1<<endl;

    //Setting Mu
    Mu = p1*p2;

    cout<<"MU: "<<Mu<<endl;

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
		t = (int)ceil(NodeNum/2.0);
		cout<<"T_all: "<<t<<endl;
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
	}
	// Randomly generate the order of nodes --> node_order
	else{
		i = EdgeFile.find_last_of("/");
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

	// cout<<"Nodeorder: "<<endl;

	// for(j=0;j<NodeNum;j++) cout<<node_order[0][j]<<endl;

	// Output the header
	i = EdgeFile.find_last_of("/");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "SampProb" + p2_s + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + to_string(Alg) + "_eps" + Eps_s + "SampProb" + p2_s + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "Triangle(true),Triangle(est),Triangle(emp-est),Triangle(rel-err),Triangle(l2-loss)\n");
		fclose(fp);
	}

	tri_num_avg_re = 0;
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
                for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) total_edges += 1;
            }

			// cout<<"Adjacency List: "<<endl;

			// for(i=0;i<NodeNum;i++){
            //     for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++){
			// 		cout<<"Node: "<<i<<" Key: "<<aitr->first<<" Value: "<<aitr->second<<endl;
			// 	};
            // }

            cout<<"Total edges: "<<total_edges<<endl;

			tri_num = 0;
			for(i=0;i<NodeNum;i++){
				// cout<<i<<" Node"<<endl;
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					// cout<<"Iter:1 "<<aitr->first<<" "<<aitr->second<<endl;
					if (i >= j) continue;
					for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
						k = aitr2->first;
						// cout<<"Iter:2 "<<aitr2->first<<" "<<aitr2->second<<endl;
						if (j >= k) continue;
						// cout<<"Count: "<<a_mat[j].count(k)<<endl;
						if(a_mat[i].count(k) > 0){
							tri_num++;
							// if(tri_num%100000 == 0){
							// 	cout<<"Heart Beat: "<<tri_num<<endl;
							// }
						} 
					}
				}
			}
        }
		cout<<"Triangle: "<<tri_num<<endl;

		/************************ Calculate sub-graph counts ************************/

		if (Alg == 1){
			CalcNLocTriARRPerm(a_mat, outfile, tri_num_ns_emp, tri_num_ns);
			cout<<"Main Triangle Noisy: "<<tri_num_ns<<endl;
			cout<<"Unbiased Triangle Noisy: "<<tri_num_ns_emp<<endl;

			tri_num_re = fabs(tri_num_ns - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);

			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			fprintf(fp, "%lld,%e,%e,%e\n", 
			tri_num, tri_num_ns, tri_num_ns_emp, tri_num_re);
			fclose(fp);
		}
		else if (Alg == 2){
			CalITriangleCountingShuffler(a_mat, outfile, tri_num_ns_unb, tri_num_ns_clip);
			cout<<"Triangle counts unbiased: "<<tri_num_ns_unb<<endl;
			cout<<"Triangle counts clipping: "<<tri_num_ns_clip<<endl;

			tri_num_re = fabs(tri_num_ns_clip - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
			
			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			fprintf(fp, "%lld,%e,%e,%e\n", 
			tri_num, tri_num_ns_unb, tri_num_ns_clip, tri_num_re);
			fclose(fp);
		}

		tri_num_avg_re += tri_num_re;

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	tri_num_avg_re /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err)\n");
	fprintf(fp, "Triangles,%e\n", tri_num_avg_re);
	fclose(fp);

	cout<<"Averaged Relative Error: "<<tri_num_avg_re<<endl;
	// free
	free2D(node_order, ItrNum);

	return 0;
}
