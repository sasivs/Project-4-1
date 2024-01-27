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
double Eps, Mu, p2, p2_flag;
string Eps_s, p2_s, p2_flag_s;

int ItrNum;
string Alg;

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
				// But actual triangles get counted twice
				// edge_3 += deg_ns[i] + deg_ns[j] + deg_ns[k];
				// if(a_mat_ns[k].count(i)>0){
				// 	edge_3 += NodeNum-3;
				// }
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
					else if(a_mat_ns[i].count(k) == 0 && a_mat_ns[j].count(l) == 0 && a_mat_ns[i].count(l) == 0){
						edge_3++;
					}
				}
			}
		}
	}

	//Adding the triangles that are not counted while counting 4-cliques
	edge_3 += tri_num;

	//Counting 2-edges with respect to 4-clique
	//Selecting any two edges gives 2-edge
	edge_2 = (total_edges_ns*(total_edges_ns-1)-st2_num) + st2_num*(NodeNum-3) - 3*edge_3 - 6*edge_4 - 10*edge_5 - 15*cliq_4_num;

	//No of 4-vertex tuples occuring as single edge
	edge_1 = total_edges_ns*(NodeNum-2)*(NodeNum-3) - 2*edge_2 - 3*edge_3 - 4*edge_4 - 5*edge_5 - 6*cliq_4_num;

	q2 = 1.0 - p2;

	alp = exp(Eps);
	alp_1_6 = pow((alp-1.0), 6);

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


int main(int argc, char *argv[])
{
	int all_node_num;
	int **node_order;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<int, int>::iterator aitr3;
	int *deg;									// degree
	long long tot_edge_num;
	long long cliq_num;
	double cliq_num_ns, cliq_num_emp;
	double cliq_re_ns, cliq_re_ns_avg, cliq_re_emp, cliq_re_emp_avg;
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
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon (default: 1)] [Sampling Probability-base ] [#itr(-1) (default: 1)] \n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon]: Parameters epsilon \n");
		printf("[Sampling Probability-base]: Parameters p2 (10^-p2/NodeNumm^-p2) \n");
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
		Eps_s = argv[3];			//String values
		Eps = atof(argv[3]);   	//converting string to float
	}

	cout<<"Epsilon: "<<Eps<<endl;

	p2_flag = 0.0;
	if (argc >= 5){
		tok  = strtok(argv[4], "-");
		p2_s = tok;
		p2 = stod(p2_s);
		if((tok  = strtok(NULL, "-")) != NULL){
			p2_flag_s = tok;
			p2_flag = stod(p2_flag_s);
		}
	}

	if (p2_flag == 0.0){
		p2 = pow(10, -p2);
	}
	else if (p2_flag == 1.0 && NodeNum != -1){
		p2 = pow(NodeNum, -1/p2);
	}

	cout<<"Sampling Probability: "<<p2<<endl;

	ItrNum = 1;			//Number of iterations (set #itr-1 to fix the permutation of nodes)
	fix_perm = 0;
	if (argc >= 6){
		tok  = strtok(argv[5], "-");
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

	Alg = "LocalARRClique";

	cout<<"Algorithm: "<<Alg<<endl;

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
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
		cout<<"NodeNum: "<<NodeNum<<endl;
		if (p2_flag == 1.0) p2 = pow(NodeNum, -1/p2);
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

	// #triplet --> triplet_num

	// Initialization
	malloc1D(&deg, NodeNum);
	cliq_re_ns_avg = cliq_re_emp_avg = 0.0;	//tri=triangle re=relative error, l2=l2 loss

	// Output the header
	i = EdgeFile.find_last_of("/");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_p2_" + p2_s + "_p2flag_" + p2_flag_s + "_itr_" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_p2_" + p2_s + "_p2flag_" + p2_flag_s + "_itr_" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp,"4-Clique(true),4-Clique(est),4-Clique(emp-est),4-Clique(rel-err(noisy)),4-Clique(rel-error(emp)))\n");
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

			// Degree --> deg
			//Find degree of each node
			for(i=0;i<NodeNum;i++) deg[i] = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg[i] += 1;
			}

			// Total number of edges --> tot_edge_num
			tot_edge_num = 0;
			for(i=0;i<NodeNum;i++) tot_edge_num += (long long)deg[i];
			tot_edge_num /= 2;

			// #triangles --> tri_num
			//counting triangles
			cliq_num = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					if (i >= j) continue;
					for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
						k = aitr2->first;
						if (j >= k) continue;
						if(a_mat[j].count(k) > 0) cliq_num++;
					}
				}
			}
		}
		/************************ Calculate sub-graph counts ************************/
		// Calculate #triangles
		CalNIFCliqCountingARR(a_mat, outfile, cliq_num_ns, cliq_num_emp);

		/**************************** Evaluate the loss *****************************/
		// relative error --> cliq_re_ns
		cliq_re_ns = fabs(cliq_num_ns - (double)cliq_num) / max((double)cliq_num, 0.001 * NodeNum);
		cliq_re_ns_avg += cliq_re_ns;
		cliq_re_emp = fabs(cliq_num_emp - (double)cliq_num) / max((double)cliq_num, 0.001 * NodeNum);
		cliq_re_emp_avg += cliq_re_emp;

		/**************************** Output the results ****************************/
		fp = FileOpen(outfile, "a");
		fprintf(fp, "%e,%e,%e,%e,%e\n", 
		(double)cliq_num, cliq_num_ns, cliq_num_emp, cliq_re_ns, cliq_re_emp);
		fclose(fp);

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	/************************* Output the results (AVG) *************************/
	cliq_re_ns_avg /= (double)ItrNum;
	cliq_re_emp_avg /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err(noisy)),AVG(rel-error(emp))\n");
	fprintf(fp, "4-Clique,%e,%e\n", cliq_re_ns_avg, cliq_re_emp_avg);
	fclose(fp);

	// free
	free2D(node_order, ItrNum);
	free1D(deg);

	return 0;
}
