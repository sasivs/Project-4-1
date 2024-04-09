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
double Eps, Eps_l, Eps_1, Eps_2, p2, p1, p2_flag;
string Eps_s, p2_s, Eps_l_s, Eps_1_s, Eps_2_s, delta_s, p2_flag_s;

int ItrNum;
string Alg;
double Mu;


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

void CalcNLocShuffleTri(map<int, int> *a_mat, int* deg, string outfile, double &tri_num_ns, double &tri_num_samp, double &tri_num_emp1, double &tri_num_emp2){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns, *deg_clip;									// noisy degree
	long long tot_edge_num_ns;
	double tri_num, st2_num, ed2_num, ed1_num, non_num;
	double tri_num_bs, ed2_num_bs, ed1_num_bs, non_num_bs;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int rnd_ind, curr_ind;
	int i, j, k;
	double murho, q2, murho_bar, Mu_bar;
	double murho_emp2, murho_bar_emp2, Mu_bar_emp2, Mu_emp2, p1_emp2;
	int* randperm;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&deg_clip, NodeNum);
	malloc1D(&randperm, NodeNum);

	MakeRndPerm(randperm, NodeNum, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
    p1 = exp(Eps_l)/(exp(Eps_l)+1);
    Mu = p1*p2;
	murho = Mu / exp(Eps_l);
    murho_bar = 1 - murho;
    Mu_bar = 1 - Mu;

    cout<<"Murhho: "<<murho<<endl;
    cout<<"Mu: "<<Mu<<endl;

	// Flip 0/1 in a_mat with probability Mu --> a_mat_ns
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
                deg_ns[i]+=1;
				deg_ns[j]+=1;
			}
			// 1 --> 1 (not flip)
			else if(rnd < Mu && a_mat[i].count(randperm[curr_ind]) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
                deg_ns[i]+=1;
				deg_ns[j]+=1;
			}
			curr_ind++;
		}
		deg_clip[i] = max((int)floor((double)deg[i] + stats::rlaplace(0.0, 1.0/Eps_1, engine)), 0);
	}

    // Graph Projection
    for (i=0; i<NodeNum; i++){
		if(i%1000==0){
			cout<<"Heart Beat: "<<i<<endl;
		}
		if(deg_ns[i]>deg_clip[i]){
			MakeRndPerm(randperm, deg_ns[i], (deg_ns[i]-deg_clip[i]));

			sort(randperm, randperm+(deg_ns[i]-deg_clip[i]));

			int iter_ind = 0, next_ind = 0, key;
			aitr = a_mat_ns[i].begin();
			while(aitr != a_mat_ns[i].end()){
				if(iter_ind == randperm[next_ind]){
					key = aitr->first;
					aitr++;
					a_mat_ns[i].erase(key);
					a_mat_ns[key].erase(i);
					deg_ns[i]-=1;
					deg_ns[key]-=1;
					next_ind++;
				}
				else{
					aitr++;
				}
				iter_ind++;
			}
		}
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
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
	tri_num_ns = (double)tri_num;

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
    non_num = (double)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num - ed2_num - ed1_num;

	// Calculate #triangles, #2-edges, #1-edge before sampling --> tri_num_bs, ed2_num_bs, ed1_num_bs
	tri_num_bs = (double)tri_num / (p2 * p2 * p2);
    tri_num_samp = tri_num_bs;

	q2 = 1.0 - p2;
	ed2_num_bs = (double)ed2_num / (p2 * p2) - 3.0 * q2 * tri_num_bs;
	ed1_num_bs = (double)ed1_num / p2 - 2.0 * q2 * ed2_num_bs - 3.0 * q2 * q2 * tri_num_bs;

	// #none --> non_num_bs
	non_num_bs = (double)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num_bs - ed2_num_bs - ed1_num_bs;

	alp = exp(Eps_2);
	alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
	q_inv_11 = (alp*alp*alp) / alp_1_3;
	q_inv_21 = - alp*alp / alp_1_3;
	q_inv_31 = alp / alp_1_3;
	q_inv_41 = - 1.0 / alp_1_3;

	tri_num_emp1 = tri_num_bs * q_inv_11 + ed2_num_bs * q_inv_21 + ed1_num_bs * q_inv_31 + non_num_bs * q_inv_41;

	// cout<<"Emp1: "<<tri_num_emp1<<endl;
	// cout<<"Triangle coeff: "<<q_inv_11<<endl;
	// cout<<"Triangle term: "<<tri_num_bs * q_inv_11<<endl;
	// cout<<"2-edges coeff: "<<q_inv_21<<endl;
	// cout<<"2-edges term: "<<ed2_num_bs * q_inv_21<<endl;
	// cout<<"1-edges coeff: "<<q_inv_31<<endl;
	// cout<<"1-edges term: "<<ed1_num_bs * q_inv_31<<endl;
	// cout<<"no-edges coeff: "<<q_inv_41<<endl;
	// cout<<"no-edges term: "<<non_num_bs * q_inv_41<<endl;

    p1_emp2 = exp(Eps_2)/(exp(Eps_2)+1);
    Mu_emp2 = p1_emp2*p2;
	murho_emp2 = Mu_emp2 / exp(Eps_2);
    murho_bar_emp2 = 1 - murho_emp2;
    Mu_bar_emp2 = 1 - Mu_emp2;

    alp_1_3 = pow((Mu_bar_emp2*murho_emp2 - Mu_emp2*murho_bar_emp2), 3);
	q_inv_11 = - (pow(murho_bar_emp2, 3) / alp_1_3);
	q_inv_21 = murho_emp2*pow(murho_bar_emp2, 2) / alp_1_3;
	q_inv_31 = - (pow(murho_emp2, 2)*murho_bar_emp2 / alp_1_3);
	q_inv_41 = pow(murho_emp2, 3) / alp_1_3;

	tri_num_emp2 = tri_num * q_inv_11 + ed2_num * q_inv_21 + ed1_num * q_inv_31 + non_num * q_inv_41;
	// cout<<"Emp2: "<<tri_num_emp2<<endl;
	// cout<<"Triangle coeff: "<<q_inv_11<<endl;
	// cout<<"Triangle term: "<<tri_num * q_inv_11<<endl;
	// cout<<"2-edges coeff: "<<q_inv_21<<endl;
	// cout<<"2-edges term: "<<ed2_num * q_inv_21<<endl;
	// cout<<"1-edges coeff: "<<q_inv_31<<endl;
	// cout<<"1-edges term: "<<ed1_num * q_inv_31<<endl;
	// cout<<"no-edges coeff: "<<q_inv_41<<endl;
	// cout<<"no-edges term: "<<non_num * q_inv_41<<endl;

	delete[] a_mat_ns;
	free(randperm);
	free1D(deg_ns);
	free1D(deg_clip);
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int **node_order;
	int *deg;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;

	long long tri_num;

	double tri_num_ns, tri_num_samp_emp, tri_num_emp1, tri_num_emp2, tri_num_ns_re, tri_num_ns_avg_re, tri_num_samp_emp_re, tri_num_emp1_re, tri_num_emp2_re,
	tri_num_samp_emp_avg_re, tri_num_emp1_avg_re, tri_num_emp2_avg_re;
	
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
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon] [epsilon-1] [epsilon-2] [epsilon-l] [delta] [Sampling Probability-base] [#itr(-1) (default: 1)] \n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
        printf("[epsilon]: Parameters epsilon \n");
        printf("[epsilon-1]: Parameters epsilon-1 (Lap Noise Parameter)\n");
		printf("[epsilon-2]: Parameters epsilon-2 (Wedge shuffle alg parameter to calculate ep-l)\n");
		printf("[epsilon-l]: Parameters epsilon-l (Wedge shuffle amplification parameter)\n");
		printf("[delta]: Parameters delta (Wedge shuffle alg parameter to calculate ep-l)\n");
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
	cout<<"RR Epsilon: "<<Eps_2<<endl;

	if (argc >= 7){
		Eps_l = atof(argv[6]);   	//converting string to float
		Eps_l_s = argv[6];
	}
	cout<<"Shuffle Epsilon: "<<Eps_1<<endl;

    if (argc >= 8){
        delta_s = argv[7];
    }

    p2_flag = 0.0;
	if (argc >= 9){
		tok  = strtok(argv[8], "-");
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

	Alg = "LocalShuffle";
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
		cout<<"NodeNum: "<<NodeNum<<endl;
        if (p2_flag == 1.0) p2 = pow(NodeNum, -1/p2);
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

	malloc1D(&deg, NodeNum);

	// Output the header
	i = EdgeFile.find_last_of("\\");
	outdir = EdgeFile.substr(0, i+1);
	for(i=0;i<3;i++){
		if(fix_perm) outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_eps_1_" + Eps_1_s + "_eps_2_" + Eps_2_s + "_eps_l_" + Eps_l_s + "_delta_" + delta_s + "_p2_" + p2_s + "_p2flag_" + p2_flag_s + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n_" + to_string(NodeNum) + "_alg_" + Alg + "_eps_" + Eps_s + "_eps_1_" + Eps_1_s + "_eps_2_" + Eps_2_s + "_eps_l_" + Eps_l_s + "_delta_" + delta_s + "_p2_" + p2_s + "_p2flag_" + p2_flag_s + "_itr" + to_string(ItrNum) + ".csv";
		fp = FileOpen(outfile, "w");
		fprintf(fp, "Triangle(true),Triangle(est),Triangle(samp-emp),Triangle(emp1),Triangle(emp2),Triangle(rel-err(noisy)),Triangle(rel-error(samp-emp)),Triangle(rel-error(emp1)),Triangle(rel-error(emp2))\n");
		fclose(fp);
	}

    tri_num_ns_avg_re = tri_num_samp_emp_avg_re = tri_num_emp1_avg_re = tri_num_emp2_avg_re = 0.0;
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

			tri_num = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
					j = aitr->first;
					if (i >= j) continue;
					for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
						k = aitr2->first;
						if (j >= k) continue;
						if(a_mat[i].count(k) > 0){
							tri_num++;
						} 
					}
				}
			}
        }
		cout<<"Triangle: "<<tri_num<<endl;

		/************************ Calculate sub-graph counts ************************/

        CalcNLocShuffleTri(a_mat, deg, outfile, tri_num_ns, tri_num_samp_emp, tri_num_emp1, tri_num_emp2);

        tri_num_ns_re = fabs(tri_num_ns - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
        tri_num_samp_emp_re = fabs(tri_num_samp_emp - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
        tri_num_emp1_re = fabs(tri_num_emp1 - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
        tri_num_emp2_re = fabs(tri_num_emp2 - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);

        /**************************** Output the results ****************************/
        fp = FileOpen(outfile, "a");
        fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e\n", 
        (double)tri_num, tri_num_ns, tri_num_samp_emp, tri_num_emp1, tri_num_emp2, tri_num_ns_re, tri_num_samp_emp_re, tri_num_emp1_re, tri_num_emp2_re);
        fclose(fp);

		tri_num_ns_avg_re += tri_num_ns_re;
		tri_num_samp_emp_avg_re += tri_num_samp_emp_re;
		tri_num_emp1_avg_re += tri_num_emp1_re;
		tri_num_emp2_avg_re += tri_num_emp2_re;

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	tri_num_ns_avg_re /= (double)ItrNum;
	tri_num_samp_emp_avg_re /= (double)ItrNum;
	tri_num_emp1_avg_re /= (double)ItrNum;
	tri_num_emp2_avg_re /= (double)ItrNum;

	fp = FileOpen(outfile, "a");
	fprintf(fp, "function,AVG(rel-err(noisy)),AVG(rel-error(samp_emp)),AVG(rel-error(emp1)),AVG(rel-error(emp2))\n");
	fprintf(fp, "Triangles,%e,%e,%e,%e\n", tri_num_ns_avg_re, tri_num_samp_emp_avg_re, tri_num_emp1_avg_re, tri_num_emp2_avg_re);
	fclose(fp);

	// free
	free2D(node_order, ItrNum);
	free1D(deg);

	return 0;
}
