// #include "bits/stdc++.h"
#include <map>
long long int CalTrue4CliqueCount(map<int, int> *a_mat, long long int NodeNum){
    long long int cliq_4_num = 0;
    int j,k,l;
    for (int i=0; i<NodeNum; i++){
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
    return cliq_4_num;
}


// Emperical estimation is pending
void CalNILocFCliq(map<int, int> *a_mat, string outfile, double &cliq_4_num_ns, int emp){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<int, int>::iterator aitr3;
	long long cliq_4_num;
	// int *deg_ns;									// noisy degree
	// long long tot_edge_num_ns;
	// long long tri_num, st2_num, ed2_num, ed1_num, non_num;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int i, j, k, l;

	a_mat_ns = new map<int, int>[NodeNum];
	// malloc1D(&deg_ns, NodeNum); //Allocate memory for storing noisy degrees

	// Flip probability --> q
    q = 1.0 / (exp(Eps) + 1.0);

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

	if(emp == 1){

	}
	else{
		cliq_4_num_ns = (double)cliq_4_num;
	}

}