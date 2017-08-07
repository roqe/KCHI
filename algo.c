#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "algo.h"
#include "func.h"

void nuclearFamily (int *first, PERSON **population, int num_people) 
{
	int i, j, c, r;
	PERSON *N, *F, *M, *K;
	int *new_children, *new_relationship;
		
	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->father == -1 && N->mother == -1) continue;
		//update data of father
		F = population[N->father];
		++(F->num_relationship);
		if (F->num_relationship == F->cap_relationship) {
			F->cap_relationship *= 2;
			new_relationship = (int *) memset1D (F->cap_relationship, sizeof(int), 0);
			for (j = 1; j < F->num_relationship; ++j) new_relationship[j] = F->relationship[j];
			F->relationship = new_relationship;
		} F->relationship[F->num_relationship] = N->mother;
		++(F->num_children);
		if (F->num_children	== F->cap_children) {
			F->cap_children *= 2;
			new_children = (int *) memset1D (F->cap_children, sizeof(int), 0);
			for (j = 1; j < F->num_children; ++j) new_children[j] = F->children[j];
			F->children	= new_children;
		} F->children[F->num_children] = i;
		//update data of mother
		M = population[N->mother];
		++(M->num_relationship);
		if (M->num_relationship == M->cap_relationship) {
			M->cap_relationship *= 2;
			new_relationship = (int *) memset1D (M->cap_relationship, sizeof(int), 0);
			for (j = 1; j < M->num_relationship; ++j) new_relationship[j] = M->relationship[j];
			M->relationship = new_relationship;
		} M->relationship[M->num_relationship] = N->father;
		++(M->num_children);
		if (M->num_children	== M->cap_children) {
			M->cap_children *= 2;
			new_children = (int *) memset1D (M->cap_children, sizeof(int), 0);
			for (j = 1; j < M->num_children; ++j) new_children[j] = M->children[j];
			M->children	= new_children;
		} M->children[M->num_children] = i;		
	}

	 for (i = 1; i <= num_people; ++i) {
		N = population[i];	
		//index first child and assign local cross edge
		first[N->children[1]] = 1;
		if (N->num_children > 1 && N->sex == MALE) {
			for (c = 2; c <= N->num_children; ++c) {	
				K = population[N->children[c]];
				for (r = 1; r < c; ++r) {
					if (first[N->children[r]] == 1 && population[N->children[r]]->mother == K->mother) {
						K->to_father_edge->type = LOCAL_CROSS;
						break;
					}
				}
				if (K->to_father_edge->type == LOCAL_CROSS) continue;
				else first[N->children[c]] = 1;
			}
		}	
  	}
}

void initialP (GENOTYPE **g, BIT_VALUE **p, PERSON **population, int num_people, int num_locus) 
{
	int i, j;
	PERSON *N;

	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		for (j = 1; j <= num_locus; ++j) {
			if (N->father != -1 && N->mother != -1) {
				if (g[i][j] == G01 && (g[N->father][j] == G00 || g[N->mother][j] == G11)) 
					p[i][j] = 0;
				else if (g[i][j] == G01 && (g[N->father][j] == G11 || g[N->mother][j] == G00)) 
					p[i][j] = 1;
			} 
			if (g[i][j] == G00) p[i][j] = 0;
			else if (g[i][j] == G11) p[i][j] = 1;
			else if (g[i][j] == H01) p[i][j] = 0;
			else if (g[i][j] == H10) p[i][j] = 1;
		}
	}
}

void initialD (GENOTYPE **g, BIT_VALUE ***d, PERSON **population, int num_people, int num_locus) 
{
	int i, j;
	PERSON *N;
	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->father == -1 && N->mother == -1) continue;
		for (j = 1; j <= num_locus; ++j) {
			d[i][j][1] = 0;
			if (g[i][j] == G00 || g[i][j] == G11) d[i][j][0] = 0;
			else if (g[i][j] == G01 || g[i][j] == H01 || g[i][j] == H10) d[i][j][0] = 1;
			else d[i][j][0] = UNKNOWN;
		}
	}
}

void spanningTree (PERSON **population, int **global, int num_people) 
{
	int i, num_global_cycle = 0, *visit;
	visit = (int *) memset1D (num_people+1, sizeof(int), 0);
	PERSON *N;
	
	span (1, visit, population, global, num_people);

	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->father != -1 && N->mother != -1) {
			if (N->to_father_edge->type == GLOBAL_CROSS) {
				++num_global_cycle;
				global[num_global_cycle][0] = N->father;
				global[num_global_cycle][1] = i;
				//printf ("global: %2d %2d\n", i, N->father);
			}
			if (N->to_mother_edge->type == GLOBAL_CROSS) {
				++num_global_cycle;
				global[num_global_cycle][0] = N->mother;
				global[num_global_cycle][1] = i;
				//printf ("global: %2d %2d\n", i, N->mother);
			}
		}
	}
	
	global[0][0] = num_global_cycle; 
	//printf ("There are %d global cycles.\n", num_global_cycle);
	free (visit);
}

void span (int i, int *visit, PERSON **population, int **global, int num_people) 
{
	int c, child;
	PERSON *N, *K;
	N = population[i];

	if (N->father != -1 && N->mother != -1) {
		if (visit[N->mother] == 0 && N->to_mother_edge->type == GLOBAL_CROSS) {
			N->to_mother_edge->type = TREE;
			visit[N->mother] = 1;
			span (N->mother, visit, population, global, num_people);
		}
		if (visit[N->father] == 0 && N->to_father_edge->type == GLOBAL_CROSS) {
			N->to_father_edge->type = TREE;
			visit[N->father] = 1;
			span (N->father, visit, population, global, num_people);
		} 
	}

	for (c = 1; c <= N->num_children; ++c) {
		child = N->children[c];
		K = population[child];
		if (visit[child] == 0 && N->sex == MALE && K->to_father_edge->type == GLOBAL_CROSS) {
			K->to_father_edge->type = TREE;
			visit[child] = 1;
		} else if (visit[child] == 0 && N->sex == FEMALE && K->to_mother_edge->type == GLOBAL_CROSS) {
			K->to_mother_edge->type = TREE;
			visit[child] = 1;
		}
		span (child, visit, population, global, num_people);
	}
}

void treeEQ (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int *first, int **global,  int num_people, int num_locus) 
{
    int *visit, *sumD, i, j, root = 0;
    visit = (int *) memset1D (num_people+1, sizeof(int), 0);
	sumD = (int *) memset1D (num_people+1, sizeof(int), 0);
	PERSON *N;

	//generate tree constraint
	for (j = 1; j <= num_locus; ++j) {
		for (i = 1; i <= num_people; ++i) if (visit[i] == 0 && g[i][j] < 5) {
			if (p[i][j] == B0 || p[i][j] == B1) root = i;
			else root = 0;
			visit[i] = 1;
			DFS (i, j, root, g, d, p, population, global[0], sumD, visit, 't');
		}
		for (i = 1; i <= num_people; ++i) visit[i] = 0;
		root = 0;	
	}
	
	//map to constraint graph G*
	for (i = 1; i <= num_people; ++i) if (visit[i] == 0) {
		N = population[i];
		if (N->num_constraint == 0) continue;
		visit[i] = 1;
		N->root = i;
		N->W = 0;
		DFSW (i, i, population, visit);
	}

	assignFounder (population, first, num_people);

	free (visit);
	free (sumD);
}

void DFS (int i, int j, int root, GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, 
	PERSON **population,  int *global, int *sumD, int *visit, char mode) 
{
	int c, child, n, *temp;
	PERSON *N;
	N = population[i];
	temp = (int *) memset1D (N->num_children+3, sizeof(int), 0);	
	++(temp[0]);

	if (N->father != -1 && N->to_father_edge->type == TREE && visit[N->father] == 0 && g[N->father][j] != G00 && g[N->father][j] != G11 && g[N->father][j] < 5) {
		if (mode == 't') root = treeConstraint (i, N->father, j, temp[0], root, g, p, population, sumD, temp, d[i][j][1]);
		else if (mode == 'c') cycleConstraint (i, N->father, j, temp[0], global, p, population, sumD, temp, d[i][j][1], d);
		else if (mode == 'p') pathConstraint (i, N->father, j, temp[0], global, p, population, sumD, temp, d[i][j][1]);
		else if (mode == 'f') pathTransform (i, N->father, j, temp[0], global, p, population, sumD, temp, d[i][j][1], d);
	}

	if (N->mother != -1 && N->to_mother_edge->type == TREE && visit[N->mother] == 0 && g[N->mother][j] != G00 && g[N->mother][j] != G11 && g[N->mother][j] < 5) {
		if (mode == 't') root = treeConstraint (i, N->mother, j, temp[0], root, g, p, population, sumD, temp, d[i][j][0]);
		else if (mode == 'c') cycleConstraint (i, N->mother, j, temp[0], global, p, population, sumD, temp, d[i][j][0], d);
		else if (mode == 'p') pathConstraint (i, N->mother, j, temp[0], global, p, population, sumD, temp, d[i][j][0]);
		else if (mode == 'f') pathTransform (i, N->mother, j, temp[0], global, p, population, sumD, temp, d[i][j][0], d);
	}

	for (c = 1; c <= N->num_children; ++c) {
		child = N->children[c];
		if (N->sex == MALE && population[child]->to_father_edge->type == TREE && visit[child] == 0 && g[i][j] != G00 && g[i][j] != G11 && g[child][j] < 5) {
			if (mode == 't') root = treeConstraint (i, child, j, temp[0], root, g, p, population, sumD, temp, d[child][j][1]);
			else if (mode == 'c') cycleConstraint (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][1], d);
			else if (mode == 'p') pathConstraint (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][1]);
			else if (mode == 'f') pathTransform (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][1], d);
		} else if (N->sex == FEMALE && population[child]->to_mother_edge->type == TREE && visit[child] == 0 && g[i][j] != G00 && g[i][j] != G11 && g[child][j] < 5) {
			if (mode == 't') root = treeConstraint (i, child, j, temp[0], root, g, p, population, sumD, temp, d[child][j][0]);
			else if (mode == 'c') cycleConstraint (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][0], d);
			else if (mode == 'p') pathConstraint (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][0]);
			else if (mode == 'f') pathTransform (i, child, j, temp[0], global, p, population, sumD, temp, d[child][j][0], d);
		}
	}

	visit[i] = 1;
	for (n = 1; n < temp[0]; ++n) DFS (temp[n], j, root, g, d, p, population, global, sumD, visit, mode);
	free (temp);
}

int treeConstraint (int i, int r, int j, int t, int root, GENOTYPE **g, BIT_VALUE **p, 
	PERSON **population, int *sumD, int *temp, BIT_VALUE dd) 
{
	int c, leaf, **extd;
	PERSON *R, *L;

	if ((p[r][j] == B0 || p[r][j] == B1) && g[r][j] < 5) {
		if (root == 0) {
			root = r;		
			sumD[r] = sumD[i] ^ dd;
			temp[t] = r;
			++(temp[0]);
		} else {	
			leaf = r;
			sumD[leaf] = sumD[i] ^ dd;
			//printf ("tree constraint root: %d, leaf %d\n", root, leaf);
			R = population[root];
			if (R->num_constraint == R->cap_constraint) {
				extd = (int **) memset2D (R->cap_constraint * 2, 2, sizeof(int), -1);
				for (c = 0; c < R->num_constraint; ++c) {
					extd[c][0] = R->constraint[c][0];
					extd[c][1] = R->constraint[c][1];
				}
				free2D(R->constraint, R->cap_constraint);
				R->constraint = extd;
				R->cap_constraint *= 2;
			}
			R->constraint[R->num_constraint][0] = leaf;
			R->constraint[R->num_constraint][1] = p[root][j] ^ p[leaf][j] ^ sumD[root] ^ sumD[leaf];
			++(R->num_constraint);
			L = population[leaf];
			if (L->num_constraint == L->cap_constraint) {
				extd = (int **) memset2D (L->cap_constraint * 2, 2, sizeof(int), -1);
				for (c = 0; c < L->num_constraint; ++c) {
					extd[c][0] = L->constraint[c][0];
					extd[c][1] = L->constraint[c][1];
				}
				free2D(L->constraint, L->cap_constraint);
				L->constraint = extd;
				L->cap_constraint *= 2;
			}
			L->constraint[L->num_constraint][0] = root;
			L->constraint[L->num_constraint][1] = p[root][j] ^ p[leaf][j] ^ sumD[root] ^ sumD[leaf];
			++(L->num_constraint);
		}
	} else {
		sumD[r] = sumD[i] ^ dd;
		temp[t] = r;
		++(temp[0]);
	}

	return root;
}

void DFSW (int i, int r, PERSON **population, int *visit) 
{
	int c, con, n, t = 0, *temp;
	PERSON *L, *N;
	N = population[i];
	temp = (int *) memset1D (N->num_constraint, sizeof(int), 0);
	
	for (c = 0; c < N->num_constraint; ++c) {
		con = N->constraint[c][0];
		L = population[con];
		if (visit[con] == 0) {
			visit[con] = 1;
			L->root = r;
			L->W = N->W ^ N->constraint[c][1];
			temp[t] = con;
			++t;
		} else if (L->W != (N->W ^ N->constraint[c][1])) {
	        fprintf (stderr, "tree constraint inconsistent\n");
		   	exit(EXIT_FAILURE);
		}
	}

	for (n = 0; n < t; ++n) DFSW(temp[n], r, population, visit);
	free (temp);
}

void assignFounder (PERSON **population, int *first, int num_people) 
{
	int i, c;
	PERSON *N, *K;

	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->father == -1 && N->mother == -1 && N->root == -1) {
			for (c = 1; c <= N->num_children; ++c) {
				K = population[N->children[c]];
				if (first[N->children[c]] == 1 && K->W != UNKNOWN) {
					if (N->sex == MALE && K->to_father_edge->h == UNKNOWN && K->to_father_edge->type == TREE) {
						N->W = K->W; 
						N->root = K->root;
						K->to_father_edge->h = B0;
					} else if (N->sex == FEMALE && K->to_mother_edge->h == UNKNOWN && K->to_mother_edge->type == TREE) {
						N->W = K->W; 
						N->root = K->root;
						K->to_mother_edge->h = B0;
					}
					break;
				}
			}
		}
	}
	//solve h of tree edge
	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->W != UNKNOWN && N->father != -1 && population[N->father]->W != UNKNOWN && N->to_father_edge->type == TREE) 
			if (N->root == population[N->father]->root) N->to_father_edge->h = N->W ^ population[N->father]->W;			
		if (N->W != UNKNOWN	&& N->mother != -1 && population[N->mother]->W != UNKNOWN && N->to_mother_edge->type == TREE) 
			if (N->root == population[N->mother]->root) N->to_mother_edge->h = N->W ^ population[N->mother]->W;
	}
}

void localEQ (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int *first, int **local, int num_people, int num_locus) 
{
	int i, j, c, r, b, k;
	BIT_VALUE bc = UNKNOWN, bp = UNKNOWN;
	PERSON *F, *K, *B;

	for (i = 1; i <= num_people; ++i) {
		F = population[i];
		if (F->sex == MALE) {
			for (c = 1; c <= F->num_children; ++c) {
				//first child with no local constraint
				if (first[F->children[c]] == 1) continue;
				k = F->children[c];
				K = population[k];
				//find the first child
				for (r = 1; r < c; ++r) {
					if (first[F->children[r]] == 1 && population[F->children[r]]->mother == K->mother) {
						b = F->children[r];
						B = population[b];
					}
				}
				//local cycle constraint
				if (B->to_father_edge->type == TREE && B->to_mother_edge->type == TREE && K->to_mother_edge->type == TREE) {
					for (j = 1; j <= num_locus; ++j) {
						//exclude missing
						if (g[i][j] > 4 || g[k][j] > 4 || g[b][j] > 4) continue;
						//dad is homozygous >> tree constraint
						if (g[K->father][j] == G00 || g[K->father][j] == G11) continue;
						//both are heterozygous >> local cycle constraint
						if (g[K->mother][j] != G00 && g[K->mother][j] != G11 && g[K->mother][j] < 5) {
							bc = d[b][j][0] ^ d[k][j][0];
							if (K->W != UNKNOWN && F->W != UNKNOWN && K->root == F->root && K->to_father_edge->h == UNKNOWN) { 
								K->to_father_edge->h = (K->W ^ F->W ^ bc);
							} else if (K->W != UNKNOWN && F->W != UNKNOWN && K->root == F->root && K->to_father_edge->h != (K->W ^ F->W ^ bc)) {
						        fprintf (stderr, "local cycle constraint inconsistent\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				} 
				//local path constraint
				if (B->to_father_edge->type == TREE) {
					for (j = 1; j <= num_locus; ++j) {
						//exclude missing
						if (g[i][j] > 4 || g[k][j] > 4 || g[b][j] > 4) continue;
						//dad is homozygous >> tree constraint
						if (g[K->father][j] == G00 || g[K->father][j] == G11) continue;
						//mom is heterozygous >> local path constraint
						if (p[b][j] != UNKNOWN && p[k][j] != UNKNOWN) {
							bp = p[b][j] ^ p[k][j];
							if (bc == UNKNOWN && B->W != UNKNOWN && K->W != UNKNOWN && B->root == K->root) 
								bc = bp ^ B->W ^ K->W;
							//transform from path to tree
							if (bc != UNKNOWN) {
								if (K->W == UNKNOWN && B->W != UNKNOWN) {
									K->W = B->W ^ bc ^ bp;
									K->root = B->root;
									if (K->to_mother_edge->type == TREE && population[K->mother]->W != UNKNOWN 
										&& K->root == population[K->mother]->root)
										K->to_mother_edge->h = K->W ^ population[K->mother]->W;
								} else if (K->W != UNKNOWN && B->W == UNKNOWN) {
									B->W = K->W ^ bc ^ bp;
									B->root = K->root;
									if (B->to_mother_edge->type == TREE && population[B->mother]->W != UNKNOWN 
										&& B->root == population[B->mother]->root)
										B->to_mother_edge->h = B->W ^ population[B->mother]->W;
									if (F->W != UNKNOWN && F->root == B->root) B->to_father_edge->h = (B->W ^ F->W);
								} else if (K->W == UNKNOWN && B->W == UNKNOWN) {
									B->W = population[B->father]->W;
									B->to_father_edge->h = FREE;
									K->W = B->W ^ bc ^ bp;
								} else if (K->root != B->root) {
									++(local[0][0]);
									local[local[0][0]][0] = k;
									local[local[0][0]][1] = b;
									K->to_father_edge->d = bc ^ bp;
								} else if (K->root == B->root && K->W != (B->W ^ bc ^ bp)) {
									fprintf (stderr, "local path constraint inconsistent\n");
									exit(EXIT_FAILURE);
								}
								if (K->W != UNKNOWN && F->W != UNKNOWN && K->root == F->root && K->to_father_edge->h == UNKNOWN) { 
									K->to_father_edge->h = K->W ^ F->W ^ bc;
								} else if (K->W != UNKNOWN && F->W != UNKNOWN && K->root == F->root 
									&& K->to_father_edge->h != (K->W ^ F->W ^ bc)) {
					        			fprintf (stderr, "local cycle constraint inconsistent\n");
									exit(EXIT_FAILURE);
								}					
							} else if (B->W != UNKNOWN && F->W != UNKNOWN && K->to_father_edge->h == UNKNOWN) { 
								K->to_father_edge->h = B->W ^ F->W ^ bp;
								if (B->root != F->root) {
									++(local[0][0]);
									local[local[0][0]][0] = b;
									local[local[0][0]][1] = i;
								}
							} else if (B->W != UNKNOWN && F->W != UNKNOWN && K->to_father_edge->h != (B->W ^ F->W ^ bp)) {
						        	fprintf (stderr, "local path constraint inconsistent\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				}
				bc = UNKNOWN;
				bp = UNKNOWN;
			}
		}
	}
	assignFounder (population, first, num_people);
}

void globalCEQ (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int **global, int num_people, int num_locus) 
{
	int j, q;
	int *sumD, *visit;

	for (q = 1; q <= global[0][0]; ++q) {
		//global cycle constraint
		for (j = 1; j <= num_locus; ++j) if (p[global[q][1]][j] == UNKNOWN && p[global[q][0]][j] == UNKNOWN) {
			sumD = (int *) memset1D (num_people+1, sizeof(int), 0);
			visit = (int *) memset1D (num_people+1, sizeof(int), 0);
			visit[global[q][1]] = 1;
			DFS (global[q][1], j, global[q][0], g, d, p, population, global[q], sumD, visit, 'c');
			free (sumD);
			free (visit);
		}
	}
}

void cycleConstraint (int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
	int *sumD, int *temp, BIT_VALUE dd, BIT_VALUE ***d) 
{
	BIT_VALUE dc = UNKNOWN;
	BIT_VALUE h = UNKNOWN;
	PERSON *G0, *G1;
	G0 = population[global[0]];
	G1 = population[global[1]];

	if (p[r][j] != B0 && p[r][j] != B1) {
		sumD[r] = sumD[i] ^ dd;
		temp[t] = r;
		++(temp[0]);
		if (global[0] == r) {
			if (r == G1->mother) {
				dc = sumD[r] ^ d[global[1]][j][0];
				G1->to_mother_edge->d = dc;  
				if (G1->W != UNKNOWN && G0->W != UNKNOWN && G0->root == G1->root) {
					h = dc ^ G1->W ^ G0->W;
					if (h != FREE && (G1->to_mother_edge->h == UNKNOWN || G1->to_mother_edge->h == FREE))
						G1->to_mother_edge->h = h;
					else if (h == FREE)
						printf ("global cycle constraint is FREE because of missing data\n");
					else if (G1->to_mother_edge->h != h)  
						printf ("global cycle constraint inconsistent 534\n");
				}
			} else if (r == G1->father) {
				dc = sumD[r] ^ d[global[1]][j][1];
				G1->to_father_edge->d = dc;  
				if (G1->W != UNKNOWN && G0->W != UNKNOWN && G0->root == G1->root) {
					h = dc ^ G1->W ^ G0->W;
					if (h != FREE && (G1->to_father_edge->h == UNKNOWN || G1->to_father_edge->h == FREE))
						G1->to_father_edge->h = h;
					else if (h == FREE)
						printf ("global cycle constraint is FREE because of missing data\n");
					else if (G1->to_father_edge->h != h) 
						printf ("global cycle constraint inconsistent 546\n");
				}
			}
		}
	} 
}

void globalPEQ (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int *first, int **global, int **odd, int num_people, int num_locus, char mode) 
{
	int j, q;
	int *sumD, *visit, **oddpc;
	PERSON *G1, *G0;
	BIT_VALUE path = UNKNOWN;

	for (q = 1; q <= global[0][0]; ++q) {
		//global path constraint 
		G0 = population[global[q][0]];
		G1 = population[global[q][1]];

		if (G0->sex == MALE && G1->to_father_edge->h != UNKNOWN) continue;
		else if (G0->sex == FEMALE && G1->to_mother_edge->h != UNKNOWN) continue;
		oddpc = (int **) memset2D (3, num_locus, sizeof(int), 0);

		for (j = 1; j <= num_locus; ++j) {
			sumD = (int *) memset1D (num_people+1, sizeof(int), 0);
			visit = (int *) memset1D (num_people+1, sizeof(int), 0);
			visit[global[q][1]] = 1;
			visit[global[q][0]] = 1;
			
			if (p[global[q][0]][j] != UNKNOWN && p[global[q][1]][j] != UNKNOWN && g[global[q][0]][j] < 5 && g[global[q][1]][j] < 5 
				&& g[global[q][0]][j] != G00 && g[global[q][0]][j] != G11) {
				if (G0->sex == MALE) {
					G1->to_father_edge->h = p[global[q][0]][j] ^ p[global[q][1]][j] ^ d[global[q][1]][j][1];
					break;
				} else if (G0->sex == FEMALE) {
					G1->to_mother_edge->h = p[global[q][0]][j] ^ p[global[q][1]][j] ^ d[global[q][1]][j][0];
					break;
				}
			} else if (p[global[q][1]][j] != UNKNOWN && p[global[q][0]][j] == UNKNOWN && g[global[q][0]][j] < 5) { 
				DFS (global[q][0], j, global[q][0], g, d, p, population, global[q], sumD, visit, mode);
				if (G0->W != UNKNOWN && G0->root == population[global[q][2]]->root) {
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN) {
						G1->to_father_edge->h = G0->W ^ G1->to_father_edge->d ^ p[global[q][1]][j] ^ d[global[q][1]][j][1];
						break;
					} else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN) {
						G1->to_mother_edge->h = G0->W ^ G1->to_mother_edge->d ^ p[global[q][1]][j] ^ d[global[q][1]][j][0];
						break;
					}
				} else if (G1->root != population[global[q][2]]->root) {	
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = G1->to_father_edge->d ^ p[global[q][1]][j] ^ d[global[q][1]][j][1];
						oddpc[1][oddpc[0][0]] = global[q][1];
						oddpc[2][oddpc[0][0]] = global[q][2];
					} else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = G1->to_mother_edge->d ^ p[global[q][1]][j] ^ d[global[q][1]][j][0];
						oddpc[1][oddpc[0][0]] = global[q][1];
						oddpc[2][oddpc[0][0]] = global[q][2];
					}
				}
			} else if (p[global[q][0]][j] != UNKNOWN && g[global[q][0]][j] == G01 && p[global[q][1]][j] == UNKNOWN && g[global[q][1]][j] < 5) {
				DFS (global[q][1], j, global[q][1], g, d, p, population, global[q], sumD, visit, mode);
				if (G1->W != UNKNOWN && G1->root == population[global[q][2]]->root) {
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN) {
						G1->to_father_edge->h = G1->W ^ G1->to_father_edge->d ^ p[global[q][0]][j] ^ d[global[q][1]][j][1];
						break;
					} else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN) {
						G1->to_mother_edge->h = G1->W ^ G1->to_mother_edge->d ^ p[global[q][0]][j] ^ d[global[q][1]][j][0];
						break;	
					}
				} else if (G0->root != population[global[q][2]]->root) {	
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = G1->to_father_edge->d ^ p[global[q][0]][j] ^ d[global[q][1]][j][1];
						oddpc[1][oddpc[0][0]] = global[q][0];
						oddpc[2][oddpc[0][0]] = global[q][2];
					} else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = G1->to_mother_edge->d ^ p[global[q][0]][j] ^ d[global[q][1]][j][0];
						oddpc[1][oddpc[0][0]] = global[q][0];
						oddpc[2][oddpc[0][0]] = global[q][2];
					}
				}
			} else if (p[global[q][0]][j] == UNKNOWN && p[global[q][1]][j] == UNKNOWN && g[global[q][0]][j] < 5 && g[global[q][1]][j] < 5) {
				DFS (global[q][0], j, global[q][0], g, d, p, population, global[q], sumD, visit, mode);
				oddpc[1][0] = global[q][2];
				if (G0->sex == MALE) path = G1->to_father_edge->d;
				else if (G0->sex == FEMALE) path = G1->to_mother_edge->d;
				DFS (global[q][1], j, global[q][1], g, d, p, population, global[q], sumD, visit, mode);
				oddpc[2][0] = global[q][2];
				if (G0->root == population[oddpc[1][0]]->root && G1->root == population[oddpc[2][0]]->root 
					&& G0->W != UNKNOWN && G1->W != UNKNOWN) {
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN && path != UNKNOWN) {
						G1->to_father_edge->h = path ^ G1->to_father_edge->d ^ G0->W ^ G1->W ^ d[global[q][1]][j][1];
						break;
					} else if (G0->sex == FEMALE && G1->to_mother_edge->h != UNKNOWN && path != UNKNOWN) {
						G1->to_mother_edge->h = path ^ G1->to_mother_edge->d ^ G0->W ^ G1->W ^ d[global[q][1]][j][0];
						break;
					}
				} else if (population[oddpc[1][0]]->root != population[oddpc[2][0]]->root) {	
					//odd paths may be combined into cycle her
					//oddpc[0][0]: number of odd paths, oddpc[1][0]: current end 1, oddpc[2][0]: current end 2
					//oddpc[0][x]: the b+W+W of the path x, oddpc[1][x]: end 1 of path x, oddpc[2][x]: end 2 of path x
					if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = path ^ G1->to_father_edge->d ^ d[global[q][1]][j][1];
						oddpc[1][oddpc[0][0]] = oddpc[1][0];
						oddpc[2][oddpc[0][0]] = oddpc[2][0];
					} else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN) {
						++(oddpc[0][0]);
						oddpc[0][oddpc[0][0]] = path ^ G1->to_mother_edge->d ^ d[global[q][1]][j][0];
						oddpc[1][oddpc[0][0]] = oddpc[1][0];
						oddpc[2][oddpc[0][0]] = oddpc[2][0];
					}
				}
				if (G0->sex == MALE && G1->to_father_edge->h != UNKNOWN) break;
				else if (G0->sex == FEMALE && G1->to_mother_edge->h != UNKNOWN) break;
			}
			if (G1->to_father_edge->h == UNKNOWN) G1->to_father_edge->d = UNKNOWN;
			if (G1->to_mother_edge->h == UNKNOWN) G1->to_mother_edge->d = UNKNOWN;	
		}
		if (oddpc[0][0] > 2) oddpathstocycle (population, oddpc, odd, num_people);
		free (sumD);
		free (visit);
		free2D (oddpc, 3);
	}
	assignFounder (population, first, num_people);
}

void oddpathstocycle (PERSON **population, int **oddpc, int **odd, int num_people) 
{
	int i, R1, R2, c1, c2, *visit, c, t;
	visit = (int *) memset1D (num_people+1, sizeof(int), 0);
	BIT_VALUE o, b;

	for (i = 1; i <= oddpc[0][0]; ++i) {
		R1 = population[oddpc[1][i]]->root;
		R2 = population[oddpc[2][i]]->root;
		if (population[R1]->num_constraint != 10000) {
			population[R1]->constraint = (int **) memset2D (oddpc[0][0]+1, 2, sizeof(int), 0); 
			population[R1]->num_constraint = 10000;
		}
		if (population[R2]->num_constraint != 10000) {
			population[R2]->constraint = (int **) memset2D (oddpc[0][0]+1, 2, sizeof(int), 0); 
			population[R2]->num_constraint = 10000;
		}
		++(population[R1]->constraint[0][0]);
		++(population[R2]->constraint[0][0]);
		c1 = population[R1]->constraint[0][0];
		c2 = population[R2]->constraint[0][0];
		population[R1]->constraint[c1][0] = R2;
		population[R1]->constraint[c1][1] = oddpc[0][i];
		population[R2]->constraint[c2][0] = R1;
		population[R2]->constraint[c2][1] = oddpc[0][i];
	}

	for (i = 1; i <= oddpc[0][0]; ++i) {
		R1 = population[oddpc[1][i]]->root;
		o = B0;
		b = B0;
		if (visit[R1] == 0) {
			visit[R1] = 1;
			c = DFSO (population, R1, R1, o, b, oddpc, visit);
		}
	}
	if (c != -1) for (i = 1; i <= oddpc[0][0]; ++i) {
		++(odd[0][0]);
		t = odd[0][0];
		odd[0][t] = oddpc[0][i] ^ c;
		odd[1][t] = population[oddpc[1][i]]->root;	
		odd[2][t] = population[oddpc[2][i]]->root;
	}
}

int DFSO (PERSON **population, int R1, int RP, BIT_VALUE o, BIT_VALUE b, int **oddpc, int *visit) 
{
	int n, t = 1, *temp, c;
	PERSON *N;

	N = population[R1];
	temp = (int *) memset1D (oddpc[0][0], sizeof(int), 0);
	visit[R1] = 1;

	for (n = 1; n <= N->constraint[0][0]; ++n) {
		if (visit[N->constraint[n][0]] == 0) {
			temp[t] = N->constraint[n][0];
			++t;
			o = o ^ B1;
			b = b ^ N->constraint[n][1];
		} else if (N->constraint[n][0] != RP && o == B0) {
			o = o ^ B1;
			b = b ^ N->constraint[n][1];
			return b;	
		}
	}

	for (n = 1; n < t; ++n) {
		c = DFSO (population, temp[n], R1, o, b, oddpc, visit);
		if (c != -1) return c;
	}

	free(temp);
	return -1;
}

void pathConstraint (int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
	int *sumD, int *temp, BIT_VALUE dd) 
{
	int leaf;
	PERSON *G0, *G1, *L;
	G0 = population[global[0]];
	G1 = population[global[1]];

	if (p[r][j] == B0 || p[r][j] == B1) {
		leaf = r;
		sumD[leaf] = sumD[i] ^ dd;
		L = population[leaf];
		if (L->W == UNKNOWN) { 
			L->W = 0;
			L->root = leaf;
		}
		if (G0->sex == MALE && G1->to_father_edge->h == UNKNOWN) {
			G1->to_father_edge->d = sumD[leaf] ^ p[r][j] ^ L->W;
		} else if (G0->sex == FEMALE && G1->to_mother_edge->h == UNKNOWN) {
			G1->to_mother_edge->d = sumD[leaf] ^ p[r][j] ^ L->W;
		}
		global[2] = leaf;
	} else {
		sumD[r] = sumD[i] ^ dd;
		temp[t] = r;
		++(temp[0]);
	}
}

void pathTransform (int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
	int *sumD, int *temp, BIT_VALUE dd, BIT_VALUE ***d) 
{
	int leaf, c, child;
	PERSON *G0, *G1, *L, *K;
	G0 = population[global[0]];
	G1 = population[global[1]];
	BIT_VALUE leafW = UNKNOWN;

	if (p[r][j] == B0 || p[r][j] == B1) {
		leaf = r;
		sumD[leaf] = sumD[i] ^ dd;
		L = population[leaf];
		if (p[global[0]][j] != UNKNOWN && G0->W != UNKNOWN) {
			if (G0->sex == MALE)
				leafW = G0->W ^ p[global[0]][j] ^ d[global[1]][j][1] ^ sumD[leaf] ^ p[leaf][j] ^ G1->to_father_edge->h;
			else if (G0->sex == FEMALE)  
				leafW = G0->W ^ p[global[0]][j] ^ d[global[1]][j][0] ^ sumD[leaf] ^ p[leaf][j] ^ G1->to_mother_edge->h;
		} else if (p[global[1]][j] != UNKNOWN && G1->W != UNKNOWN) {
			if (G0->sex == MALE)
				leafW = G1->W ^ p[global[1]][j] ^ d[global[1]][j][1] ^ sumD[leaf] ^ p[leaf][j] ^ G1->to_father_edge->h;
			else if (G0->sex == FEMALE)  
				leafW = G1->W ^ p[global[1]][j] ^ d[global[1]][j][0] ^ sumD[leaf] ^ p[leaf][j] ^ G1->to_mother_edge->h;
		} else if (p[global[0]][j] == UNKNOWN && p[global[1]][j] == UNKNOWN) {
			if (G0->sex == MALE && G1->to_father_edge->d != UNKNOWN)
				leafW = G1->to_father_edge->d ^ d[global[1]][j][1] ^ sumD[leaf] ^ G1->to_father_edge->h;
			else if (G0->sex == FEMALE && G1->to_mother_edge->d != UNKNOWN)
				leafW = G1->to_mother_edge->d ^ d[global[1]][j][0] ^ sumD[leaf] ^ G1->to_mother_edge->h;
		}
		if (L->to_father_edge->type == TREE && L->to_father_edge->h == UNKNOWN && population[L->father]->W != UNKNOWN && leafW != UNKNOWN)  
			L->to_father_edge->h = leafW ^ population[L->father]->W;
		if (L->to_mother_edge->type == TREE && L->to_mother_edge->h == UNKNOWN && population[L->mother]->W != UNKNOWN && leafW != UNKNOWN) 
			L->to_mother_edge->h = leafW ^ population[L->mother]->W;

		for (c = 1; c <= L->num_children; ++c) if (leafW != UNKNOWN) {
			child = L->children[c]; 
			K = population[child];
			if (K->to_father_edge->type == TREE && K->to_father_edge->h == UNKNOWN && L->sex == MALE && K->W != UNKNOWN) 
				K->to_father_edge->h = leafW ^ population[child]->W;
			else if (K->to_mother_edge->type == TREE && K->to_mother_edge->h == UNKNOWN && L->sex == FEMALE && K->W != UNKNOWN) 
				K->to_mother_edge->h = leafW ^ population[child]->W;
		}
	} else {
		sumD[r] = sumD[i] ^ dd;
		temp[t] = r;
		++(temp[0]);
	}
}

void linkTree (GENOTYPE **g, BIT_VALUE **p, PERSON **population, int *first, int **local, int **global, 
	int **odd, int num_people, int num_locus) 
{
	int i, j, b, k, h, c, r, m;
	PERSON *N, *F, *M, *K, *B; 

	for (i = 1; i <= num_people; ++i) {
		N = population[i];
		if (N->father != -1) {
			F = population[N->father];
			if (N->to_father_edge->type == TREE && N->to_father_edge->h == UNKNOWN && N->root != -1 && F->root != -1) {
				h = isLink (population, local, global, odd, N->root, F->root); 
				if (h != -1 && h != UNKNOWN && h != FREE) N->to_father_edge->h = N->W ^ F->W ^ h;
			}
		}	
		if (N->mother != -1) {
			M = population[N->mother];
			if (N->to_mother_edge->type == TREE && N->to_mother_edge->h == UNKNOWN && N->root != -1 && M->root != -1) {
				h = isLink (population, local, global, odd, N->root, M->root);
				if (h != -1 && h != UNKNOWN && h != FREE) N->to_mother_edge->h = N->W ^ M->W ^ h;
			}
		}
		if (N->sex == MALE) {
			for (c = 1; c <= N->num_children; ++c) {
				//first child with no local constraint
				if (first[N->children[c]] == 1) continue;
				k = N->children[c];
				K = population[k];
				//find the first child
				for (r = 1; r < c; ++r) {
					if (first[N->children[r]] == 1 && population[N->children[r]]->mother == K->mother) {
						m = K->mother;
						M = population[m];
						b = N->children[r];
						B = population[b];
					}
				}
				for (j = 1; j <= num_locus; ++j) {
					if ((g[i][j] == G01 || g[i][j] == H01 || g[i][j] == H10) && p[k][j] != UNKNOWN && p[b][j] != UNKNOWN 
							&& B->to_father_edge->h != UNKNOWN && K->to_father_edge->h == UNKNOWN) {
						K->to_father_edge->h = B->to_father_edge->h ^ p[k][j] ^ p[b][j];
						break;
					}
				}
			}
		}
	}
}

int isLink (PERSON **population, int **local, int **global, int **odd, int n, int m) {
	int i;
	PERSON *L0, *L1, *G0, *G1;

	for (i = 1; i <= local[0][0]; i++) {
		L0 = population[local[i][0]];
		L1 = population[local[i][1]];
		if ((L0->root == n && L1->root == m) || (L0->root == m && L1->root == n)) 
			return L0->to_father_edge->h;
	}
	for (i = 1; i <= global[0][0]; i++) {
		G0 = population[global[i][0]];
		G1 = population[global[i][1]];
		if (G0->root == G1->root) continue;
		if ((G0->root == n && G1->root == m) || (G0->root == m && G1->root == n)) {
			if (G0->sex == MALE) return G1->to_father_edge->h;
			else return G1->to_mother_edge->h;
		}
	}
	for (i = 1; i <= odd[0][0]; i++) {
		if ((odd[1][i] == n && odd[2][i] == m) || (odd[1][i] == m && odd[2][i] == n)) 
			return odd[0][i];
	}
	return -1;
}

void calculateP (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int num_people, int num_locus) 
{
	int i, j, *visit;
	visit = (int *) memset1D (num_people+1, sizeof(int), 0);
	
	for (j = 1; j <= num_locus; ++j) {
		for (i = 1; i <= num_people; ++i) {
			if (p[i][j] != UNKNOWN) {
				if (g[i][j] == G01 && p[i][j] == B0) g[i][j] = H01;
				if (g[i][j] == G01 && p[i][j] == B1) g[i][j] = H10;
				visit[i] = 1;
				solveP (g, d, p, population, i, j, visit);
			}
		}
		for (i = 1; i <= num_people; ++i) visit[i] = 0;
	}
	free (visit);
}

void solveP (GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
	int i, int j, int *visit) 
{
	int c, child, t = 1, n;
	int *temp;
	BIT_VALUE w;
	PERSON *N, *K;

	N = population[i];
	temp = (int *) memset1D ((3 + N->num_children), sizeof(int), 0);

	if (N->father != -1 && visit[N->father] == 0) {
		visit[N->father] = 1;
		if (p[N->father][j] == UNKNOWN && p[i][j] != UNKNOWN && g[N->father][j] < 5) {
			if (g[N->father][j] == G00 || g[N->father][j] == G11) w = B0;	
			else w = B1;
			p[N->father][j] = p[i][j] ^ (N->to_father_edge->h * w) ^ d[i][j][1]; 
			if (g[N->father][j] == G01 && p[N->father][j] == B0) g[N->father][j] = H01;
			if (g[N->father][j] == G01 && p[N->father][j] == B1) g[N->father][j] = H10;
		}
		if (p[N->father][j] != UNKNOWN) {
			temp[t] = N->father;
			++t;
		} 
	}
	if (N->mother != -1 && visit[N->mother] == 0) {
		visit[N->mother] = 1;
		if (p[N->mother][j] == UNKNOWN && p[i][j] != UNKNOWN && g[N->mother][j] < 5) {
			if (g[N->mother][j] == G00 || g[N->mother][j] == G11) w = B0;	
			else w = B1;
			p[N->mother][j] = p[i][j] ^ (N->to_mother_edge->h* w) ^ d[i][j][0];
			if (g[N->mother][j] == G01 && p[N->mother][j] == B0) g[N->mother][j] = H01;
			if (g[N->mother][j] == G01 && p[N->mother][j] == B1) g[N->mother][j] = H10;
		}
		if (p[N->mother][j] != UNKNOWN) {
			temp[t] = N->mother;
			++t;
		} 
	} 

	if (p[i][j] != UNKNOWN) {
		for (c = 1; c <= N->num_children; ++c) {
			child = N->children[c];
			K = population[child];
			if (visit[child] == 0) {
				visit[child] = 1;
				if (p[child][j] == UNKNOWN) {
					if (g[child][j] < 5) {
						if (g[i][j] == G00 || g[i][j] == G11) w = B0;	
						else w = B1;
						if (N->sex == MALE) p[child][j] = p[i][j] ^ (K->to_father_edge->h * w) ^ d[child][j][1];
						else p[child][j] = p[i][j] ^ (K->to_mother_edge->h * w) ^ d[child][j][0];
						if (g[child][j] == G01 && p[child][j] == B0) g[child][j] = H01;
						if (g[child][j] == G01 && p[child][j] == B1) g[child][j] = H10;
					}	 
				} 
				if (p[child][j] != UNKNOWN) {
					temp[t] = child;
					++t;
				} 
			}
		}
	}

	for (n = 1; n < t; ++n) solveP(g, d, p, population, temp[n], j, visit);
	free (temp);
}

void statisitics (GENOTYPE **g, BIT_VALUE **p, PERSON **population, BIT_VALUE **realh, BIT_VALUE **realp, 
	int *first, int num_missing, int num_people, int num_locus) 
{
	int i, j, *visit, *error, nosol, nofix, ptfix;
	char hap[4];
	FILE *res, *sum;
	visit = (int *) memset1D (num_people+1, sizeof(int), 0);
	
	res = fopen("res", "a");
	fprintf (res, "num_people: %d, num_locus: %d\n\n", num_people, num_locus);
	sum = fopen("sum", "a");
	fprintf (sum, "num_people: %d, num_locus: %d, num_missing: %d\n", num_people, num_locus, num_missing);
	
	assignFounder (population, first, num_people);
	error = (int *) memset1D (2, sizeof(int), 0); 	
	//error[0] = h error, error[1] = p error, 
	//confirm h 
	visit[1] = 1;
	DFSH (population, 1, realh, visit, error);
	for (j = 1; j <= num_locus; ++j) {
		for (i = 1; i <= num_people; ++i) visit[i] = 0;
		visit[1] = 1;
		if (p[1][j] != realp[1][j]) ++(error[1]);
		DFSP (population, 1, j, p, realp, visit, error);
	}
	//results print to file
	nosol = 0;
	nofix = 0;
	ptfix = 0;
	for (i = 1; i <= num_people; ++i) {
		for (j = 1; j <= num_locus; ++j) {
			if (g[i][j] == G00) strcpy (hap, "00");
			else if (g[i][j] == G11) strcpy (hap, "11");
			else if (g[i][j] == G01) strcpy (hap, "ns");
			else if (g[i][j] == H01) strcpy (hap, "01");
			else if (g[i][j] == H10) strcpy (hap, "10");
			else if (g[i][j] == G1x) strcpy (hap, "1?");
			else if (g[i][j] == G0x) strcpy (hap, "0?");
			else if (g[i][j] == Gxx) strcpy (hap, "xx");
			else if (g[i][j] == H0x) strcpy (hap, "0x");
			else if (g[i][j] == H1x) strcpy (hap, "1x");
			else if (g[i][j] == Hx0) strcpy (hap, "x0");
			else if (g[i][j] == Hx1) strcpy (hap, "x1");
			else strcpy (hap, "er");
			fprintf (res, "%s ", hap);
			if (g[i][j] == G01) { 
				++nosol; //no hap cuz the link breaks by missing data
			}
			if (g[i][j] == Gxx) ++nofix; //not fix
			else if	(g[i][j] != Gxx && g[i][j] > 4 && g[i][j] < 12) ++ptfix; //partial fix
		}
		fprintf (res, "\n");
	}
	fprintf (res, "\n");
	if (num_missing != 0) {
		fprintf (sum, "\ttotally recover:\t%f%%\n", (double)(num_missing-nofix-ptfix)/num_missing*100.0);
		fprintf (sum, "\tpartially recover:\t%f%%\n", (double)(ptfix)/num_missing*100.0);
		fprintf (sum, "\tmissing recover:\t%f%%\n", (double)(num_missing-nofix)/num_missing*100.0);
		fprintf (sum,"\n");	
	} 
	fprintf(sum, "\ttotal h mistakes: %d\n", error[0]);
	fprintf(sum, "\ttotal p mistakes: %d\n", error[1]);
	if (num_missing != 0) {
		fprintf(sum, "\t\t G01 break by missing: %d\n", nosol);
		fprintf(sum, "\t\t Gxx missing no fix: %d\n", nofix);
		fprintf(sum, "\t\t missing partial fix: %d\n", ptfix);
	}
	fprintf(sum, "\tcorrectly phased markers: %f%%\n", 
		(double)(num_people*num_locus-error[1])/(num_people*num_locus)*100.0);
	fclose (sum);
	fclose (res);

	free (error);
	free (visit);
}

void DFSH (PERSON **population, int i, BIT_VALUE **realh, int *visit, int *error) 
{
	int c, child, t = 1, n;
	int *temp;
	PERSON *N, *K;

	N = population[i];
	temp = (int *) memset1D (2+N->num_children, sizeof(int), 0);
	visit[i] = 1;

	if (N->father != -1 && visit[N->father] == 0) {
		temp[t] = N->father;
		++t;
		if (N->to_father_edge->h != realh[i][1]) {
			++(error[0]);
			printf ("error h[%d][f] is %d.\n", i, N->to_father_edge->h);
		}
	}
	if (N->mother != -1 && visit[N->mother] == 0) {
		temp[t] = N->mother;
		++t;
		if (N->to_mother_edge->h != realh[i][0]) {
			++(error[0]);
			printf ("error h[%d][m] is %d.\n", i, N->to_mother_edge->h);
		}
	}
	if (N->num_children != 0) for (c = 1; c <= N->num_children; ++c) {
		child = N->children[c];
		K = population[child];
		if (visit[child] == 0) {
			temp[t] = child;
			++t;
			if (N->sex == MALE && K->to_father_edge->h != realh[child][1]) {
				++(error[0]);
				printf ("error h[%d][f] is %d.\n", child, K->to_father_edge->h);
			}
			else if (N->sex == FEMALE && K->to_mother_edge->h != realh[child][0]) {
				++(error[0]);
				printf ("error h[%d][m] is %d.\n", child, K->to_mother_edge->h);
			}
		}
	}

	for (n = 1; n < t; ++n) DFSH (population, temp[n], realh, visit, error);
	free(temp);
}

void DFSP (PERSON **population, int i, int j, BIT_VALUE **p, BIT_VALUE **realp, int *visit, int *error) 
{
	int c, child, t = 1, n;
	int *temp;
	PERSON *N;
	N = population[i];
	temp = (int *) memset1D (2+N->num_children, sizeof(int), 0);

	if (N->father != -1 && visit[N->father] == 0) {
		visit[N->father] = 1;
		temp[t] = N->father;
		++t;
		if (p[N->father][j] != realp[N->father][j]) ++(error[1]);
	}
	if (N->mother != -1 && visit[N->mother] == 0) {
		visit[N->mother] = 1;
		temp[t] = N->mother;
		++t;
		if (p[N->mother][j] != realp[N->mother][j]) ++(error[1]);
	}
	for (c = 1; c <= N->num_children; ++c) {
		child = N->children[c];
		if (visit[child] == 0) {
			visit[child] = 1;
			temp[t] = child;
			++t;
			if (p[child][j] != realp[child][j]) ++(error[1]);
		}
	}

	for (n = 1; n < t; ++n) DFSP (population, temp[n], j, p, realp, visit, error);
	free (temp);
}

