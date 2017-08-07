#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "func.h"

void *memset1D (int l, int size, int a) {
	int i;
	int *x = malloc (l * size);
	for (i = 0; i < l; ++i) x[i] = a;
	return x;
}

void *memset2D (int l, int w, int size, int a) {
	int i, j;
	int **x = malloc (l * sizeof(void *));
	for (i = 0; i < l; ++i) {
		x[i] = malloc (w * size);
		for (j = 0; j < w; ++j) x[i][j] = a;
	}
	return x;
}

void *memset3D (int l, int w, int h, int size, int a) {
	int i, j, k;
	int ***x = malloc (l * sizeof(void **));
	for (i = 0; i < l; ++i)
		x[i] = malloc (w * sizeof(void *));
	for (i = 0; i < l; ++i) for (j = 0; j < w; ++j) {
		x[i][j] = malloc (h * size);
		for (k = 0; k < h; ++k) x[i][j][k] = a;
	}
	return x;
}

void free2D(void *x, int m) {
	int i, **y = (int **)x;
	for (i = 0; i < m; ++i)
		free(y[i]);
	free(x);
}

void free3D(void *x, int m, int n) {
	int i, j, ***y = (int ***)x;
	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			free(y[i][j]);
	for (i = 0; i < m; ++i)
		free(y[i]);
	free(x);
}

EDGE *new_edge () {
	EDGE *edge = (EDGE *) malloc (sizeof(EDGE));
	edge->type = GLOBAL_CROSS;
	edge->h = UNKNOWN;
	edge->d = UNKNOWN;
	return edge;
}

PERSON *new_person () {
	PERSON *person = (PERSON *) malloc (sizeof(PERSON));
	person->sex = BABY;
	person->W = UNKNOWN;
	person->relationship = (int *) memset1D (INIT_RELATIONSHIP_NUM, sizeof(int), 0);
	person->num_relationship = 0;
	person->cap_relationship = INIT_RELATIONSHIP_NUM;
	person->children = (int *) memset1D (INIT_CHILDREN_NUM, sizeof(int), 0);
	person->num_children = 0;
	person->cap_children = INIT_CHILDREN_NUM;
	person->constraint = (int **) memset2D (INIT_CONSTRAINT_NUM, 2, sizeof(int), -1); 
	//constraint[c][0] = another node, constraint[c][1] = b
	person->num_constraint = 0;
	person->cap_constraint = INIT_CONSTRAINT_NUM;
	person->father = -1;
	person->mother = -1;
	person->to_father_edge = new_edge();
	person->to_mother_edge = new_edge();
	person->root = -1;
	return person;
}

void readPedigree (FILE *ped, PERSON **population, int num_people) {
	int i, flag = 0;
	char line[MAX_LINE], *str;

	for (i = 1; i <= num_people; ++i) {
		fgets(line, MAX_LINE, ped);
		str = strtok(line, "\t");   // the first one is index, skip it
		
		/* Second column: gender */
		str = strtok(NULL, "\t");
		if (str[0] == 'M' || str[0] == 'm')
			population[i]->sex = MALE;
		else if (str[0] == 'F' || str[0] == 'f')
			population[i]->sex = FEMALE;
		else {
			++flag;
			fprintf (stderr, "The gender of individual %d is wrong\n", i);
		}

		/* Third and fourth column: father and mother */
		str = strtok(NULL, "\t");
		if (str[0] != 'N' && str[0] != 'n')
			population[i]->father = atoi(str);

		str = strtok(NULL, "\t");
		if (str[0] != 'N' && str[0] != 'n')
			population[i]->mother = atoi(str);

	}
}

int readGenotype (FILE *gen, GENOTYPE **g, int num_people, int num_locus) {
	int i, j, tmp, flag = 0;

	for (i = 1; i <= num_people; ++i) {
		for (j = 1; j <= num_locus; ++j) {
			fscanf(gen, "%d", &tmp);
			if (tmp == 0) g[i][j] = G00;
			else if (tmp == 1) g[i][j] = G11;
			else if (tmp == 2) g[i][j] = G01;
			else {
				++flag;
				g[i][j] = Gxx;
			}
		}
	}
 	return flag;
}

void readSolution
	(FILE *solp, FILE *solh, BIT_VALUE **realh, BIT_VALUE **realp, int num_people, int num_locus) {
	int i, j, tmp;

	for (i = 1; i <= num_people; ++i) {
		for (j = 1; j <= num_locus; ++j) {
			fscanf(solp, "%d", &tmp);
			if (tmp == 0) realp[i][j] = B0;
			else if (tmp == 1) realp[i][j] = B1;
			else realp[i][j] = FREE;
		}
	}
	for (i = 1; i <= num_people; ++i) {
		fscanf(solh, "%d", &tmp);
		if (tmp == 0) realh[i][1] = B0;
		else if (tmp == 1) realh[i][1] = B1;
		else realh[i][1] = FREE;	
		fscanf(solh, "%d", &tmp);
		if (tmp == 0) realh[i][0] = B0;
		else if (tmp == 1) realh[i][0] = B1;
		else realh[i][0] = FREE;
	}
}

void preRecover
	(GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C, PERSON **population, GENOTYPE **g, 
	 int num_locus, int num_people) {
	int i, j, c, child, stable, mod;
	GENOTYPE temp;
	PERSON *N;

	//missing alleles assignment from parents and children
	for (j = 1; j <= num_locus; ++j) {
		stable = num_people;
		mod = 0;
		while (stable != mod) {
			stable = mod;
			for (i = 1; i <= num_people; ++i) if (g[i][j] > 4 || g[i][j] == 2) {
				N = population[i];
				for (c = 1; c <= N->num_children; ++c) {
					child = N->children[c];
					temp = g[i][j];
					if (N->sex == MALE) g[i][j] = C2F[g[child][j]][2][g[i][j]];
					else g[i][j] = C2M[g[child][j]][2][g[i][j]];
					if (temp != g[i][j]) ++mod;
					if (g[i][j] < 5 && g[i][j] != 2) break;	
				} 
				if (N->father == -1 || N->mother == -1) continue;
				temp = g[i][j];
				if (g[i][j] == 2 || g[i][j] > 4)
					g[i][j] = F2C[g[N->father][j]][2][g[i][j]];
				if (g[i][j] == 2 || g[i][j] > 4)
					g[i][j] = M2C[g[N->mother][j]][2][g[i][j]];
				if (temp != g[i][j]) ++mod;
			}
		}
	}
}

void recover (GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C, 
	PERSON **population, GENOTYPE **g, BIT_VALUE **p, int num_people, int num_locus) 
{
	int i, j, c, child, stable, fix;
	BIT_VALUE h;
	GENOTYPE temp;
	PERSON *N, *K;

	for (j = 1; j <= num_locus; ++j) {
		stable = 1;
		fix = 0;
		while (stable != fix) {
			stable = fix;
			for (i = 1; i <= num_people; ++i) if (g[i][j] > 4 || g[i][j] == 2) {
				N = population[i];
				if (N->father != -1) {
					temp = g[i][j];
					h = N->to_father_edge->h;
					if (h == FREE) h = UNKNOWN;
					g[i][j] = F2C[g[N->father][j]][h][g[i][j]];
					fix = errorfix (g[i], temp, fix, j);	
					if (g[i][j] == G00 || g[i][j] == H01) p[i][j] = 0;
					else if (g[i][j] == G11 || g[i][j] == H10) p[i][j] = 1;
					if (g[i][j] < 5 && g[i][j] != 2) continue;
				}
				if (N->mother != -1) {
					temp = g[i][j];
					h = N->to_mother_edge->h;
					if (h == FREE) h = UNKNOWN;
					g[i][j] = M2C[g[N->mother][j]][h][g[i][j]];
					fix = errorfix (g[i], temp, fix, j);	
					if (g[i][j] == G00 || g[i][j] == H01) p[i][j] = 0;
					else if (g[i][j] == G11 || g[i][j] == H10) p[i][j] = 1;
					if (g[i][j] < 5 && g[i][j] != 2) continue;
				}
				for (c = 1; c <= N->num_children; ++c) {
					temp = g[i][j];
					child = N->children[c];
					K = population[child];
					if (N->sex == MALE) {
						h = K->to_father_edge->h;
						if (h == FREE) h = UNKNOWN;
						g[i][j] = C2F[g[child][j]][h][g[i][j]];
					} else {
						h = K->to_mother_edge->h;
						if (h == FREE) h = UNKNOWN;
						g[i][j] = C2M[g[child][j]][h][g[i][j]];
					}
					fix = errorfix (g[i], temp, fix, j);	
					if (g[i][j] == G00 || g[i][j] == H01) p[i][j] = 0;
					else if (g[i][j] == G11 || g[i][j] == H10) p[i][j] = 1;
					if (g[i][j] < 5 && g[i][j] != 2) break;
				}
			}
		}
	}
}

int errorfix (GENOTYPE *g, GENOTYPE temp, int fix, int j) {

	if (g[j] == ERR) {
		printf ("propogation error?\n");
		g[j] = temp;
	} else if (temp != g[j]) ++fix;

	return fix;
}

void haploReference(GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C) {
	//table[G = G00,G11,H01,H10,H1x,Hx1,H0x,Hx0,G01,Gxx,G1x,G0x][H = 0,1,2(unknown)]
	//	   [G = G01,Gxx,G1x,G0x,H1x,Hx1,H0x,Hx0]
	int g;
	//reference child to father
	for (g = 0; g < 12; ++g) {
		if (g == G00 || g == H01 || g == H0x) {
			//child get 0 from dad
			C2F[g][0][G01] = H01;		C2F[g][1][G01] = H10;		C2F[g][2][G01] = G01;
			C2F[g][0][G0x] = H0x;		C2F[g][1][G0x] = Hx0;		C2F[g][2][G0x] = G0x;
			C2F[g][0][G1x] = H01;		C2F[g][1][G1x] = H10;		C2F[g][2][G1x] = G01;
			C2F[g][0][Gxx] = H0x;		C2F[g][1][Gxx] = Hx0;		C2F[g][2][Gxx] = G0x;
			C2F[g][0][H1x] = ERR;		C2F[g][1][H1x] = H10;		C2F[g][2][H1x] = G01;
			C2F[g][0][Hx1] = H01;		C2F[g][1][Hx1] = ERR;		C2F[g][2][Hx1] = G01;
			C2F[g][0][H0x] = H0x;		C2F[g][1][H0x] = G00;		C2F[g][2][H0x] = H0x;	
			C2F[g][0][Hx0] = G00;		C2F[g][1][Hx0] = Hx0;		C2F[g][2][Hx0] = Hx0;	
		} else if (g == G11 || g == H10 || g == H1x) {
			//child get 1 from dad
			C2F[g][0][G01] = H10;		C2F[g][1][G01] = H01;		C2F[g][2][G01] = G01;
			C2F[g][0][G0x] = H10;		C2F[g][1][G0x] = H01;		C2F[g][2][G0x] = G01;
			C2F[g][0][G1x] = H1x;		C2F[g][1][G1x] = Hx1;		C2F[g][2][G1x] = G1x;
			C2F[g][0][Gxx] = H1x;		C2F[g][1][Gxx] = Hx1;		C2F[g][2][Gxx] = G1x;
			C2F[g][0][H1x] = H1x;		C2F[g][1][H1x] = G11;		C2F[g][2][H1x] = H1x;
			C2F[g][0][Hx1] = G11;		C2F[g][1][Hx1] = Hx1;		C2F[g][2][Hx1] = Hx1;
			C2F[g][0][H0x] = ERR;		C2F[g][1][H0x] = H01;		C2F[g][2][H0x] = G01;	
			C2F[g][0][Hx0] = H10;		C2F[g][1][Hx0] = ERR;		C2F[g][2][Hx0] = G01;	
		} else {
			//no information
			C2F[g][0][G01] = G01;		C2F[g][1][G01] = G01;		C2F[g][2][G01] = G01;
			C2F[g][0][G0x] = G0x;		C2F[g][1][G0x] = G0x;		C2F[g][2][G0x] = G0x;
			C2F[g][0][G1x] = G1x;		C2F[g][1][G1x] = G1x;		C2F[g][2][G1x] = G1x;
			C2F[g][0][Gxx] = Gxx;		C2F[g][1][Gxx] = Gxx;		C2F[g][2][Gxx] = Gxx;
			C2F[g][0][H1x] = H1x;		C2F[g][1][H1x] = H1x;		C2F[g][2][H1x] = H1x;
			C2F[g][0][Hx1] = Hx1;		C2F[g][1][Hx1] = Hx1;		C2F[g][2][Hx1] = Hx1;
			C2F[g][0][H0x] = H0x;		C2F[g][1][H0x] = H0x;		C2F[g][2][H0x] = H0x;	
			C2F[g][0][Hx0] = Hx0;		C2F[g][1][Hx0] = Hx0;		C2F[g][2][Hx0] = Hx0;	
		}
	}
	//reference child to mother
	for (g = 0; g < 12; ++g) {
		if (g == G00 || g == H10 || g == Hx0) {
			//child get 0 from mom
			C2M[g][0][G01] = H01;		C2M[g][1][G01] = H10;		C2M[g][2][G01] = G01;
			C2M[g][0][G0x] = H0x;		C2M[g][1][G0x] = Hx0;		C2M[g][2][G0x] = G0x;
			C2M[g][0][G1x] = H01;		C2M[g][1][G1x] = H10;		C2M[g][2][G1x] = G01;
			C2M[g][0][Gxx] = H0x;		C2M[g][1][Gxx] = Hx0;		C2M[g][2][Gxx] = G0x;
			C2M[g][0][H1x] = ERR;		C2M[g][1][H1x] = H10;		C2M[g][2][H1x] = G01;
			C2M[g][0][Hx1] = H01;		C2M[g][1][Hx1] = ERR;		C2M[g][2][Hx1] = G01;
			C2M[g][0][H0x] = H0x;		C2M[g][1][H0x] = G00;		C2M[g][2][H0x] = H0x;	
			C2M[g][0][Hx0] = G00;		C2M[g][1][Hx0] = Hx0;		C2M[g][2][Hx0] = Hx0;	
		} else if (g == G11 || g == H01 || g == Hx1) {
			//child get 1 from mom
			C2M[g][0][G01] = H10;		C2M[g][1][G01] = H01;		C2M[g][2][G01] = G01;
			C2M[g][0][G0x] = H10;		C2M[g][1][G0x] = H01;		C2M[g][2][G0x] = G01;
			C2M[g][0][G1x] = H1x;		C2M[g][1][G1x] = Hx1;		C2M[g][2][G1x] = G1x;
			C2M[g][0][Gxx] = H1x;		C2M[g][1][Gxx] = Hx1;		C2M[g][2][Gxx] = G1x;
			C2M[g][0][H1x] = H1x;		C2M[g][1][H1x] = G11;		C2M[g][2][H1x] = H1x;
			C2M[g][0][Hx1] = G11;		C2M[g][1][Hx1] = Hx1;		C2M[g][2][Hx1] = Hx1;
			C2M[g][0][H0x] = ERR;		C2M[g][1][H0x] = H01;		C2M[g][2][H0x] = G01;	
			C2M[g][0][Hx0] = H10;		C2M[g][1][Hx0] = ERR;		C2M[g][2][Hx0] = G01;	
		} else {
			//no information
			C2M[g][0][G01] = G01;		C2M[g][1][G01] = G01;		C2M[g][2][G01] = G01;
			C2M[g][0][G0x] = G0x;		C2M[g][1][G0x] = G0x;		C2M[g][2][G0x] = G0x;
			C2M[g][0][G1x] = G1x;		C2M[g][1][G1x] = G1x;		C2M[g][2][G1x] = G1x;
			C2M[g][0][Gxx] = Gxx;		C2M[g][1][Gxx] = Gxx;		C2M[g][2][Gxx] = Gxx;
			C2M[g][0][H1x] = H1x;		C2M[g][1][H1x] = H1x;		C2M[g][2][H1x] = H1x;
			C2M[g][0][Hx1] = Hx1;		C2M[g][1][Hx1] = Hx1;		C2M[g][2][Hx1] = Hx1;
			C2M[g][0][H0x] = H0x;		C2M[g][1][H0x] = H0x;		C2M[g][2][H0x] = H0x;	
			C2M[g][0][Hx0] = Hx0;		C2M[g][1][Hx0] = Hx0;		C2M[g][2][Hx0] = Hx0;	
		}
	}
	//reference father/mother to child
	for (g = 0; g < 12; ++g) {
		if (g == G00 || g == G11 || g == H01 || g == H10 || 
			g == H1x || g == Hx1 || g == H0x || g == Hx0) {
			if (g == G00 || g == H01 || g == H0x) {
				//trans 0 from left
				F2C[g][0][G01] = H01;		F2C[g][0][G0x] = H0x;		
				F2C[g][0][G1x] = H01;		F2C[g][0][Gxx] = H0x;		
				F2C[g][0][H1x] = ERR;		F2C[g][0][Hx1] = H01;		
				F2C[g][0][H0x] = H0x;		F2C[g][0][Hx0] = G00;			

				M2C[g][0][G01] = H10;		M2C[g][0][G0x] = Hx0;		
				M2C[g][0][G1x] = H10;		M2C[g][0][Gxx] = Hx0;		
				M2C[g][0][H1x] = H10;		M2C[g][0][Hx1] = ERR;		
				M2C[g][0][H0x] = G00;		M2C[g][0][Hx0] = Hx0;
			} 
			if (g == G11 || g == H10 || g == H1x) {
				//trans 1 from left
				F2C[g][0][G01] = H10;		F2C[g][0][G0x] = H10;		
				F2C[g][0][G1x] = H1x;		F2C[g][0][Gxx] = H1x;		
				F2C[g][0][H1x] = H1x;		F2C[g][0][Hx1] = G11;		
				F2C[g][0][H0x] = ERR;		F2C[g][0][Hx0] = H10;	

				M2C[g][0][G01] = H01;		M2C[g][0][G0x] = H01;		
				M2C[g][0][G1x] = Hx1;		M2C[g][0][Gxx] = Hx1;		
				M2C[g][0][H1x] = G11;		M2C[g][0][Hx1] = Hx1;		
				M2C[g][0][H0x] = H01;		M2C[g][0][Hx0] = ERR;
			} 
			if (g == G00 || g == H10 || g == Hx0) {
				//trans 0 from right
				F2C[g][1][G01] = H01;		F2C[g][1][G0x] = H0x;		
				F2C[g][1][G1x] = H01;		F2C[g][1][Gxx] = H0x;		
				F2C[g][1][H1x] = ERR;		F2C[g][1][Hx1] = H01;		
				F2C[g][1][H0x] = H0x;		F2C[g][1][Hx0] = G00;			

				M2C[g][1][G01] = H10;		M2C[g][1][G0x] = Hx0;		
				M2C[g][1][G1x] = H10;		M2C[g][1][Gxx] = Hx0;		
				M2C[g][1][H1x] = H10;		M2C[g][1][Hx1] = ERR;		
				M2C[g][1][H0x] = G00;		M2C[g][1][Hx0] = Hx0;
			}
			if (g == G11 || g == H01 || g == Hx1) {
				//trans 1 from right
				F2C[g][1][G01] = H10;		F2C[g][1][G0x] = H10;		
				F2C[g][1][G1x] = H1x;		F2C[g][1][Gxx] = H1x;		
				F2C[g][1][H1x] = H1x;		F2C[g][1][Hx1] = G11;		
				F2C[g][1][H0x] = ERR;		F2C[g][1][Hx0] = H10;	

				M2C[g][1][G01] = H01;		M2C[g][1][G0x] = H01;		
				M2C[g][1][G1x] = Hx1;		M2C[g][1][Gxx] = Hx1;		
				M2C[g][1][H1x] = G11;		M2C[g][1][Hx1] = Hx1;		
				M2C[g][1][H0x] = H01;		M2C[g][1][Hx0] = ERR;
			} 
			if (g == G00) {
				//trans 0 
				F2C[g][2][G01] = H01;		F2C[g][2][G0x] = H0x;		
				F2C[g][2][G1x] = H01;		F2C[g][2][Gxx] = H0x;		
				F2C[g][2][H1x] = ERR;		F2C[g][2][Hx1] = H01;		
				F2C[g][2][H0x] = H0x;		F2C[g][2][Hx0] = G00;			

				M2C[g][2][G01] = H10;		M2C[g][2][G0x] = Hx0;		
				M2C[g][2][G1x] = H10;		M2C[g][2][Gxx] = Hx0;		
				M2C[g][2][H1x] = H10;		M2C[g][2][Hx1] = ERR;		
				M2C[g][2][H0x] = G00;		M2C[g][2][Hx0] = Hx0;
			}
			if (g == G11) {
				//trans 1
				F2C[g][2][G01] = H10;		F2C[g][2][G0x] = H10;		
				F2C[g][2][G1x] = H1x;		F2C[g][2][Gxx] = H1x;		
				F2C[g][2][H1x] = H1x;		F2C[g][2][Hx1] = G11;		
				F2C[g][2][H0x] = ERR;		F2C[g][2][Hx0] = H10;	

				M2C[g][2][G01] = H01;		M2C[g][2][G0x] = H01;		
				M2C[g][2][G1x] = Hx1;		M2C[g][2][Gxx] = Hx1;		
				M2C[g][2][H1x] = G11;		M2C[g][2][Hx1] = Hx1;		
				M2C[g][2][H0x] = H01;		M2C[g][2][Hx0] = ERR;
			} 
			if (g == H0x || g == H1x) {
				//no information from right
				F2C[g][1][G01] = G01;		F2C[g][1][G0x] = G0x;		
				F2C[g][1][G1x] = G1x;		F2C[g][1][Gxx] = Gxx;		
				F2C[g][1][H1x] = H1x;		F2C[g][1][Hx1] = Hx1;		
				F2C[g][1][H0x] = H0x;		F2C[g][1][Hx0] = Hx0;	

				M2C[g][1][G01] = G01;		M2C[g][1][G0x] = G0x;		
				M2C[g][1][G1x] = G1x;		M2C[g][1][Gxx] = Gxx;		
				M2C[g][1][H1x] = H1x;		M2C[g][1][Hx1] = Hx1;		
				M2C[g][1][H0x] = H0x;		M2C[g][1][Hx0] = Hx0;
			}
			if (g == Hx0 || g == Hx1) {
				//no information from left
				F2C[g][0][G01] = G01;		F2C[g][0][G0x] = G0x;		
				F2C[g][0][G1x] = G1x;		F2C[g][0][Gxx] = Gxx;		
				F2C[g][0][H1x] = H1x;		F2C[g][0][Hx1] = Hx1;		
				F2C[g][0][H0x] = H0x;		F2C[g][0][Hx0] = Hx0;	

				M2C[g][0][G01] = G01;		M2C[g][0][G0x] = G0x;		
				M2C[g][0][G1x] = G1x;		M2C[g][0][Gxx] = Gxx;		
				M2C[g][0][H1x] = H1x;		M2C[g][0][Hx1] = Hx1;		
				M2C[g][0][H0x] = H0x;		M2C[g][0][Hx0] = Hx0;			
			}
			if (g == H10 || g == H01 || g == H1x || g == Hx1 || g == H0x || g == Hx0) {
				//no information
				F2C[g][2][G01] = G01;		F2C[g][2][G0x] = G0x;		
				F2C[g][2][G1x] = G1x;		F2C[g][2][Gxx] = Gxx;		
				F2C[g][2][H1x] = H1x;		F2C[g][2][Hx1] = Hx1;		
				F2C[g][2][H0x] = H0x;		F2C[g][2][Hx0] = Hx0;					

				M2C[g][2][G01] = G01;		M2C[g][2][G0x] = G0x;		
				M2C[g][2][G1x] = G1x;		M2C[g][2][Gxx] = Gxx;		
				M2C[g][2][H1x] = H1x;		M2C[g][2][Hx1] = Hx1;		
				M2C[g][2][H0x] = H0x;		M2C[g][2][Hx0] = Hx0;
			}
		} else {
			//no information
			F2C[g][0][G01] = G01;		F2C[g][1][G01] = G01;		F2C[g][2][G01] = G01;
			F2C[g][0][G0x] = G0x;		F2C[g][1][G0x] = G0x;		F2C[g][2][G0x] = G0x;
			F2C[g][0][G1x] = G1x;		F2C[g][1][G1x] = G1x;		F2C[g][2][G1x] = G1x;
			F2C[g][0][Gxx] = Gxx;		F2C[g][1][Gxx] = Gxx;		F2C[g][2][Gxx] = Gxx;
			F2C[g][0][H1x] = H1x;		F2C[g][1][H1x] = H1x;		F2C[g][2][H1x] = H1x;
			F2C[g][0][Hx1] = Hx1;		F2C[g][1][Hx1] = Hx1;		F2C[g][2][Hx1] = Hx1;
			F2C[g][0][H0x] = H0x;		F2C[g][1][H0x] = H0x;		F2C[g][2][H0x] = H0x;	
			F2C[g][0][Hx0] = Hx0;		F2C[g][1][Hx0] = Hx0;		F2C[g][2][Hx0] = Hx0;	

			M2C[g][0][G01] = G01;		M2C[g][1][G01] = G01;		M2C[g][2][G01] = G01;
			M2C[g][0][G0x] = G0x;		M2C[g][1][G0x] = G0x;		M2C[g][2][G0x] = G0x;
			M2C[g][0][G1x] = G1x;		M2C[g][1][G1x] = G1x;		M2C[g][2][G1x] = G1x;
			M2C[g][0][Gxx] = Gxx;		M2C[g][1][Gxx] = Gxx;		M2C[g][2][Gxx] = Gxx;
			M2C[g][0][H1x] = H1x;		M2C[g][1][H1x] = H1x;		M2C[g][2][H1x] = H1x;
			M2C[g][0][Hx1] = Hx1;		M2C[g][1][Hx1] = Hx1;		M2C[g][2][Hx1] = Hx1;
			M2C[g][0][H0x] = H0x;		M2C[g][1][H0x] = H0x;		M2C[g][2][H0x] = H0x;	
			M2C[g][0][Hx0] = Hx0;		M2C[g][1][Hx0] = Hx0;		M2C[g][2][Hx0] = Hx0;			
		}
	}
}

