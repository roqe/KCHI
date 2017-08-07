//June 100803 version 1.0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "kchi.h"
#include "algo.h"
#include "func.h"

int main (int argc, char **argv) {

	int  i; 
	int  num_people, num_locus, num_missing, num_free = 0;
	int  *first, **local, **global, **odd; 
	FILE *ped, *gen, *solp, *solh;
	char line[MAX_LINE];

	GENOTYPE **g;
	GENOTYPE ***C2F, ***C2M, ***F2C, ***M2C;		
	BIT_VALUE **p, ***d, **realh, **realp;
	PERSON **population, *G0, *G1;

	if (argc != 5) {
		fprintf (stderr, "Usage:\n");
		fprintf (stderr, "> kchi <pedigree file> <genotype file> <p solution file> <h solution file>\n");
		exit(EXIT_FAILURE);
	}
	if ((ped = fopen (argv[1], "r")) == NULL) {
		fprintf (stderr, "Bad pedigree file name: %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	if ((gen = fopen (argv[2], "r")) == NULL) {
		fprintf (stderr, "Bad genotype file name: %s\n", argv[2]);
		exit(EXIT_FAILURE);
	}
	if ((solp = fopen (argv[3], "r")) == NULL) {
		fprintf (stderr, "Bad solution file name: %s\n", argv[3]);
		exit(EXIT_FAILURE);
	}
	if ((solh = fopen (argv[4], "r")) == NULL) {
		fprintf (stderr, "Bad solution file name: %s\n", argv[4]);
		exit(EXIT_FAILURE);
	}

	fgets(line, MAX_LINE, ped);
	fscanf(gen, "%d %d", &num_people, &num_locus);

	//-----initialization-----//
	g = (GENOTYPE **) memset2D (num_people+1, num_locus+1, sizeof(GENOTYPE), (int)Gxx);
	p = (BIT_VALUE **) memset2D (num_people+1, num_locus+1, sizeof(BIT_VALUE), (int)UNKNOWN);
	d = (BIT_VALUE ***) memset3D (num_people+1, num_locus+1, 2, sizeof(BIT_VALUE), (int)UNKNOWN);
	realh = (BIT_VALUE **) memset2D (num_people+1, 2, sizeof(BIT_VALUE), (int)UNKNOWN);
	realp = (BIT_VALUE **) memset2D (num_people+1, num_locus+1, sizeof(BIT_VALUE), (int)UNKNOWN);
	//d[n][m][0:mother/1:father]
	first = (int *) memset1D (num_people+1, sizeof(int), 0);
	local = (int **) memset2D (num_people+1, 2, sizeof(int), 0);
	//local[c][0:parent/1:child]
	global = (int **) memset2D (INIT_CYLCE_NUM, 3, sizeof(int), 0);
	//global[c][0:parent/1:child/2:leaf]
	odd = (int **) memset2D (INIT_CYLCE_NUM, 3, sizeof(int), 0);
	//odd[c][0:ccp 1/1:ccp 2/2:tree link]
	population = (PERSON **) malloc ((num_people+1) * sizeof(PERSON *));
	for (i = 1; i <= num_people; ++i) population[i] = new_person();

	//-----input data-----//
	readPedigree (ped, population, num_people);
	fclose (ped);		//output: person-sex,father,mother
	nuclearFamily(first, population, num_people);
						//output: relationship, children
	num_missing = readGenotype (gen, g, num_people, num_locus);
	fclose (gen);		//output: return the number of missing genotype
	readSolution(solp, solh, realh, realp, num_people, num_locus);
	fclose (solp);
	fclose (solh);	

	//-----pre-recover missing genotype-----//
	if (num_missing != 0) {
		//table[G = G00,G11,G01,H01,H10,Gxx,G1x,G0x,H1x,Hx1,H0x,Hx0][H = 0,1,2(unknown)][G]
		C2F = (GENOTYPE ***) memset3D (12, 3, 12, sizeof(GENOTYPE), (int)UNKNOWN);
		C2M = (GENOTYPE ***) memset3D (12, 3, 12, sizeof(GENOTYPE), (int)UNKNOWN);
		F2C = (GENOTYPE ***) memset3D (12, 3, 12, sizeof(GENOTYPE), (int)UNKNOWN);
		M2C = (GENOTYPE ***) memset3D (12, 3, 12, sizeof(GENOTYPE), (int)UNKNOWN);
		haploReference 	(C2F, C2M, F2C, M2C);	
						//output: haplotype reference table
		preRecover 		(C2F, C2M, F2C, M2C, population, g, num_locus, num_people);
						//output: partial recovered genotype
	}

	//-----main algorithm-----//
	initialP 	(g, p, population, num_people, num_locus);
						//output: predetermined nodes and their p values
	initialD 	(g, d, population, num_people, num_locus);
						//output: d values for on edges
	spanningTree	(population, global, num_people);		
						//output: type of edges, TREE, LOCAL CROSS or GLOBAL CROSS	
	treeEQ 		(g, d, p, population, first, global, num_people, num_locus);
						//output: connected components and tree constraints
	localEQ		(g, d, p, population, first, local, num_people, num_locus);
						//output: local cross constraints (path, cycle)

	if (global[0][0] != 0) {	
		/* deal with global cycles */
		globalCEQ	(g, d, p, population, global, num_people, num_locus); 	
						//output: global cycle constraints
		globalPEQ	(g, d, p, population, first, global, odd, num_people, num_locus, 'p'); 	
						//output: global path constraints and intrinsic cycle update
		globalPEQ	(g, d, p, population, first, global, odd, num_people, num_locus, 'f'); 	
						//output: global path constraints transform into tree constraint
		for (i = 1; i <= global[0][0]; ++i) {
			G0 = population[global[i][0]];
			G1 = population[global[i][1]];
			if ((G0->sex == MALE && G1->to_father_edge->h > 1) || 
					(G0->sex == FEMALE && G1->to_mother_edge->h > 1)) {
				++num_free;
/*				printf ("test in kchi.c! G0=%d, G1=%d, G0->root=%d, G1->root=%d, hf=%d, hm=%d.\n", 
						global[i][0], global[i][1], G0->root, G1->root, G1->to_father_edge->h, G1->to_mother_edge->h);
*/			}	
		}
	}
	if (num_missing != 0) printf ("number of free global h (m): %d\n", num_free);
	else printf ("number of free global h: %d\n", num_free);
	
	linkTree	(g, p, population, first, local, global, odd, num_people, num_locus);
						//output: link tree components, assign h
	calculateP 	(g, d, p, population, num_people, num_locus);
						//output: p solution
	
	if (num_missing != 0) {		//propagation algorithm
		recover	(C2F, C2M, F2C, M2C, population, g, p, num_people, num_locus);
						//output: fixed genotype data
	}
	
	statisitics (g, p, population, realh, realp, first, num_missing, num_people, num_locus);
						//output: confirm h confirm p

	//-----free-----//
	if (num_missing != 0){
		free3D(C2F,12, 3);
		free3D(C2M,12, 3);
		free3D(F2C,12, 3);
		free3D(M2C,12, 3);
	}

	free3D(d, num_people+1, num_locus+1);

	free2D(g, num_people+1);
	free2D(p, num_people+1);
	free2D(realh, num_people+1);
	free2D(realp, num_people+1);
	free2D(local, num_people+1);
	free2D(global, INIT_CYLCE_NUM);
	free2D(odd, INIT_CYLCE_NUM);

	free(first);
	free(population);
	
	return 0;
}
