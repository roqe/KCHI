#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "func.h"

typedef enum {A00 = 0, A11 = 1, A10 = 2, A01 = 3, Axx = 4} GENOTYPE_ANS;

int main (int argc, char **argv) {

	if (argc != 3 && argc != 5) {
		printf ("Usage:\n");
		printf ("> randomGen <file: pedigree> <int: locus number> <m: missing, optional> <d: default missing rate, 10%% > or <int: missing rate, the denominator is 100>\n");
		exit(EXIT_FAILURE);
	}
			//random generate founder answer from 00 01 10 11
			//generate children, output pedigree structure
			//get the whole g value and p value
			//output g input file and p solution file
			//a function in kchi to check p answers, ignore free p

	int node_num;
	int fodr_num = 0;
	int flag = 0;
	char line[1000], *str;
	FILE *ped;
	if ((ped = fopen (argv[1], "r")) == NULL) {
		fprintf (stderr, "Bad pedigree file name: %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	fscanf(ped, "%d ", &node_num);
	
	int locu_num = atoi(argv[2]);
	int *father = (int *) memset1D (node_num+1,  sizeof(int), -1);
	int *mother = (int *) memset1D (node_num+1,  sizeof(int), -1);
	int **first = (int **) memset2D (node_num+1, node_num+1, sizeof(int), 0);
	GENDER *gender = (GENDER *) malloc ((node_num+1) * sizeof(GENDER));	
	BIT_VALUE **h = (BIT_VALUE **) memset2D (node_num+1, 2, sizeof(BIT_VALUE), (int)UNKNOWN);
	GENOTYPE_ANS **a = (GENOTYPE_ANS **) memset2D (node_num+1, locu_num+1, sizeof(GENOTYPE_ANS), (int)Axx);
	FILE *outg, *outp, *outh;
	int i, j, assign, g, p, fix = 0, stable = node_num;

	//read pedigree
	srand (time(0));

	for (i = 1; i <= node_num; ++i) {
		fgets(line, 1000, ped);
		str = strtok(line, "\t"); 
		str = strtok(NULL, "\t");
		if (str[0] == 'M' || str[0] == 'm')
			gender[i] = MALE;
		else if (str[0] == 'F' || str[0] == 'f')
			gender[i] = FEMALE;
		else {
			++flag;
			fprintf (stderr, "The gender of individual %d is wrong\n", i);
		}

		/* Third and fourth column: father and mother */
		str = strtok(NULL, "\t");
		if (str[0] == 'N' || str[0] == 'n'){
			father[i] = -1;
			++fodr_num;
		}
		else
			father[i] = atoi(str);

		str = strtok(NULL, "\t");
		if (str[0] == 'N' || str[0] == 'n')
			mother[i] = -1;
		else
			mother[i] = atoi(str);

		if (father[i] != -1 && first[father[i]][mother[i]] == 0) 
			first[father[i]][mother[i]] = i;
	}
	
	//work out haplotype
	for (j = 1; j <= locu_num; ++j) {
		for (i = 1; i <= node_num ; ++i) {
			if (father[i] == -1) {
				assign = rand () %4;
				if (assign == 0) a[i][j] = A00;
				else if (assign == 1) a[i][j] = A11;
				else if (assign == 2) a[i][j] = A10;
				else if (assign == 3) a[i][j] = A01;
				//random set founder
			}
		}
	}
			
			
	for (j = 1; j <= locu_num; ++j) {
		stable = node_num;
		fix = 0;
		while (stable != fix) {
			stable = fix;		
			for (i = 1; i <= node_num ; ++i) 
				if (a[i][j] == Axx && a[father[i]][j] != Axx && a[mother[i]][j] != Axx) {
				++fix;
				//first child of founder give h = 0 
				if (father[father[i]] == -1 && i == first[father[i]][mother[i]]) 
					h[i][1] = 0;
				if (father[mother[i]] == -1 && i == first[father[i]][mother[i]]) 
					h[i][0] = 0;
				//homolygous
	 			if (a[father[i]][j] == A00 && a[mother[i]][j] == A00) a[i][j] = A00;
				else if (a[father[i]][j] == A11 && a[mother[i]][j] == A11) a[i][j] = A11;
				else if (a[father[i]][j] == A00 && a[mother[i]][j] == A11) a[i][j] = A01;
				else if (a[father[i]][j] == A11 && a[mother[i]][j] == A00) a[i][j] = A10;
				//heterolygous
				else if (a[father[i]][j] == A00) {
					if (h[i][0] == UNKNOWN) h[i][0] = rand() %2;
					if ((a[mother[i]][j] == A01 && h[i][0] == 0)||
						(a[mother[i]][j] == A10 && h[i][0] == 1)) a[i][j] = A00;
					else a[i][j] = A01;
				} else if (a[father[i]][j] == A11) {
					if (h[i][0] == UNKNOWN) h[i][0] = rand() %2;
					if ((a[mother[i]][j] == A01 && h[i][0] == 0)||
						(a[mother[i]][j] == A10 && h[i][0] == 1)) a[i][j] = A10;
					else a[i][j] = A11;
				} else if (a[mother[i]][j] == A00) {
					if (h[i][1] == UNKNOWN) h[i][1] = rand() %2;
					if ((a[father[i]][j] == A01 && h[i][1] == 0)||
						(a[father[i]][j] == A10 && h[i][1] == 1)) a[i][j] = A00;
					else a[i][j] = A10;
				} else if (a[mother[i]][j] == A11) {
					if (h[i][1] == UNKNOWN) h[i][1] = rand() %2;
					if ((a[father[i]][j] == A01 && h[i][1] == 0)||
						(a[father[i]][j] == A10 && h[i][1] == 1)) a[i][j] = A01;
					else a[i][j] = A11;
				} else {
					if (h[i][0] == UNKNOWN) h[i][0] = rand() %2;
					if (h[i][1] == UNKNOWN) h[i][1] = rand() %2;
					if (h[i][0] == 0 && h[i][1] == 0) {
						if (a[father[i]][j] == A01 && a[mother[i]][j] == A10) 
							a[i][j] = A01;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A01) 
							a[i][j] = A10;
						else if (a[father[i]][j] == A01 && a[mother[i]][j] == A01)
							a[i][j] = A00;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A10) 
							a[i][j] = A11;
					} else if (h[i][0] == 1 && h[i][1] == 1) {
						if (a[father[i]][j] == A01 && a[mother[i]][j] == A10) 
							a[i][j] = A10;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A01) 
							a[i][j] = A01;
						else if (a[father[i]][j] == A01 && a[mother[i]][j] == A01)
							a[i][j] = A11;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A10) 
							a[i][j] = A00;
					} else if (h[i][0] == 0 && h[i][1] == 1) {
						if (a[father[i]][j] == A01 && a[mother[i]][j] == A10)
							a[i][j] = A11;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A01) 
							a[i][j] = A00;
						else if (a[father[i]][j] == A01 && a[mother[i]][j] == A01)
							a[i][j] = A10;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A10) 
							a[i][j] = A01;
					} else if (h[i][0] == 1 && h[i][1] == 0) {
						if (a[father[i]][j] == A01 && a[mother[i]][j] == A10) 
							a[i][j] = A00;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A01) 
							a[i][j] = A11;
						else if (a[father[i]][j] == A01 && a[mother[i]][j] == A01)
							a[i][j] = A01;
						else if (a[father[i]][j] == A10 && a[mother[i]][j] == A10) 
							a[i][j] = A10;
					}
				} 
			}
		}
	}

	if ((outg = fopen("gdata", "w")) == NULL) {
		fprintf (stderr, "Bad filename: %s", argv[1]);
		return 0;
	}
	fprintf (outg, "%d %d\n", node_num, locu_num);
	for (i = 1; i <= node_num; ++i) {
		for (j = 1; j <= locu_num; ++j) {
			if (a[i][j] == A00) g = 0;
			else if (a[i][j] == A11) g = 1;
			else g = 2;
			fprintf (outg, "%d ", g);
		}
		fprintf (outg, "\n");
	}
	fclose(outg);

	if ((outp = fopen("psdata", "w")) == NULL) {
		fprintf (stderr, "Bad filename: %s", argv[1]);
		return 0;
	}
	for (i = 1; i <= node_num; ++i) {
		for (j = 1; j <= locu_num; ++j) {
			if (a[i][j] == A00 || a[i][j] == A01) p = 0;
			else if (a[i][j] == A11 || a[i][j] == A10) p = 1;
			else p = 2;
			fprintf (outp, "%d ", p);
		}
		fprintf (outp, "\n");
	}
	fclose(outp);

	if ((outh = fopen("hsdata", "w")) == NULL) {
		fprintf (stderr, "Bad filename: %s", argv[1]);
		return 0;
	}
	for (i = 1; i <= node_num; ++i) {
        fprintf(outh, "%d %d",h[i][1], h[i][0]);
		fprintf(outh, "\n");
	}
	fclose(outh);

	if (argc == 5 && *argv[3] == 'm') {
		FILE *outm;
		int m, rate;
		if (*argv[4] == 'd') rate = 10;
		else rate = atoi(argv[4]);
		if ((outm = fopen("mgdata", "w")) == NULL) {
			fprintf (stderr, "Bad filename: %s", argv[1]);
			return 0;
		}
		fprintf (outm, "%d %d\n", node_num, locu_num);
		for (i = 1; i <= node_num; ++i) {
			for (j = 1; j <= locu_num; ++j) {
				if (a[i][j] == A00) m = 0;
				else if (a[i][j] == A11) m = 1;
				else m = 2;
				if (rand()%100 < rate) m = 3; //Gxx
				fprintf (outm, "%d ", m);
			}
			fprintf (outm, "\n");
		}
		fclose(outm);	
	}

	free(father);
	free(mother);
	free(gender);
	free2D(first, node_num+1);
	free2D(h, node_num+1);
	free2D(a, node_num+1);

	return 0;
}
