#ifndef	_ALGO_H_
#define	_ALGO_H_

#include "kchi.h"

void nuclearFamily		(int *first, PERSON **population, int num_people);
void initialP			(GENOTYPE **g, BIT_VALUE **p, PERSON **population, int num_people, int num_locus);
void initialD			(GENOTYPE **g, BIT_VALUE ***d, PERSON **population, int num_people, int num_locus);
void spanningTree		(PERSON **population, int **global, int num_people);
void span 				(int i, int *visit, PERSON **population, int **global, int num_people);
void DFS 				(int i, int j, int root, GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p,
						PERSON **population, int *global, int *sumD, int *visit, char mode);
void DFSW				(int r, int i, PERSON **population, int *visit);
void assignFounder		(PERSON **population, int *first, int num_people);
void treeEQ				(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
						int *first, int **global, int num_people, int num_locus);
void localEQ 			(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
						int *first, int **local, int num_people, int num_locus);
void globalCEQ			(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population,  
						int **global, int num_people, int num_locus); 
void globalPEQ			(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population,  
						int *first, int **global, int **odd, int num_people, int num_locus, char mode); 
int  treeConstraint		(int i, int r, int j, int t, int root, GENOTYPE **g, BIT_VALUE **p, 
						PERSON **population, int *sumD, int *temp, BIT_VALUE dd);
void cycleConstraint	(int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
						int *sumD, int *temp, BIT_VALUE dd, BIT_VALUE ***d);
void pathConstraint		(int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
						int *sumD, int *temp, BIT_VALUE dd);
void pathTransform		(int i, int r, int j, int t, int *global, BIT_VALUE **p, PERSON **population, 
						int *sumD, int *temp, BIT_VALUE dd, BIT_VALUE ***d);
void oddpathstocycle 	(PERSON **population, int **oddpc, int **odd, int num_people); 
int  DFSO 				(PERSON **population, int R1, int RP, BIT_VALUE o, BIT_VALUE b, int **oddpc, int *visit);
void linkTree			(GENOTYPE **g, BIT_VALUE **p, PERSON **population, int *first, int **local, 
						int **global, int **odd, int num_people, int num_locus);
int  isLink				(PERSON **population, int **local, int **global, int **odd, int n, int m);
void calculateP 		(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
						int num_people, int num_locus);
void solveP 			(GENOTYPE **g, BIT_VALUE ***d, BIT_VALUE **p, PERSON **population, 
						int i, int j, int *visit);						
void statisitics 		(GENOTYPE **g, BIT_VALUE **p, PERSON **population, BIT_VALUE **realh, BIT_VALUE **realp, 
						int *first, int num_missing, int num_people, int num_locus);
void DFSH 				(PERSON **population, int i, BIT_VALUE **realh, int *visit, int *error);
void DFSP 				(PERSON **population, int i, int j, BIT_VALUE **p, BIT_VALUE **realp, int *visit, int *error);
#endif
