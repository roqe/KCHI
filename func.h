#ifndef	_FUNC_H_
#define	_FUNC_H_

#define INIT_BUCKET_SIZE		50
#define INIT_CHILDREN_NUM		1
#define INIT_RELATIONSHIP_NUM	5
#define INIT_CONSTRAINT_NUM		1

#include "kchi.h"

void *memset1D		(int l, int size, int a);
void *memset2D		(int l, int w, int size, int a);
void *memset3D		(int l, int w, int h, int size, int a);

void free2D			(void *x, int m);
void free3D			(void *x, int m, int n);

void readPedigree	(FILE *ped, PERSON **population, int num_people);
int  readGenotype	(FILE *gen, GENOTYPE **g, int num_people, int num_locus);
void readSolution	(FILE *solp, FILE *solh, BIT_VALUE **realh, BIT_VALUE **realp, int num_people, int num_locus);

void preRecover		(GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C, PERSON **population, 
					GENOTYPE **g, int num_locus, int num_people);
void recover 		(GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C, 
					PERSON **population, GENOTYPE **g, BIT_VALUE **p, int num_people, int num_locus);
int  errorfix 		(GENOTYPE *ansg, GENOTYPE temp, int fix, int j);
void haploReference	(GENOTYPE ***C2F, GENOTYPE ***C2M, GENOTYPE ***F2C, GENOTYPE ***M2C);

EDGE		*new_edge();
PERSON		*new_person();

#endif
