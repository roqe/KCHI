#ifndef	_KCHI_H_
#define	_KCHI_H_

#define	MAX_LINE		1000
#define INIT_CYLCE_NUM	20

typedef enum {
	G00 = 0, G11 = 1, G01 = 2, H01 = 3, 
	H10 = 4, G1x = 5, G0x = 6, Gxx = 7,
	H0x = 8, H1x = 9, Hx0 =10, Hx1 =11,
	ERR =12
} GENOTYPE;

typedef enum {
	MALE = 0, FEMALE = 1, BABY = 2
} GENDER;

typedef enum {
	B0 = 0, B1 = 1, UNKNOWN = 2, FREE = 3
} BIT_VALUE;

typedef enum {
	TREE = 0, LOCAL_CROSS = 1, GLOBAL_CROSS = 2
} EDGE_TYPE;

typedef struct {
	EDGE_TYPE type;
	BIT_VALUE h;
	BIT_VALUE d; //if both ends are in the G*, then the h will not be used.
} EDGE;

typedef struct {
	GENDER sex;
	BIT_VALUE W;
	int *relationship, num_relationship, cap_relationship;
	int *children, num_children, cap_children;
	int **constraint, num_constraint, cap_constraint;
	int father, mother;
	EDGE *to_father_edge, *to_mother_edge; 	
	int root;
} PERSON;

#endif
