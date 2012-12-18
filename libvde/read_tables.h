/*
 * read_tables.h : routines to read and store data 
 * 		   from a user-defined table
 * */

int get_lines(FILE*);

void invertA(void);
void invertH(void);

void printVectors(double,double,int);

void read_all_interpolation_tables(void);
void read_interp_table(char*);
