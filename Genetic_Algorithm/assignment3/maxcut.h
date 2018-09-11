#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

//structures 
typedef struct edge {
    int v1_num;
    int v2_num;
    int value;
}Edges;
Edges *in_edge;

typedef struct neighbor {
    int num;
    int weight;
}Neighbor;

typedef struct vertex{
    int index;
    int num;
    Neighbor *n_vertex;
}Vertex;
Vertex *ver_n;

typedef struct sol {
    int *genetic;
    int score;
}genetic_solution;

genetic_solution *gene_sol;
genetic_solution selection_1;
genetic_solution selection_2;
genetic_solution optimal_sol;
//parameter
int population_num = 1000;
int select_mode = 2; //0,1,2
int cross_mode = 1;  //1,2
int cut_num;
float uniform_threshold;
float mutation_threshold = 0.03;
int exec_time = 495;
int replacement_mode = 1;
char *input_file_name = "maxcut.in";
char *output_file_name = "maxcut.out";
// global variables
int k_loop = 3;
int k_index;
int score_1;
int score_2;
int *genetic_1;
int *genetic_2;
int vertex_num;
int edge_num;
int loop_num = 0;
int start_time = 0;
int sol1=-1,sol2=-1;
int *score_list;
float *fit;
int *index_sort;
int *cut_points;
int min_index,max_index;
int min2_index,max2_index;
// functions 
void file_read(char *);
void initialize_gene(int);
void set_vertex_neighbor();
int calculate_maxcut_score(genetic_solution *);
int calculate_maxcut_score_2(int *);
void array_index_sort();
void array_sort();
int randomize(int);
void maxcut_solver();
void free_all(char *,char *);
void selection(int);
void cross_over(int);
int check_cutpoint(int);
void mutation();
void replacement(int);
int max_score_maxcut();
void file_write(char *);
void print_gene(genetic_solution *);
void print_gene_2(int *);
void local_optimization();
void local_optimization_with_2flip();
