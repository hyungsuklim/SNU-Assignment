#include "maxcut.h"


int main(void) 
{
    int first;
    int last;
    first = clock();
    srand((unsigned)time(NULL));
    file_read(input_file_name);
    start_time = clock();
    initialize_gene(population_num);
    maxcut_solver();
    file_write(output_file_name);
    //printf("Write Complete\n");
    last = clock();
    //printf("Time is : %.0f\n",(float)(last-first)/CLOCKS_PER_SEC);
    return 0;
}

void file_read(char *file_name)
{
    int i;
    int V_num,E_num;
    int warn_1,warn_2;
    FILE *maxcut_in;
    maxcut_in = fopen(file_name,"r");
    if(maxcut_in)
    {
        warn_1 = fscanf(maxcut_in,"%d %d",&V_num,&E_num);
        vertex_num = V_num; 
        edge_num = E_num;
        in_edge = (Edges *)malloc(sizeof(Edges)*E_num);
        for(i=0;i<E_num;i++)
        {
           warn_2 = fscanf(maxcut_in,"%d %d %d",&(in_edge[i].v1_num),&(in_edge[i].v2_num),&(in_edge[i].value));
        }
    }
    else
    {
        printf("File Error\n");
        return ;
    }
    set_vertex_neighbor();
    fclose(maxcut_in);
}

void initialize_gene(int pop_num)
{
    int i,j;
    int cnt_1=0,cnt_0=0;
    gene_sol = (genetic_solution *)malloc(sizeof(genetic_solution)*pop_num);
    score_list = (int *)malloc(sizeof(int)*pop_num);
    fit = (float *)malloc(sizeof(float)*pop_num);
    index_sort = (int *)malloc(sizeof(int)*pop_num);
    for(i=0;i<pop_num;i++)
    {
        gene_sol[i].genetic = (int *)malloc(sizeof(int)*vertex_num);
        for(j=0;j<vertex_num;j++)
        {
            gene_sol[i].genetic[j] = randomize(2);
            /* equal 0 : 1 rate
            if(j>vertex_num * 0.3)
            {
                if(cnt_1 > cnt_0)
                    gene_sol[i].genetic[j] = 0;
                else gene_sol[i].genetic[j] = 1;

            }
            else 
            {
                // basic initialize
                gene_sol[i].genetic[j] = randomize(2);
            }
            if(gene_sol[i].genetic[j] == 1)
                cnt_1++;
            else cnt_0++;
          */
        }
        gene_sol[i].score = calculate_maxcut_score(&gene_sol[i]);
        score_list[i] = gene_sol[i].score;
    }
}

void set_vertex_neighbor() 
{
    int i;
    int *index;
    int e1_num;
    int e2_num;
    ver_n = (Vertex*)malloc(sizeof(Vertex)*vertex_num);
    index = (int *)malloc(sizeof(int)*vertex_num);
    
    for(i=0;i<vertex_num;i++)
    {
        ver_n[i].num =0;
        index[i] =0;
    }
    
    for(i=0;i<edge_num;i++)
    {
        ver_n[in_edge[i].v1_num-1].num += 1;
        ver_n[in_edge[i].v2_num-1].num += 1;
    }
    
    for(i=0;i<vertex_num;i++)
    {
        ver_n[i].n_vertex = (Neighbor*)malloc(sizeof(Neighbor)*ver_n[i].num);
    }
    
    for(i=0;i<edge_num;i++)
    {
        e1_num = in_edge[i].v1_num-1;
        e2_num = in_edge[i].v2_num-1;
        ver_n[e1_num].n_vertex[index[e1_num]].num = e2_num;
        ver_n[e2_num].n_vertex[index[e2_num]].num = e1_num;
        ver_n[e1_num].n_vertex[index[e1_num]].weight = in_edge[i].value;
        ver_n[e2_num].n_vertex[index[e2_num]].weight = in_edge[i].value;
        index[e1_num] += 1;
        index[e2_num] += 1;
    }

    free(index);

}
int calculate_maxcut_score(genetic_solution *solution)
{
    int i;
    int maxcut_score = 0;
    for(i=0;i<edge_num;i++)
    {
        if(solution->genetic[in_edge[i].v1_num-1] != solution->genetic[in_edge[i].v2_num-1]) 
            maxcut_score += in_edge[i].value;
    }
    return maxcut_score;
}

void array_index_sort()
{
    int i,j;
    int temp;
    for(i=0;i<population_num;i++)
    {
        index_sort[i] = i;
    }
    for(i=0;i<population_num;i++)
    {
        for(j=i+1;j<population_num;j++)
        {
            if(score_list[index_sort[i]] > score_list[index_sort[j]])
            {
                temp = index_sort[i];
                index_sort[i] = index_sort[j];
                index_sort[j] = temp;
            }
        }
    }       
}

void array_sort()
{
    int i,j;
    int temp;
    for(i=0;i<cut_num;i++)
    {
        for(j=i+1;j<cut_num;j++)
        {
            if(cut_points[i] > cut_points[j])
            {
                temp = cut_points[i];
                cut_points[i] = cut_points[j];
                cut_points[j] = temp;
            }
        }
    }
}

int randomize(int bit)
{
    int rand_num;
    rand_num = rand() % bit;
    return rand_num;
}

void maxcut_solver()
{
    int i = 0;
    int max_score,max_score_2;
    while((float)(clock() - start_time)/(CLOCKS_PER_SEC)<exec_time)
    {
        selection(select_mode);
        cross_over(cross_mode);
        mutation();
        /* for hw2 */
        local_optimization();
        //local_optimization(&selection_2);
        replacement(replacement_mode);
        max_score = max_score_maxcut();
        if(loop_num %10000 == 0)
        {
           file_write(output_file_name);
           //printf("Check process %d: %d\n",loop_num,max_score);
           //max_score_2 = calculate_maxcut_score(&optimal_sol);
            //if(max_score != max_score_2 )
                //printf("max_score error %d %d %d\n",max_score,max_score_2,optimal_sol.score);
        }
        loop_num += 1;
    }
    //check max score
    //printf("max score is %d\n",max_score);
}

void print_gene(genetic_solution *solution)
{
    int i;
    for(i=0;i<vertex_num;i++)
        printf("%d",solution->genetic[i]);
    printf("\n");
    return ;
}
void selection(int selection_mode)
{
    int i,j;
    float sum_fit = 0;
    float sum;
    float point;
    
    for(i=0;i<population_num;i++)
        score_list[i] = gene_sol[i].score;
    array_index_sort();
    min_index = index_sort[0];
    min2_index = index_sort[1];
    max_index = index_sort[population_num-1];
    max2_index = index_sort[population_num-2];
    if(selection_mode == 0) //random
    {
        sol1 = randomize(population_num);
        sol2 = randomize(population_num);
    }
    else if(selection_mode == 1) // select max and max-second score
    {
        sol1 = max_index;
        sol2 = max2_index;
    }
    else if(selection_mode == 2) //roulette wheel method
    {
        for(i=0;i<population_num;i++)
        {
            fit[i] = (score_list[i] - score_list[min_index]) + (score_list[max_index] - score_list[min_index]) / 4;
            sum_fit += fit[i];
        }
        for(j=0;j<2;j++)
        {
            sum = 0;
            point = ((float)rand()/(float)RAND_MAX) *sum_fit;
            for(i=0;i<population_num;i++)
            {
                sum += fit[i];
                if(point < sum)
                {
                    if(j==0)
                    {
                        sol1 = i;
                        break;
                    }
                    else if(j == 1)
                    {
                        sol2 = i;
                        break;
                    }
                }
            }
            if(sol1 == sol2)
            {
                j -= 1;
            }
        }
    }
    selection_1 = gene_sol[sol1];
    selection_2 = gene_sol[sol2];
    score_1 = gene_sol[sol1].score;
    score_2 = gene_sol[sol2].score;
    genetic_1 = (int *)malloc(sizeof(int)*vertex_num);
    genetic_2 = (int *)malloc(sizeof(int)*vertex_num);
    for(i=0;i<vertex_num;i++)
    {
        genetic_1[i] = gene_sol[sol1].genetic[i];
        genetic_2[i] = gene_sol[sol2].genetic[i];
    }
    /* To Check 
    printf("-----------------------Selection--------------\n");
    print_gene(&selection_1);
    print_gene(&selection_2);*/
}

void cross_over(int crossover_mode)
{
    int i,j;
    int temp;
    float rand_prob;
    if(crossover_mode == 1) //k-point-crossover 
    { 
        if(vertex_num < 10)
           cut_num = 2;
        else if(vertex_num <= 30)
            cut_num = 3;
        else if(vertex_num <= 100)
            cut_num = 4;
        else
            cut_num = 5;
            
        cut_points = (int *)malloc(sizeof(int)*(cut_num+1));
        cut_points[0] = 0;
        cut_points[cut_num] = vertex_num-1;
        for(i=1;i<cut_num;i++)
        {
            while(1)
            {
                cut_points[i] = randomize(vertex_num-1);
                if(cut_points[i] != 0 && cut_points[i] != vertex_num-1)
                {  
                   if(check_cutpoint(i) == 1)
                       break;
                }
            }
        }

        array_sort();
        for(i=0;i<cut_num;i++)
        {
            if(i%2 == 0)
            {
                for(j=cut_points[i];j<cut_points[i+1];j++)
                {
                    temp = selection_1.genetic[j];
                    selection_1.genetic[j] = selection_2.genetic[j];
                    selection_2.genetic[j] = temp;

                }
            }
        }
        

    }
    else if(crossover_mode == 2) // uniform-crossover
    {
        if(selection_1.score > selection_2.score)
            uniform_threshold = 0.55;
        else if(selection_1.score == selection_2.score)
            uniform_threshold = 0.5;
        else uniform_threshold = 0.45;

        for(i=0;i<population_num;i++)
        {
            rand_prob = (float)(rand() % 100) / (float)100;
            if(uniform_threshold < rand_prob)
            {
                temp = selection_1.genetic[i];
                selection_1.genetic[i] = selection_2.genetic[i];
                selection_2.genetic[i] = temp;
            }
        }
    }
/* To Check 
    printf("----------------------------crossover----------------\n");
    print_gene(&selection_1);
    print_gene(&selection_2);*/

}

int check_cutpoint(int num)
{
    int i;
    for(i=0;i<num;i++)
    {
        if(cut_points[num] == cut_points[i])
            return -1;
    }
    return 1;
}

void mutation()
{
    float rand_prob,rand_prob2;
    int i;

    for(i=0;i<vertex_num;i++) 
    {
        rand_prob = (float)rand() / (float)RAND_MAX;
        rand_prob2 = (float)rand() / (float)RAND_MAX;
        if(mutation_threshold > rand_prob)
        {
            selection_1.genetic[i] = 1 - selection_1.genetic[i];
        }
        if(mutation_threshold > rand_prob2)
        {
           selection_2.genetic[i] = 1 - selection_2.genetic[i];
        }
    }
/* To Check 
    printf("-------------------mutation----------------\n");
    print_gene(&selection_1);
    print_gene(&selection_2); */

}

void replacement(int replace_mode)
{
   
    selection_1.score = calculate_maxcut_score(&selection_1);
    selection_2.score = calculate_maxcut_score(&selection_2);
    gene_sol[sol1].score = selection_1.score;
    gene_sol[sol2].score = selection_2.score;
    if(replace_mode == 0)  //exchange parent & child
    {
        if(selection_1.score < score_1 )
        {
            gene_sol[sol1].genetic = genetic_1;
            gene_sol[sol1].score = score_1;
        }
        if(selection_2.score < score_2 )
        {
            gene_sol[sol2].genetic = genetic_2;
            gene_sol[sol2].score = score_2;
        }
    }
    else if(replace_mode == 1) //exchange lowest score gene
    {
        if(selection_1.score < score_1 )
        {
            gene_sol[min_index].genetic = genetic_1;
            gene_sol[min_index].score = score_1;
        }
        if(selection_2.score < score_2 )
        {
            gene_sol[min2_index].genetic = genetic_2;
            gene_sol[min2_index].score = score_2;
        }
    }
/* To Check 
    printf("------------------replacement--------------------\n");
    print_gene_2(genetic_1);
    print_gene_2(genetic_2);
    print_gene(&gene_sol[sol1]);
    print_gene(&gene_sol[sol2]);
*/
}

void print_gene_2(int *genetic)
{
    int i;
    for(i=0;i<vertex_num;i++)
        printf("%d",genetic[i]);
    printf("\n");
}
int max_score_maxcut()
{
    int i;
    int max_score = gene_sol[0].score;
    optimal_sol = gene_sol[0];
    for(i=0;i<population_num;i++)
    {
        if(max_score < gene_sol[i].score) 
        {
            max_score = gene_sol[i].score;
            optimal_sol = gene_sol[i];
        }
    }
    return max_score;
}
void file_write(char *file_name)
{
    FILE *maxcut_out;
    int i = 0;
    maxcut_out = fopen(file_name,"w");
    for(i=0;i<vertex_num;i++)
    {
        if(optimal_sol.genetic[i] == 1)
        {
            fprintf(maxcut_out,"%d ",i+1);
        }
    }
    fclose(maxcut_out);
}

void local_optimization() {
    int improved; 
    int i,j,k=0;
    int *idx;
    int temp,temp2;
    int df,neighbor_num,neighbor_idx,neighbor_weight;
    
    idx = (int *)malloc(sizeof(int)*vertex_num);
    for(i=0;i<vertex_num;i++)
        idx[i] = i;
     
    for(i=vertex_num-1;i>=0;i--)
    {
        temp = randomize(vertex_num);
        temp2 = idx[i];
        idx[i] = idx[temp];
        idx[temp] = temp2;
    }
    improved = 1;
    while(improved)
    {
        improved = 0;
        for(i=0;i<vertex_num;i++)
        {
           df = 0;
           neighbor_num = ver_n[idx[i]].num;
           for(j=0;j<neighbor_num;j++)
           {
                neighbor_idx = ver_n[idx[i]].n_vertex[j].num;
                neighbor_weight = ver_n[idx[i]].n_vertex[j].weight;
                if(selection_1.genetic[idx[i]] == selection_1.genetic[neighbor_idx])
                    df += neighbor_weight;
                else 
                    df -= neighbor_weight;
           }
           if(df > 0)
           {
               selection_1.genetic[idx[i]] = 1 - selection_1.genetic[idx[i]];
               improved = 1;
           }
        }
      
    }
    free(idx);
}

