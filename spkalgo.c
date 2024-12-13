/*TODO do we need _GNU_SOURCE*/
#define  _GNU_SOURCE
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include "spkalgo.h"

#define MEMERROR(X)                                                 \
    if (X == NULL)                                                  \
    {                                                               \
        return NULL;                                                \
    }                                                               \

double const eps = 0;
int const iter = 300;
int k, N, dim = 0;
double **init_centroids;


struct cord
{
    double value;
    struct cord *next;
};

struct vector
{
    struct vector *next;
    struct cord *cords;
};



int min_dist(struct vector *curr_point, struct vector *head_vec){
    /*determens the closest cluster for curr_point*/
    
    int index = 0;
    struct vector *curr_clust;
    struct cord *head_cord, *curr_cord, *curr_point_cord;
    double sum = 0;
    double dist, p1, p2;
    double min_dist = INFINITY;
    int i;

    curr_clust = head_vec;
    head_cord = curr_point->cords;
    curr_cord = curr_clust->cords;
    curr_point_cord = head_cord;

    for(i = 0; i < k; i++){

        /*calculating dist of point and current cluster*/
        
        while(curr_cord != NULL){
            p1 = (curr_cord->value);
            p2 = (curr_point_cord->value);
            sum += (p1 - p2)*(p1-p2);
            curr_cord = curr_cord->next;
            curr_point_cord = curr_point_cord->next;
        }

        /*check if dist is smaller than min_dist*/
        dist = sqrt(sum);
        if(dist < min_dist){
            index = i;
            min_dist = dist;
        }

        curr_point_cord = head_cord;
        curr_clust = curr_clust->next;
        curr_cord = curr_clust->cords;
        sum = 0;
    }

    return index;
}

double** kmeans(struct vector *head)
{
    /*performs kmeans clustering algorithm according to HW2 description*/

    int *clust_assignment;
    double **coord_sum, **final_centroids;
    struct vector *head_vec, *curr_vec, *curr_point, *prev_vec;
    struct cord *head_cord, *curr_cord, *prev_cord;
    int index;
    int counter = 0;
    double sum = 0, tmp = 0, dist = 0;
    int terminate = 1;
    int i, j, c, m, n = 0;


    clust_assignment = calloc(k, sizeof(int));
    MEMERROR(clust_assignment)

    coord_sum = (double **) calloc(k,sizeof(double*));
    MEMERROR(coord_sum)
    for(i =0; i<k; i++){
        coord_sum[i] = (double *) calloc(dim, sizeof(double));
        MEMERROR(coord_sum[i])
    } 

    /*create a linked list for the initial centroids*/
    head_cord = malloc(sizeof(struct cord));
    MEMERROR(head_cord)
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector));
    MEMERROR(head_vec)
    curr_vec = head_vec;
    curr_vec->next = NULL;

    for(i=0; i<k; i++){
        
        /*adding point coordinates*/
        for(j=0; j<dim; j++){
            curr_cord->value = init_centroids[i][j];
            curr_cord->next = malloc(sizeof(struct cord));
            MEMERROR(curr_cord->next)
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }
        /*creating a new point*/
        curr_vec->cords = head_cord;
        curr_vec->next = malloc(sizeof(struct vector));
        MEMERROR(curr_vec->next)
        curr_vec = curr_vec->next;
        curr_vec->next = NULL;
        head_cord = malloc(sizeof(struct cord));
        MEMERROR(head_cord)
        curr_cord = head_cord;
        curr_cord->next = NULL;
    }

    /*set curr_point to beginning of datapoints linked list*/
    curr_point = head;
    for(i = 0; i < iter; i++){

        /*assign every point to the closest cluster*/
        while(curr_point->next != NULL){
            index = min_dist(curr_point, head_vec);
            /*count how many points are in each cluster*/
            clust_assignment[index]++;
            curr_cord = curr_point->cords;
            while(curr_cord != NULL){
                coord_sum[index][counter] += curr_cord->value;
                curr_cord = curr_cord->next;
                counter++;
            }
           
            counter = 0;
            curr_point = curr_point->next;
           
        }

        curr_point = head;

        /*update the centroids*/
        curr_vec = head_vec;
        for(j = 0; j < k; j++){

            curr_cord = curr_vec->cords;
            
            /*assing new vealue for each coordinate*/
            for (c = 0; c < dim; c++)
            {
                tmp = coord_sum[j][c] / clust_assignment[j];
                sum += pow((tmp - curr_cord->value),2);
                curr_cord->value = tmp;
                curr_cord = curr_cord->next;
            }
            
            /*calculate dist to previous centroid*/
            dist = sqrt(sum);
            if(dist > eps){
                terminate = 0;
            }

            sum = 0;
            curr_vec = curr_vec->next;

        }

        /*check for convergence*/
        if(terminate){
            break;
        }
        terminate = 1;

        /*reset arrays to zero*/
        for (m = 0; m < k; m++)
        {
            clust_assignment[m] = 0;
            for (n = 0; n < dim; n++)
            {
                coord_sum[m][n] = 0;
            }
        }
    }


    /*convert the linked list of centroids to a 2D array*/
    final_centroids = malloc(k*sizeof(double*));
    MEMERROR(final_centroids);
    for ( i = 0; i < k; i++){
        final_centroids[i] = malloc(dim*sizeof(double));
        MEMERROR(final_centroids[i]);
        }
    
    curr_vec = head_vec;

    for (i = 0; i < k; i++)
    {
        j=0;
        curr_cord = curr_vec->cords;
        while (curr_cord->next != NULL)
        {
            final_centroids[i][j] = curr_cord->value;
            curr_cord = curr_cord->next;
            j++;
        }
        curr_vec = curr_vec->next;

    }



    /*Free Centroids*/
    curr_vec = head_vec;
    while (curr_vec->next != NULL){
        curr_cord = curr_vec->cords;
        while (curr_cord->next != NULL){
            prev_cord = curr_cord;
            curr_cord = curr_cord->next;
            free(prev_cord);
        }
        free(curr_cord);
        prev_vec = curr_vec;
        curr_vec = curr_vec->next;
        free(prev_vec);
        }
    free(curr_vec);
    /*free(head_vec);*/

    for (m = 0; m < k; m++)
    {
        free(coord_sum[m]);
    }

    free(coord_sum);
    free(clust_assignment);
    free(head_cord);

    /*Return the 2D array of centroids*/
    return final_centroids;
}

double** spk_calc(double **datapoints, double **initial_centroids, int n, int K){
    /*
        Converts data to linked lists and call on kmeans method
        returns final centroids
    */
    int i,j;
    double **final_centroids;

    struct vector *head_vec, *curr_vec, *prev_vec;
    struct cord *head_cord, *curr_cord, *prev_cord;

    k = K;
    dim = K;
    N = n;
    init_centroids = initial_centroids;
    
    /*create a linked list for the datapoints*/
    head_cord = malloc(sizeof(struct cord));
    MEMERROR(head_cord)
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector));
    MEMERROR(head_vec)
    curr_vec = head_vec;
    curr_vec->next = NULL;

    for(i=0; i<N; i++){
        
        /*adding point coordinates*/
        for(j=0; j<dim; j++){
            curr_cord->value = datapoints[i][j];
            curr_cord->next = malloc(sizeof(struct cord));
            MEMERROR(curr_cord->next)
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }
        /*creating a new point*/
        curr_vec->cords = head_cord;
        curr_vec->next = malloc(sizeof(struct vector));
        MEMERROR(curr_vec->next)
        curr_vec = curr_vec->next;
        curr_vec->next = NULL;
        head_cord = malloc(sizeof(struct cord));
        MEMERROR(head_cord)
        curr_cord = head_cord;
        curr_cord->next = NULL;


    }

    /*call kmeans*/
    final_centroids = kmeans(head_vec);

    /*freeing the datapoints*/
    /*FREE MEMORY*/
    curr_vec = head_vec;
    while (curr_vec->next != NULL)
    {
        curr_cord = curr_vec->cords;
        while (curr_cord->next != NULL)
        {
            prev_cord = curr_cord;
            curr_cord = curr_cord->next;
            free(prev_cord);
        }
        free(curr_cord);
        
        prev_vec = curr_vec;
        curr_vec = curr_vec->next;
        free(prev_vec);
    }
    
    free(curr_vec);
    free(head_cord);


    return final_centroids;
}