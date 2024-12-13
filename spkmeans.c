#define  _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "spkmeans.h"



#define MEMERROR(X)                                                 \
    if (X == NULL)                                                  \
    {                                                               \
        printf("%s", "An Error Has Occured\n");                     \
        exit(EXIT_FAILURE);                                         \
    }                                                               \
    
#define MEMERRORFUNC(X)                                             \
    if (X == NULL)                                                  \
    {                                                               \
        return NULL;                                                \
    }                                                               \

double** wam_matrix;
double** ddg_matrix;
double** gl_matrix;
double** combined_values;
double const epsilon = 0.0;

/*HELPER FUNCTIONS*/

double SquaredED(double* a, double* b, int dim){
    /*Calculate the square euclidean distance according to the project description*/
    double sum = 0;
    int i;
    for(i=0; i<dim; i++){
        sum += pow((a[i]-b[i]),2);
    }
    return sum;
}

double rowSum(double* r, int length){
    /*calculate the sum of matrix row*/
    int i;
    double sum = 0;
    for(i=0; i<length; i++){
        sum += r[i];
    }
    return sum;
}

int* pivot(double** A, int n, int* ij){
    /*find index i and j with largest absolute off-diagonal value*/
    int i,j,max_i = 0, max_j = 0;
    double max = -1;
    for(i=0; i<n; i++){
        /*only look at the right side of the diagonal of matrix A because it's semetric*/
        for(j=i+1;j<n;j++){
            if(fabs(A[i][j]) > max){
                max = fabs(A[i][j]);
                max_i = i;
                max_j = j;
             }
        }
    }
    ij[0] = max_i;
    ij[1] = max_j;
    return ij;

}

double* rotationMatrix(double** A, int i, int j, double* sc){
    /*calculates the s and c values of the Jacobi rotation matrix*/
    double teta, t, c, s;

    /*calculating teta*/
    teta = (A[j][j]-A[i][i])/(2*A[i][j]);

    /*calculating t*/
    t = (teta < 0) ? -1 : 1;
    t = t/((teta*t)+ sqrt(pow(teta,2)+1));

    /*calculating c*/
    c = 1/(sqrt(pow(t,2)+1));

    /*calculating s*/
    s = t*c;

    sc[0] = s;
    sc[1] = c;

    return sc;
}

void transformA(double** A,double c,double s, int i, int j, int n){
    /*calculate Matrix A' according to project description*/
    int r;
    double a_ri, a_rj;
    double a_ij = A[i][j];
    double a_ii = A[i][i];
    double a_jj = A[j][j];
    for(r=0; r<n; r++){
        if((r!=i) & (r!=j)){
            a_ri = A[r][i];
            a_rj = A[r][j];
            A[r][i] = c*a_ri - s*a_rj;
            A[i][r] = c*a_ri - s*a_rj;
            A[r][j] = c*a_rj + s*a_ri;
            A[j][r] = c*a_rj + s*a_ri;

        }
            
    
    }
    A[i][i] = pow(c,2)*a_ii + pow(s,2)*a_jj - 2*s*c*a_ij;
    A[j][j] = pow(s,2)*a_ii + pow(c,2)*a_jj + 2*s*c*a_ij;
    A[i][j] = (pow(c,2) - pow(s,2))*a_ij + s*c*(a_ii - a_jj);
    A[j][i] = A[i][j];
}



double off(double** A, int n){
    /*calculating sum of squares of all off diagonal elements*/
    double sum = 0;
    int i,j;
    for( i = 0; i<n ; i++){
        for(j = 0; j<n; j++){
            if(i!=j){
                sum += pow(A[i][j],2);
            }
        }
    }

    return sum;
}

double** updateEGV(double** matrix_vals, int i, int j, double s, double c, int n){
    /*calculate the product of multiplying two rotation matrices*/
    int row,col,ind;
    double sum_i = 0, sum_j = 0;
    double *new_i_col, *new_j_col, *i_col, *j_col;
    /*create ith and jth columns in new rotation matrix
      and new i and j columns for the undated columns in the multiplied matrix*/
    new_i_col = malloc(n*sizeof(double));
    MEMERRORFUNC(new_i_col);
    new_j_col = malloc(n*sizeof(double));
    MEMERRORFUNC(new_j_col);
    i_col = malloc(n*sizeof(double));
    MEMERRORFUNC(i_col);
    j_col = malloc(n*sizeof(double));
    MEMERRORFUNC(j_col);
    for(ind = 0; ind<n; ind++){
        new_i_col[ind] = 0;
        new_j_col[ind] = 0;
        j_col[ind] = 0;
        i_col[ind] = 0;
    }
    i_col[i] = c;
    i_col[j] = -s;
    j_col[j] = c;
    j_col[i] = s;

    /*calculate the ith and jth column in eigenvectors matrix*/
    for(row = 1; row<n+1; row++){
        for(col = 0; col<n; col++){
            sum_i = sum_i+(matrix_vals[row][col] * i_col[col]);
            sum_j = sum_j+(matrix_vals[row][col]*j_col[col]);
        }
        new_i_col[row-1] = sum_i;
        new_j_col[row-1] = sum_j;
        sum_i = 0;
        sum_j = 0;
    }

    /*replace ith and jth columns in eignevectors matrix with new values*/
    for(row = 1; row<n+1; row++){
        matrix_vals[row][i] = new_i_col[row-1];
        matrix_vals[row][j] = new_j_col[row-1];
    }

    /*free memory used*/
    free(new_i_col);
    free(new_j_col);
    free(i_col);
    free(j_col);

    return matrix_vals;
}

void printMatrix(double** matrix, int rows,int cols){
    int i,j;
    for(i=0; i<rows; i++){
        for(j=0; j<cols-1; j++){
            printf("%.4f%c",matrix[i][j],',');
        }
        printf("%.4f\n",matrix[i][j]);
    }
}

/*END OF HELPER FUNCTIONS*/

double** wam_calc(double** datapoints, int n, int dim){
    /*
        datapoints as 2d array, n: number ofdatapoints
        calculate the weighted adjacency matrix
    */
    int i, j;
    double sed; /*squared euclidean distance*/

    /* allocate memory for 2d array of size nxn for WAM */
    wam_matrix = malloc(n* sizeof(double *));
    MEMERRORFUNC(wam_matrix);
    for(i=0; i<n; i++){
        wam_matrix[i] = malloc(n*sizeof(double));
        MEMERRORFUNC(wam_matrix[i]);
    }
    /*Iterating over data points to calculate WAM*/
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i!=j){
                sed = SquaredED(datapoints[i],datapoints[j],dim);
                wam_matrix[i][j] = exp(-(sed)/2);
            }

            else{
                wam_matrix[i][i] = 0;
            }
        }
    }

    return wam_matrix;
}

double** ddg_calc(double** datapoints, int n, int dim){
    /*
        datapoints as 2d array, n: number ofdatapoints
        calculate the diagonal degree matrix
    */

    int i,j;

    /* intitialize wam*/
    wam_matrix = wam_calc(datapoints, n, dim);
    /* allocate memory for 2d array of size nxn and fill the values for ddg */
    ddg_matrix = malloc(n* sizeof(double *));
    MEMERRORFUNC(ddg_matrix);
    for(i=0; i<n; i++){
        ddg_matrix[i] = malloc(n*sizeof(double));
        MEMERRORFUNC(ddg_matrix[i]);
        /*calculating ddg values*/
        for(j=0; j<n; j++){
            ddg_matrix[i][j] = (i!=j)? 0: rowSum(wam_matrix[i],n);
        }
    }

    return ddg_matrix;

}

double** gl_calc(double** datapoints, int n,int dim){
    /*
        datapoints as 2d array, n: number ofdatapoints
        calculate the graph laplacian
    */
    int i,j;

    /* initialize ddg*/
    ddg_matrix = ddg_calc(datapoints, n, dim);
    

    /* allocate memory for 2d array of size nxn for WAM */
    gl_matrix = malloc(n* sizeof(double *));
    MEMERRORFUNC(gl_matrix);
    for(i=0; i<n; i++){
        gl_matrix[i] = malloc(n*sizeof(double));
        MEMERRORFUNC(gl_matrix[i]);
        /*initializeing the values of gl*/
        for(j=0; j<n;j++){
            gl_matrix[i][j] = ddg_matrix[i][j]-wam_matrix[i][j];
        }
    }

    return gl_matrix;

}

double** jacobi_calc(double** matrix, int n){
    /*
        matrix as 2d array
        calculate the eigenvectors and eigenvalues
    */
    int cnt = 0,i = 0,j = 0,l = 0, ind =0;
    double off_A, off_TA;
    double* sc;
    int* ij;
    /*Allocate memorry for eigenvectors matrix with diagonal equals 1s*/
    combined_values = malloc((n+1)* sizeof(double *));
    MEMERRORFUNC(combined_values);
    for(i=0; i<n+1; i++){
        combined_values[i] = malloc(n*sizeof(double));
        MEMERRORFUNC(combined_values[i]);
        /*fill matrix with zeros and ones on diagonal*/
        for(j=0; j<n; j++){
            combined_values[i][j] = ((i-1)!=j)? 0 : 1;
        }
    }

    /*Allocate memory for ij and sc*/
    ij = malloc(2*sizeof(int));
    MEMERRORFUNC(ij);
    sc = malloc(2*sizeof(double));
    MEMERRORFUNC(sc);

    /*implement the iterative Jacobi algorithm*/
    while(cnt<100){
        off_A = off(matrix,n);
        ij = pivot(matrix,n,ij); /*choose pivot indices*/
        i = ij[0];
        j = ij[1];

        sc = rotationMatrix(matrix, i, j, sc); /*calculate c & s in Rotation Matrix P*/
        transformA(matrix, sc[1], sc[0],i,j, n); /*calculate A'=P^(T)AP*/
        /*multiply rotation matrices*/
        /*fill s and c values in the first iteration*/
        if(cnt == 0){
            combined_values[i+1][i] = sc[1];
            combined_values[j+1][j] = sc[1];
            combined_values[i+1][j] = sc[0];
            combined_values[j+1][i] = -sc[0];
        }
        
        else{
            combined_values = updateEGV(combined_values,i,j,sc[0],sc[1],n);
        }
        
        /*check for convergence*/
        off_TA = off(matrix,n);
        if(off_A-off_TA <= epsilon){
            break;
        }

        cnt++;
    }

    /*save eignevalues*/
    for(l=0;l<n;l++){
        if(combined_values[0][l] == 0){
            /*check for -0.000 eigenvalues, if found: multiply eigenvector with -1*/
            if(signbit(combined_values[0][l])){
                combined_values[0][l]=0.0;
                for(ind=1;ind<n+1;ind++){
                    combined_values[ind][l]= -1 * combined_values[ind][l];
                }
            }

            else{
                combined_values[0][l]=matrix[l][l];
            }
        }

        else{
            combined_values[0][l]=matrix[l][l];
        }
    }

    /*free memory allocated in helper function*/
    free(ij);
    free(sc);

    return combined_values;

}





int main(int argc, char **argv){

    int i, j,n = 0,dim = 0;
    size_t len = 0;
    char *line = NULL;
    char *head;
    char *cont;
    char *ptr;
    double **datapoints;
    double data;
    FILE *ifp = NULL;
    

    /* check for a valid number of CMD arguments */
    if (argc != 3){
        printf("An Error Has Occured\n");
        exit(EXIT_FAILURE);
    }

    /* parse file data */ 
    if((ifp = fopen(argv[2], "r"))!=NULL){
        
        errno = 0;
        
        /* go through file to determine dim and number of datapoints*/
        while(getline(&line, &len, ifp) != -1) {
            /*get the dim from the first line of the file*/
            if(n == 0){
                head = line;
                for (i=0; line[i]; line[i]==',' ? i++ : *line++);
                dim = i + 1;
                line = head;
            }
            /* jacobi matrix is n x n*/
            if(strcmp(argv[1], "jacobi") == 0){
                n = dim;
                break;              
            }
            n ++;

            /*reset values of line and len for memory reallocation by getline()*/
            free(line);
            line  = NULL;
            len = 0;
            
            
        }
        
        /*reset values of line and len for memory reallocation by getline()*/
        free(line);
        line  = NULL;
        len = 0;

        /* errno is 0 for EOF */
        if (errno != 0){
            printf("An Error Has Occured\n");
            exit(EXIT_FAILURE);
        }
        
        /* point back to start of the file */
        rewind(ifp);

        /* allocate memory for 2d array of the data */
        datapoints = malloc(n* sizeof(double *));
        MEMERROR(datapoints);
        for(i=0; i<n; i++){
            datapoints[i] = malloc(dim*sizeof(double));
            MEMERROR(datapoints[i]);
        }
        
        /* iterate over file and store data in datapoints, i as row index, j as column index */
        i = 0;
        while(getline(&line, &len, ifp) != -1) {
            
            j = 0; 
            data = strtod(line, &ptr);
            datapoints[i][j] = data;

            /* fill all coordinates of the i-th point*/
            while (*ptr != '\n')
            {
                j++;
                cont = ptr;
                cont++; /*discard comma after each pint*/
                data = strtod(cont, &ptr);
                datapoints[i][j] = data;
            }   
            i++;

            /*reset values of line and len for memory reallocation by getline()*/
            free(line);
            line  = NULL;
            len = 0;

        }


        /*reset values of line and len for memory reallocation by getline()*/
            free(line);
            line  = NULL;
            len = 0;

        /* errno is 0 for EOF */
        if (errno != 0){
            printf("An Error Has Occured\n");
            exit(EXIT_FAILURE);
        }
        /*free line pointer*/
        free(line);


        fclose(ifp);
    }
    else{
        printf("An Error Has Occured\n");
        exit(EXIT_FAILURE);
    }
    

    /* perform computation according to "goal" */
    if (strcmp(argv[1], "wam") == 0){
        wam_matrix = wam_calc(datapoints,n,dim);
        MEMERROR(wam_matrix);
        printMatrix(wam_matrix,n,n);
        
        /*free memory*/
        for(i=0; i<n; i++){
            free(datapoints[i]);
            free(wam_matrix[i]);
        }
        free(datapoints);
        free(wam_matrix);
    }
     

    else if (strcmp(argv[1], "ddg") == 0) {
        ddg_matrix = ddg_calc(datapoints,n,dim);
        MEMERROR(ddg_matrix);
        printMatrix(ddg_matrix,n,n);
        /*free memory*/
        for(i=0; i<n; i++){
            free(datapoints[i]);
            free(wam_matrix[i]);
            free(ddg_matrix[i]);
        }
        free(datapoints);
        free(wam_matrix);
        free(ddg_matrix);
    } else if (strcmp(argv[1], "gl") == 0) {
        gl_matrix = gl_calc(datapoints,n,dim);
        MEMERROR(gl_matrix);
        printMatrix(gl_matrix,n,n);
        /*free memory*/
        for(i=0; i<n; i++){
            free(datapoints[i]);
            free(wam_matrix[i]);
            free(ddg_matrix[i]);
            free(gl_matrix[i]);
        }
        free(datapoints);
        free(wam_matrix);
        free(ddg_matrix);
        free(gl_matrix);
    } else if (strcmp(argv[1], "jacobi") == 0) {
        combined_values =  jacobi_calc(datapoints,n);
        MEMERROR(combined_values);
        printMatrix(combined_values,n+1, n);
        printf("\n");
        /*free memory*/
        for(i=0; i<n; i++){
            free(datapoints[i]);
            free(combined_values[i]);
        }

        free(datapoints);
        free(combined_values[n]);
        free(combined_values);
    } else { 
        printf("An Error Has Occured\n");
        exit(EXIT_FAILURE);
    }

    return 0;

}