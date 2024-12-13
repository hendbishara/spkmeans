# ifndef SPKMEANS_H_
# define SPKMEANS_H_

double** wam_calc(double** datapoints, int n, int dim);
double** ddg_calc(double** datapoints, int n, int dim);
double** gl_calc(double** datapoints, int n, int dim);
double** jacobi_calc(double** matrix, int n);

# endif