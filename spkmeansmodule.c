# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include "spkmeans.h"
#include "spkalgo.h"

#define MEMERROR(X)                                                     \
    if (X == NULL)                                                      \
    {                                                                   \
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occured\n");   \
        return NULL;                                                    \
    }                                                                   \


/*HELPER FUNCTIONS*/

double** convertPy2C(PyObject *datapoints_py,int n,int dim){
    int i,j;
    PyObject *coord_py, *point;
    double coord_c;
    double** datapoints_c;
    
    /*Allocating memory for datapoints*/
    datapoints_c = malloc(n* sizeof(double *));
    MEMERROR(datapoints_c);
    for(i=0; i<n; i++){
        datapoints_c[i] = malloc(dim*sizeof(double));
        MEMERROR(datapoints_c[i]);
    }

    /*creating the datapoints as a c array*/
    for(i = 0; i<n; i++){
        point = PyList_GetItem(datapoints_py, i);
        
        for(j=0; j< dim; j++){
            coord_py = PyList_GetItem(point,j);
            coord_c = PyFloat_AsDouble(coord_py);
            datapoints_c[i][j] = coord_c;

        }
    }

    return datapoints_c;
}

PyObject *convertC2Py(double** final_c,int n, int dim){
    /*convert c array into a PyObject List*/
    int i,j;
    PyObject *final_py;
    
    final_py = PyList_New(n);
    for(i=0; i<n; i++){
        PyList_SetItem(final_py, i, PyList_New(dim));
        for(j=0; j< dim; j++){
            PyList_SetItem(PyList_GetItem(final_py,i),j, PyFloat_FromDouble(final_c[i][j]));
        }
    }

    return final_py;

}



/*WRAPING FUNCTIONS*/

/*wraping wam function*/
static PyObject *wam_wrap(PyObject *self, PyObject *args){
    int i, n, dim;
    PyObject *datapoints_py, *final_py;
    double **datapoints_c, **final_c;


    /* This parses the Python arguments into a PyObject and two integers*/
    if(!PyArg_ParseTuple(args, "Oii", &datapoints_py, &n, &dim)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*Allocating memory for datapoints_c and filling its values*/
    datapoints_c = convertPy2C(datapoints_py,n,dim);

    /*calling on wam_calc function*/
    final_c = wam_calc(datapoints_c,n,dim);
    MEMERROR(final_c);

    /*Free dataponts memory*/
    for(i=0; i<n; i++){
        free(datapoints_c[i]);
    }
    free(datapoints_c);

    /*convert c array into a PyObject List*/
    final_py = convertC2Py(final_c, n, n);

    for (i = 0; i < n; i++) {
        free(final_c[i]);
    }
    free(final_c);


    return final_py;

}


/*wraping ddg function*/
static PyObject *ddg_wrap(PyObject *self, PyObject *args){
    int i, n, dim;
    PyObject *datapoints_py, *final_py;
    double **datapoints_c, **final_c;


    /* This parses the Python arguments into a PyObject and two integers*/
    if(!PyArg_ParseTuple(args, "Oii", &datapoints_py, &n, &dim)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*Allocating memory for datapoints_c and filling its values*/
    datapoints_c = convertPy2C(datapoints_py,n,dim);

    /*calling on ddg_calc function*/
    final_c = ddg_calc(datapoints_c,n,dim);
    MEMERROR(final_c);

    /*Free dataponts memory*/
    for(i=0; i<n; i++){
        free(datapoints_c[i]);
    }
    free(datapoints_c);

    /*convert c array into a PyObject List*/
    final_py = convertC2Py(final_c, n, n);

    for (i = 0; i < n; i++) {
        free(final_c[i]);
    }
    free(final_c);


    return final_py;

}


/*wraping gl function*/
static PyObject *gl_wrap(PyObject *self, PyObject *args){
    int i, n, dim;
    PyObject *datapoints_py, *final_py;
    double **datapoints_c, **final_c;


    /* This parses the Python arguments into a PyObject and two integers*/
    if(!PyArg_ParseTuple(args, "Oii", &datapoints_py, &n, &dim)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*Allocating memory for datapoints_c and filling its values*/
    datapoints_c = convertPy2C(datapoints_py,n,dim);

    /*calling on gl_calc function*/
    final_c = gl_calc(datapoints_c,n,dim);
    MEMERROR(final_c);

    /*Free dataponts memory*/
    for(i=0; i<n; i++){
        free(datapoints_c[i]);
    }
    free(datapoints_c);

    /*convert c array into a PyObject List*/
    final_py = convertC2Py(final_c, n, n);

    for (i = 0; i < n; i++) {
        free(final_c[i]);
    }
    free(final_c);


    return final_py;

}


/*wraping jacobi function*/
static PyObject *jacobi_wrap(PyObject *self, PyObject *args){
    int i, n;
    PyObject *datapoints_py, *final_py;
    double **datapoints_c, **final_c;


    /* This parses the Python arguments into a PyObject and two integers*/
    if(!PyArg_ParseTuple(args, "Oi", &datapoints_py, &n)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*Allocating memory for datapoints_c and filling its values*/
    datapoints_c = convertPy2C(datapoints_py,n,n);

    /*calling on jacobi_calc function*/
    final_c = jacobi_calc(datapoints_c,n);
    MEMERROR(final_c);

    /*Free dataponts memory*/
    for(i=0; i<n; i++){
        free(datapoints_c[i]);
    }
    free(datapoints_c);

    /*convert c array into a PyObject List*/
    final_py = convertC2Py(final_c, n+1, n);

    for (i = 0; i < n+1; i++) {
        free(final_c[i]);
    }
    free(final_c);


    return final_py;

}


/*wraping the spk function*/
static PyObject *spk_wrap(PyObject *self, PyObject *args)
{
    int i, n ,k;
    PyObject *centroids_py,*datapoints_py, *final_py;
    double **datapoints_c, **centroids_c, **final_c;


    /* This parses the Python arguments into a double (d)  variable named eps and int (i) variable named iter and two PyObjects*/
    if(!PyArg_ParseTuple(args, "OOii", &datapoints_py, &centroids_py, &n, &k)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*Allocating memory for datapoints_c and filling its values*/
    datapoints_c = convertPy2C(datapoints_py,n,k);
    /*Allocating memory for centroids_c and filling its values*/
    centroids_c = convertPy2C(centroids_py,k,k);


    /*call spk_calc function*/
    final_c = spk_calc(datapoints_c, centroids_c, n, k);
    MEMERROR(final_c);

    /*Free datapoints and centroids memory*/
    for(i=0; i<n; i++){
        free(datapoints_c[i]);
    }

    for ( i = 0; i < k; i++){
        free(centroids_c[i]);
    }

    /*convert c array into a PyObject List*/
    final_py = convertC2Py(final_c, k, k);

    /*free final_c*/

    for (i = 0; i < k; i++) {
        free(final_c[i]);
    }

    free(datapoints_c);
    free(centroids_c);
    free(final_c);


    return final_py;
}

/*METHOD FUNCTION*/
static PyMethodDef spkmeansMethods[] = {
    {"wam", (PyCFunction) wam_wrap, METH_VARARGS, PyDoc_STR("calculating the weighted adjacency matrix")},
    {"ddg", (PyCFunction) ddg_wrap, METH_VARARGS, PyDoc_STR("calculating the diagonal degree matrix")},
    {"gl", (PyCFunction) gl_wrap, METH_VARARGS, PyDoc_STR("calculating the graph laplacian matrix")},
    {"jacobi", (PyCFunction) jacobi_wrap, METH_VARARGS, PyDoc_STR("calculating the eigen values and eigen vectors using the jacobi algorithm")},
    {"spk", (PyCFunction) spk_wrap, METH_VARARGS, PyDoc_STR("performing the spk clustering algorithm")},
    {NULL, NULL, 0, NULL}    
};

/*initiating the Module*/

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */  
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    spkmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}