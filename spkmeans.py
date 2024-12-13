import numpy as np
import pandas as pd
import sys
from mykmeanssp import spk, wam, ddg, gl, jacobi

error_msg = "An Error Has Occured"
np.random.seed(0)


def euclidean_dist(point1, point2):
    """calculate the euclidean distance"""
    return np.sqrt(np.sum(np.square(point1 - point2)))

    

def centroid_initialization(K, datapoints):
    """calculate the first k centroids"""
    N, dim = datapoints.shape
    centroids =np.zeros((K, dim))
    idx = np.random.choice(N)
    centroids[0, :] = datapoints[idx, :]
    indices = np.array(idx)

    distance = np.zeros(N)
    for n in range(N):
        distance[n] = euclidean_dist(centroids[0, :], datapoints[n, :])

    for cluster in range(1, K):

        idx = np.random.choice(N, p=(distance / np.sum(distance)))
        while idx in indices:
            idx = np.random.choice(N, p=(distance / np.sum(distance)))
        indices = np.append(indices, idx)
        new_centroid = datapoints[idx, :]
        for n in range(N):
            dist = euclidean_dist(new_centroid, datapoints[n, :])
            if dist < distance[n]: distance[n] = dist
        centroids[cluster, :] = new_centroid
    
    return (centroids, indices)

def sortValues(eigenvalues, eigenvectors):
    """sorts eigenvalues in ascending order and the eigenvectors by column"""
    indices_sort = eigenvalues.argsort()
    eigenvalues_sorted = eigenvalues[indices_sort]
    eigenvectors_sorted = eigenvectors[:,indices_sort]
    return eigenvalues_sorted , eigenvectors_sorted

def heuristic(eigenvalues):
    """calculate the number of clusters using Eigengap Heuristic"""
    diff = np.diff(eigenvalues)
    n = len(eigenvalues)
    k = np.argmax(diff[:n//2])
    return k+1 #add 1 because indexing starts at zero

def printMatrix(matrix):
    for row in matrix:
        pstr = ""
        for col in row:
            pstr = pstr + "{num:0.4f},".format(num=col)
        print(pstr[:-1])

def printArray(arr):
    idx_str = ""
    for idx in arr:
            idx_str = idx_str + "{},".format(idx)
    print(idx_str[:-1])

if __name__ == '__main__':
    '''
        process command-line arguments
        k (optional), goal, file_name
    '''

# process arguments from command line

# too few or too many arguments
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print(error_msg)
        exit()

# setting k
    if len(sys.argv) == 4:
        k = int(sys.argv[1])
    else:
        k = -1

# read data from file
    df = pd.read_csv(sys.argv[-1], header=None)
    datapoints = df.values
    n, dim = datapoints.shape

# perform computation according to "goal"

    goal = sys.argv[-2]
    
    if goal == 'spk': 
        #calculate the symetric graph laplacian mattrix from datapoints
        symMatrix = gl(datapoints.tolist(), n, dim)
        #calculate eigenvalues and eigenvectors and sort them
        jacobi_values = jacobi(symMatrix, n)
        eigenvalues = jacobi_values[0]
        eigenvectors = jacobi_values[1:]
        eigenvalues, eigenvectors = sortValues(np.array(eigenvalues),np.array(eigenvectors))
        #if k is not provided calsulate k
        if k == -1:
            k = heuristic(eigenvalues)
        #calculating initial centroids of the rows of the first k eigenvectors 
        init_centroids, indices = centroid_initialization(k, eigenvectors[:,:k])
        final_result = spk(eigenvectors[:,:k].tolist(), init_centroids.tolist(), n, k)
        printArray(indices)
        printMatrix(final_result)
    elif goal == 'wam': 
        final_result = wam(datapoints.tolist(), n, dim)
        printMatrix(final_result)
    elif goal == 'ddg': 
        final_result = ddg(datapoints.tolist(), n, dim)
        printMatrix(final_result)
    elif goal == 'gl': 
        final_result = gl(datapoints.tolist(), n, dim)
        printMatrix(final_result)
    elif goal == 'jacobi': 
        final_result = jacobi(datapoints.tolist(), n)
        printMatrix(final_result)
    else:
        print(error_msg)
        exit()

    