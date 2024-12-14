# Spectral Clustering Project - Unnormalized Spectral Clustering

This project implements the **Unnormalized Spectral Clustering Algorithm** as part of the **Software Project (0368-2161)** course. The implementation involves creating both Python and C modules for efficient computation and integration.

## 📋 Project Overview

Spectral clustering is a popular method for clustering data using the eigenvalues and eigenvectors of graph Laplacians. This project implements the following steps:
1. Create the Weighted Adjacency Matrix.
2. Compute the Graph Laplacian.
3. Calculate eigenvalues and eigenvectors using the Jacobi method.
4. Cluster the data using K-means++ initialized with the eigenvectors.
5. Determine the optimal number of clusters using the eigengap heuristic.

## 🚀 Features

This project includes:
- A Python interface for spectral clustering (`spkmeans.py`).
- A C implementation for computationally intensive tasks (`spkmeans.c`).
- A Python C API for seamless integration (`spkmeansmodule.c`).
- Support for various goals:
  - **wam**: Compute the Weighted Adjacency Matrix.
  - **ddg**: Compute the Diagonal Degree Matrix.
  - **gl**: Compute the Graph Laplacian.
  - **jacobi**: Calculate eigenvalues and eigenvectors.
  - **spk**: Perform full spectral clustering.

## 📂 Repository Structure

- `spkmeans.py`: Python interface.
- `spkmeans.c`: Core C implementation.
- `spkmeansmodule.c`: Python C API wrapper.
- `spkmeans.h`: Header file for the C implementation.
- `setup.py`: Build script for the Python C extension.
- `Makefile`: Script to compile the C code into an executable.
- Additional `.c` and `.h` files as necessary for modular design.

