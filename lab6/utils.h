#include <stdlib.h>
#include <stdio.h>

double **createMatrix(const int sizeA, const int sizeB){
    double **matrix = (double **) malloc(sizeA * sizeof(double *));
    for (int i = 0; i < sizeA; ++i) {
        matrix[i] = (double *) malloc(sizeB * sizeof(double));
    }
    return matrix;
}

void printMatrix(double **matrix, const int sizeA, const int sizeB) {
    for (int i = 0; i < sizeA; i++) {
        for (int j = 0; j < sizeB; j++) {
            printf("%8.2f", matrix[i][j]);
        }
        printf("\n");
    }
}

void freeMatrix(double **matrix) {
    free(matrix);
}

bool isCorrect(double **src, double **matrix, const int sizeA, const int sizeB){
    for (int i = 0; i < sizeA; ++i) {
        for (int j = 0; j < sizeB; ++j) {
            if(src[i][j] != matrix[i][j]){
                return false;
            }
        }
    }
    return true;
}