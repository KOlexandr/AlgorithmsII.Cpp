#include <iostream>
#include <conio.h>
#include <stdlib.h>

namespace winograd {
    void multiplySerial(double **result, double **first, double **second, const int size) {
        const int d = size / 2;
        double *rowFactor = (double *) malloc(size * sizeof(double));
        for (int i = 0; i < size; ++i) {
            rowFactor[i] = first[i][0] * first[i][1];
            for (int j = 1; j < d; ++j) {
                rowFactor[i] += first[i][2 * j - 1] * first[i][2 * j];
            }
        }

        double *columnFactor = (double *) malloc(size * sizeof(double));
        for (int i = 0; i < size; ++i) {
            columnFactor[i] = second[0][i] * second[1][i];
            for (int j = 1; j < d; ++j) {
                columnFactor[i] += second[2 * j - 1][i] * second[2 * j][i];
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result[i][j] = -rowFactor[i] - columnFactor[j];
                for (int k = 0; k < d; ++k) {
                    result[i][j] += (first[i][2 * k] + second[2 * k + 1][j]) * (first[i][2 * k + 1] + second[2 * k][j]);
                }
            }
        }
    }

    void multiplyParallel(double **result, double **first, double **second, const int size) {
        const int d = size / 2;
        double *rowFactor = (double *) malloc(size * sizeof(double));
        #pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            rowFactor[i] = first[i][0] * first[i][1];
            for (int j = 1; j < d; ++j) {
                rowFactor[i] += first[i][2 * j - 1] * first[i][2 * j];
            }
        }

        double *columnFactor = (double *) malloc(size * sizeof(double));
        #pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            columnFactor[i] = second[0][i] * second[1][i];
            for (int j = 1; j < d; ++j) {
                columnFactor[i] += second[2 * j - 1][i] * second[2 * j][i];
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                result[i][j] = -rowFactor[i] - columnFactor[j];
                for (int k = 0; k < d; ++k) {
                    result[i][j] += (first[i][2 * k] + second[2 * k + 1][j]) * (first[i][2 * k + 1] + second[2 * k][j]);
                }
            }
        }
    }
}