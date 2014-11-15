#include <stdio.h>
#include <omp.h>
#include "lab6/winogradMultiplication.cpp"
#include "lab6/recursiveMultiplication.cpp"
#include "lab6/recursiveMultiplicationInPlace.cpp"

#define sizesSize 15
const int sizes[sizesSize] = {50,80,100,200,300,400,500,600,700,800,900,1000,1500,2000,3000};

void testMultiplicationWithPrint(const int size, void (*multiply)(double **, double **, double **, const int)){
    double **first = createMatrix(size, size);
    double **second = createMatrix(size, size);
    double **result = createMatrix(size, size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            first[i][j] = i + 1;
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            second[i][j] = j + 1;
        }
    }

    printf("A\n");
    printMatrix(first, size, size);
    printf("B\n");
    printMatrix(second, size, size);

    multiply(result, first, second, size);
    printf("Result\n");
    printMatrix(result, size, size);

    freeMatrix(first);
    freeMatrix(second);
    freeMatrix(result);
}

void testMultiplication(
        void (*multiplySerial)(double **, double **, double **, const int),
        void (*multiplyParallel)(double **, double **, double **, const int)
){
    double t_start, t_finish;
    for(int k = 0; k < sizesSize; k++){
        const int size = sizes[k];
        std::cout << size;
        double **first = createMatrix(size, size);
        double **second = createMatrix(size, size);
        double **resultSerial = createMatrix(size, size);
        double **resultParallel = createMatrix(size, size);

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                first[i][j] = i + 1;
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                second[i][j] = j + 1;
            }
        }

        t_start = omp_get_wtime();
        multiplyParallel(resultSerial, first, second, size);
        t_finish = omp_get_wtime();
        std::cout << "\t" << (t_finish-t_start);

        t_start = omp_get_wtime();
        multiplySerial(resultParallel, first, second, size);
        t_finish = omp_get_wtime();
        std::cout << "\t" << t_finish-t_start;

        std::cout << "\t" << isCorrect(resultSerial, resultParallel, size, size) << std::endl;

        freeMatrix(first);
        freeMatrix(second);
        freeMatrix(resultSerial);
        freeMatrix(resultParallel);
    }

}

int main() {
    testMultiplication(&recursiveInPlace::multiplySerial, &recursiveInPlace::multiplyParallel);
    getch();
    return 0;
}