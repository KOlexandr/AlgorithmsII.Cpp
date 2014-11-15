#include <iostream>
#include <conio.h>

/* Cormen, Leiserson, Rivest, Stein. Introduction to Algorithms, 2nd Ed.
* Chapter 28.  Matrix Operations
*       28.2 Strassen's algorithm for createMatrix multiplication
* Recursive (in place) Algorithms for Matrix Multiplication
*       result = first * second
*/
namespace recursiveInPlace {
    /****************************************************************************/
    /*      result      =     first     *    second                             */
    /*                                                                          */
    /*   | r  s |     | a  b |   | e  f |                                       */
    /*   |      |  =  |      | * |      |                                       */
    /*   | t  u |     | c  d |   | g  h |                                       */
    /*                                                                          */
    /*   r = ae + bg                             -------> j                     */
    /*   s = af + bh                             |                              */
    /*   t = ce + dg                             |                              */
    /*   u = cf + dh                             | i                            */
    /****************************************************************************/
    void multParallel(double **result, double **first, double **second,
            int iR, int jR, int iF, int jF, int iS, int jS, int size) {
        if (size == 2) {
            double a, b, e, f, c, d, g, h;

            a = first[iF][jF];
            b = first[iF][jF+1];
            c = first[iF+1][jF];
            d = first[iF+1][jF+1];

            e = second[iS][jS];
            f = second[iS][jS+1];
            g = second[iS+1][jS];
            h = second[iS+1][jS+1];


            result[iR][jR] += a*e + b*g;        // r = ae + bg
            result[iR][jR+1] += a*f + b*h;      // s = af + bh
            result[iR+1][jR] += c*e + d*g;      // t = ce + dg
            result[iR+1][jR+1] += c*f + d*h;    // u = cf + dh
        } else {
            size >>= 1;    // size = size div 2
            // r=a=e=[0][0] s=b=f=[0][size]
            // t=c=g=[size][0] u=d=h=[size][size]
            #pragma omp task
            multParallel(result, first, second, iR, jR, iF, jF, iS, jS, size);            // r = ae +
            #pragma omp task
            multParallel(result, first, second, iR, jR, iF, jF+size, iS+size, jS, size);  //        + bg

            #pragma omp task
            multParallel(result, first, second, iR, jR+size, iF, jF, iS, jS+size, size);            // s = af +
            #pragma omp task
            multParallel(result, first, second, iR, jR+size, iF, jF+size, iS+size, jS+size, size);  //        + bh

            #pragma omp task
            multParallel(result, first, second, iR+size, jR, iF+size, jF, iS, jS, size);            // t = ce +
            #pragma omp task
            multParallel(result, first, second, iR+size, jR, iF+size, jF+size, iS+size, jS, size);  //        + dg

            #pragma omp task
            multParallel(result, first, second, iR+size, jR+size, iF+size, jF, iS, jS+size, size);            // u = cf +
            #pragma omp task
            multParallel(result, first, second, iR+size, jR+size, iF+size, jF+size, iS+size, jS+size, size);  //        + dh

            #pragma omp taskwait
        }
    }

    void multiplyParallel(double **result, double **first, double **second, int size) {
        #pragma omp parallel
        {
            #pragma omp single nowait
            multParallel(result, first, second, 0, 0, 0, 0, 0, 0, size);
        }
    }

    void multSerial(double **result, double **first, double **second,
            int iR, int jR, int iF, int jF, int iS, int jS, int size) {
        if (size == 2) {
            double a, b, e, f, c, d, g, h;

            a = first[iF][jF];
            b = first[iF][jF+1];
            c = first[iF+1][jF];
            d = first[iF+1][jF+1];

            e = second[iS][jS];
            f = second[iS][jS+1];
            g = second[iS+1][jS];
            h = second[iS+1][jS+1];


            result[iR][jR] += a*e + b*g;        // r = ae + bg
            result[iR][jR+1] += a*f + b*h;      // s = af + bh
            result[iR+1][jR] += c*e + d*g;      // t = ce + dg
            result[iR+1][jR+1] += c*f + d*h;    // u = cf + dh
        } else {
            size >>= 1;    // size = size div 2
            // r=a=e=[0][0] s=b=f=[0][size]
            // t=c=g=[size][0] u=d=h=[size][size]
            multSerial(result, first, second, iR, jR, iF, jF, iS, jS, size);            // r = ae +
            multSerial(result, first, second, iR, jR, iF, jF+size, iS+size, jS, size);  //        + bg

            multSerial(result, first, second, iR, jR+size, iF, jF, iS, jS+size, size);            // s = af +
            multSerial(result, first, second, iR, jR+size, iF, jF+size, iS+size, jS+size, size);  //        + bh

            multSerial(result, first, second, iR+size, jR, iF+size, jF, iS, jS, size);            // t = ce +
            multSerial(result, first, second, iR+size, jR, iF+size, jF+size, iS+size, jS, size);  //        + dg

            multSerial(result, first, second, iR+size, jR+size, iF+size, jF, iS, jS+size, size);            // u = cf +
            multSerial(result, first, second, iR+size, jR+size, iF+size, jF+size, iS+size, jS+size, size);  //        + dh
        }
    }

    void multiplySerial(double **result, double **first, double **second, int size) {
        multSerial(result, first, second, 0, 0, 0, 0, 0, 0, size);
    }
}