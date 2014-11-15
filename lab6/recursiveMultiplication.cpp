#include <iostream>
#include <conio.h>
#include "utils.h"

/* Cormen, Leiserson, Rivest, Stein. Introduction to Algorithms, 2nd Ed.
* Chapter 28.  Matrix Operations
*       28.2 Strassen's algorithm for createMatrix multiplication
* Recursive Algorithms for Matrix Multiplication
*       result = first * second
*/
namespace recursive {
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
    void multParallel(double **result, double **first, double **second, int size) {
        if (size == 2) {
            double a, b, e, f, c, d, g, h;

            a = first[0][0];
            b = first[0][1];
            c = first[1][0];
            d = first[1][1];

            e = second[0][0];
            f = second[0][1];
            g = second[1][0];
            h = second[1][1];

            result[0][0] = a*e + b*g;   // r = ae + bg
            result[0][1] = a*f + b*h;   // s = af + bh
            result[1][0] = c*e + d*g;   // t = ce + dg
            result[1][1] = c*f + d*h;   // u = cf + dh
        } else {
            int i, j;
            double **a, **b, **e, **f, **ae, **bg, **af, **bh;
            double **c, **d, **g, **h, **ce, **dg, **cf, **dh;

            size >>= 1;
            // create subMatrix
            a = createMatrix(size, size);
            b = createMatrix(size, size);
            c = createMatrix(size, size);
            d = createMatrix(size, size);

            e = createMatrix(size, size);
            f = createMatrix(size, size);
            g = createMatrix(size, size);
            h = createMatrix(size, size);

            ae = createMatrix(size, size);
            bg = createMatrix(size, size);
            af = createMatrix(size, size);
            bh = createMatrix(size, size);

            ce = createMatrix(size, size);
            dg = createMatrix(size, size);
            cf = createMatrix(size, size);
            dh = createMatrix(size, size);

            // initialize subMatrix
            #pragma omp parallel for
            for (i=0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    a[i][j] = first[i][j];                  // first
                    b[i][j] = first[i][j + size];           //   | a  b |
                    c[i][j] = first[i + size][j];           //   |      |
                    d[i][j] = first[i + size][j + size];    //   | c  d |

                    e[i][j] = second[i][j];                 // second
                    f[i][j] = second[i][j + size];          //   | e  f |
                    g[i][j] = second[i + size][j];          //   |      |
                    h[i][j] = second[i + size][j + size];   //   | g  h |
                }
            }
            #pragma omp wait

            #pragma omp task firstprivate(ae, a, e, size)
            multParallel(ae, a, e, size);   // ae = a x e
            #pragma omp task firstprivate(bg, b, g, size)
            multParallel(bg, b, g, size);   // bg = b x g

            #pragma omp task firstprivate(af, a, f, size)
            multParallel(af, a, f, size);   // af = a x f
            #pragma omp task firstprivate(bh, b, h, size)
            multParallel(bh, b, h, size);   // bh = b x h

            #pragma omp task firstprivate(ce, c, e, size)
            multParallel(ce, c, e, size);   // ce = c x e
            #pragma omp task firstprivate(dg, d, g, size)
            multParallel(dg, d, g, size);   // dg = d x g

            #pragma omp task firstprivate(cf, c, f, size)
            multParallel(cf, c, f, size);   // cf = c x f
            #pragma omp task firstprivate(dh, d, h, size)
            multParallel(dh, d, h, size);   // dh = d x h

            #pragma omp taskwait

            #pragma omp parallel for
            for (i=0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    result[i][j]                = ae[i][j] + bg[i][j];  // r = ae + bg
                    result[i][j + size]         = af[i][j] + bh[i][j];  // s = af + bh
                    result[i + size][j]         = ce[i][j] + dg[i][j];  // t = ce + dg
                    result[i + size][j + size]  = cf[i][j] + dh[i][j];  // u = cf + dh
                }
            }

            freeMatrix(a);
            freeMatrix(b);
            freeMatrix(c);
            freeMatrix(d);

            freeMatrix(e);
            freeMatrix(f);
            freeMatrix(g);
            freeMatrix(h);

            freeMatrix(ae);
            freeMatrix(bg);
            freeMatrix(af);
            freeMatrix(bh);

            freeMatrix(ce);
            freeMatrix(dg);
            freeMatrix(cf);
            freeMatrix(dh);
        }
    }

    void multiplyParallel(double **result, double **first, double **second, int size) {
        #pragma omp parallel
        {
            #pragma omp single nowait
            multParallel(result, first, second, size);
        }
    }

    void multiplySerial(double **result, double **first, double **second, int size) {
        if (size == 2) {
            double a, b, e, f, c, d, g, h;

            a = first[0][0];
            b = first[0][1];
            c = first[1][0];
            d = first[1][1];

            e = second[0][0];
            f = second[0][1];
            g = second[1][0];
            h = second[1][1];

            result[0][0] = a*e + b*g;   // r = ae + bg
            result[0][1] = a*f + b*h;   // s = af + bh
            result[1][0] = c*e + d*g;   // t = ce + dg
            result[1][1] = c*f + d*h;   // u = cf + dh
        } else {
            int i, j;
            double **a, **b, **e, **f, **ae, **bg, **af, **bh;
            double **c, **d, **g, **h, **ce, **dg, **cf, **dh;

            size >>= 1;
            // create subMatrix
            a = createMatrix(size, size);
            b = createMatrix(size, size);
            c = createMatrix(size, size);
            d = createMatrix(size, size);

            e = createMatrix(size, size);
            f = createMatrix(size, size);
            g = createMatrix(size, size);
            h = createMatrix(size, size);

            ae = createMatrix(size, size);
            bg = createMatrix(size, size);
            af = createMatrix(size, size);
            bh = createMatrix(size, size);

            ce = createMatrix(size, size);
            dg = createMatrix(size, size);
            cf = createMatrix(size, size);
            dh = createMatrix(size, size);

            // initialize subMatrix
            for (i=0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    a[i][j] = first[i][j];                  // first
                    b[i][j] = first[i][j + size];           //   | a  b |
                    c[i][j] = first[i + size][j];           //   |      |
                    d[i][j] = first[i + size][j + size];    //   | c  d |

                    e[i][j] = second[i][j];                 // second
                    f[i][j] = second[i][j + size];          //   | e  f |
                    g[i][j] = second[i + size][j];          //   |      |
                    h[i][j] = second[i + size][j + size];   //   | g  h |
                }
            }

            multiplySerial(ae, a, e, size);   // ae = a x e
            multiplySerial(bg, b, g, size);   // bg = b x g

            multiplySerial(af, a, f, size);   // af = a x f
            multiplySerial(bh, b, h, size);   // bh = b x h

            multiplySerial(ce, c, e, size);   // ce = c x e
            multiplySerial(dg, d, g, size);   // dg = d x g

            multiplySerial(cf, c, f, size);   // cf = c x f
            multiplySerial(dh, d, h, size);   // dh = d x h

            for (i=0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    result[i][j]                = ae[i][j] + bg[i][j];  // r = ae + bg
                    result[i][j + size]         = af[i][j] + bh[i][j];  // s = af + bh
                    result[i + size][j]         = ce[i][j] + dg[i][j];  // t = ce + dg
                    result[i + size][j + size]  = cf[i][j] + dh[i][j];  // u = cf + dh
                }
            }

            freeMatrix(a);
            freeMatrix(b);
            freeMatrix(c);
            freeMatrix(d);

            freeMatrix(e);
            freeMatrix(f);
            freeMatrix(g);
            freeMatrix(h);

            freeMatrix(ae);
            freeMatrix(bg);
            freeMatrix(af);
            freeMatrix(bh);

            freeMatrix(ce);
            freeMatrix(dg);
            freeMatrix(cf);
            freeMatrix(dh);
        }
    }
}