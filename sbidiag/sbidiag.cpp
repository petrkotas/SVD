#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <ctime>
#include "mkl.h"

using namespace std;

#define sign(a)  (a>0) ? 1 : -1
#define max(a,b) (a>b) ? a : b
#define min(a,b) (a>b) ? b : a
#define null 0

void bidiag(double **a, double **u, double **v, int row, int col);
double** acumulateU(double **u, int m);
double** acumulateV(double **v, int n);
double** readMatrix(char *file, int &m, int &n);
bool writeMatrix(char *file, double **a, int m, int n);
bool writeMatrixBidiagonal(char *file, double **a, int m, int n);

int main(int argc, char *argv[]){

    char *file;
    double **a = null, **u = null, **v = null;
    double **U = null, **V = null;
    int m, n;
    long start = 0, stop = 0;


    if (argc != 2) {
        cerr << "Bad number of input parameters";
        return 1;
    } else {
        file = argv[1];
    }

    a = readMatrix(file, m, n); // nacteni matice a alokace pameti!!!
    cout << "A loaded." << endl;
    cout << "m = " << m << ", n = " << n << endl;
   /* for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
    cin.get();*/
    u = new double*[m];
    for (int i = 0; i < m; i++)
        u[i] = new double[m];
    v = new double*[n];
    for (int i = 0; i < n; i++) {
        v[i] = new double[n];
        for (int j = 0; j < n; j++) {
            v[i][j] = 0.0;
        }
    }

    start = clock();
    bidiag(a, u, v, m, n);
    stop = clock();
    cout << "Bidiag done in " << ((double)(stop-start))/CLOCKS_PER_SEC << " s" << endl;
    /*start = clock();
    U = acumulateU(u, m);
    stop = clock();
    cout << "U accumulated in " << ((double)(stop-start))/CLOCKS_PER_SEC << " s" << endl;
    start = clock();
    V = acumulateV(v, n);
    stop = clock();
    cout << "V accumulated in " << ((double)(stop-start))/CLOCKS_PER_SEC << " s" << endl;*/
    //writeMatrixBidiagonal("B.bin", a, m, n);
    writeMatrix("B.bin", a, m, n);
    //writeMatrix("U.bin", U, m, m);
    //writeMatrix("V.bin", V, n, n);

    for (int i = 0; i < m; i++)
        delete [] a[i];
    delete [] a;
    for (int i = 0; i < m; i++)
        delete [] u[i];
    delete [] u;
    for (int i = 0; i < n; i++)
        delete [] v[i];
    delete [] v;
   /* for (int i = 0; i < m; i++)
        delete [] U[i];
    delete [] U;
    for (int i = 0; i < n; i++)
        delete [] V[i];
    delete [] V;*/
    cin.get();
    return 0;
}

void bidiag(double **a, double **u, double **v, int row, int col) {
    double alpha = 0.0;
    double beta  = 0.0;
    double sbeta = 0.0;
    double gamma = 0.0;
    double norm  = 0.0;
    int    maxMN = max(row,col);
    double *vk   = new double[maxMN];
    double eps   = 1e-7;
    int    minMN = min(row, col);

    for (int k = 0; k < minMN; k++) {
        alpha = (double)(sign(a[k][k]));
        norm  = a[k][k]*a[k][k];
        for (int i = (k+1); i < row; i++) {
            norm += a[i][k]*a[i][k];
        }
        alpha *= sqrt(norm);
        for (int i = 0; i < row; i++) {
            vk[i] = (i < k) ? 0 : a[i][k];
        }
        beta   = 2*(norm + vk[k]*alpha);
        sbeta  = (beta > eps) ? sqrt(beta) : 1.0;
        vk[k] += alpha;
        for (int j = 0; j < row; j++) {
            u[k][j] = vk[j]/sbeta;
        }
        if (beta < eps)
            continue;
        for (int j = k; j < col; j++) {
            gamma = 0.0;
            for (int i = 0; i < row; i++) {
                gamma += a[i][j]*vk[i];
            }
            for (int i = 0; i < row; i++) {
                a[i][j] = a[i][j] - (2*gamma/beta)*vk[i];
            }
        }

        if (k >= (col-2))
            continue;

        alpha = (double)(sign(a[k][k+1]));
        norm = a[k][k+1]*a[k][k+1];
        for (int j = (k+2); j < col; j++) {
            norm += a[k][j]*a[k][j];
        }
        alpha *= sqrt(norm);
        for (int j = 0; j < col; j++) {
            vk[j] = (j < (k+1)) ? 0 : a[k][j];
        }
        beta     = 2*(norm + vk[k+1]*alpha);
        sbeta  = (beta > eps) ? sqrt(beta) : 1.0;
        vk[k+1] += alpha;
        for (int j = 0; j < col; j++) {
            v[k][j] = vk[j]/sbeta;
        }
        if (beta < eps)
            continue;
        for (int i = k; i < row; i++) {
            gamma = 0.0;
            for (int j = 0; j < col; j++) {
                gamma += a[i][j]*vk[j];
            }
            for (int j = 0; j < col; j++) {
                a[i][j] -= (2*gamma/beta)*vk[j];
            }
        }
    }
    delete [] vk;
}

double** acumulateU(double **u, int m) {
    double **U = new double*[m];
    for (int i = 0; i < m; i++) {
        U[i] = new double[m];
        for (int j = 0; j < m; j++) {
            U[i][j] = (i==j) ? 1.0 : 0.0;
        }
    }

    double t = 0.0;
    double *tsum = new double[m];

    for (int k = 0; k < m; k++) {

        for (int j = 0; j < m; j++) {

            for (int i = 0; i < m; i++)
                tsum[i] = 0.0;

            for (int i = 0; i < m; i++) {
                
                for (int l = 0; l < m; l++) {
                    t = ((i==l) ? 1.0 : 0.0) - 2*u[k][i]*u[k][l];
                    tsum[i] += t*U[l][j];
                }

            }

            for (int i = 0; i < m; i++)
                U[i][j] = tsum[i];

        }

    }

    return U;
}

double** acumulateV(double **v, int n) {
    double **V = new double*[n];
    for (int i = 0; i < n; i++) {
        V[i] = new double[n];
        for (int j = 0; j < n; j++) {
            V[i][j] = (i==j) ? 1.0 : 0.0;
        }
    }

    double t = 0.0;
    double *tsum = new double[n];

    for (int k = 0; k < n; k++) {

        for (int i = 0; i < n; i++) {
        
            for (int j = 0; j < n; j++)
                tsum[j] = 0.0;

            for (int j = 0; j < n; j++) {
            
                for (int l = 0; l < n; l++) {
                    t = ((j==l) ? 1.0 : 0.0) - 2*v[k][l]*v[k][j];
                    tsum[j] += V[i][l]*t;
                }

            }

            for (int j = 0; j < n; j++)
                V[i][j] = tsum[j];

        }

    }

    return V;
}

double** readMatrix(char *file, int &m, int &n) {

    double **a;
    fstream f(file, ios::in|ios::binary);
    int size[2];

    if (!f.is_open()) {
        return null;
    }

    f.read(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));

    m = size[0]; n = size[1];
    a = new double*[m];
    
    for (int i = 0; i < m; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++)
            a[i][j] = 0.0;
        f.read(reinterpret_cast<char *>(a[i]), n*sizeof(double));
    }

    f.close();
    return a;
}

bool writeMatrix(char *file, double **a, int m, int n) {

    ofstream f(file, ios::out|ios::binary);
    int size[2] = {m, n};

    if (!f.is_open()) {
        return false;
    }

    f.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));
    
    for (int i = 0; i < m; i++) {
        f.write(reinterpret_cast<char *>(a[i]), n*sizeof(double));
    }

    f.close();

    return true;
}

bool writeMatrixBidiagonal(char *file, double **a, int m, int n) {

    int minMN = min(m,n);

    ofstream f(file, ios::out|ios::binary);
    int size[2] = {minMN, 2};

    if (!f.is_open()) {
        return false;
    }

    f.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));
    
    double *tmp = new double[2*minMN];

    for (int i = 0; i < minMN; i++)
        tmp[i] = a[i][i];

    for (int i = 0; i < minMN-1; i++)
        tmp[m+i] = a[i][i+1];

    tmp[2*minMN-1] = 0.0;

    f.write(reinterpret_cast<char *>(tmp), (2*minMN)*sizeof(double));

    delete [] tmp;
    f.close();

    return true;
}

