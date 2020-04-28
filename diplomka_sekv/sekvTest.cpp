#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

using namespace std;

#define sign(a) (a>0) ? 1 : -1
#define null 0


void dqd(double *a, double *b, int len);
string itoa (int n);
void givens(double &cs, double &sn, double &r, double f, double g);
void svdB(double *a, double *b, int len, int &it);
void tLU(double *lp, double *up, double *bp, double *lm, double *um, double *bm, double *a, double *b, double s, int n);
void tsUV(double *a, double *b, double *s, double *u, double *v, int n);
void tFwrd(double *x, double *l, double *b, int n);
void tBckwd(double *x, double *u, double *b, int n);
void normalize(double *v, int n);
void bidiag(double **a, double **u, double **v, int row, int col);
double** acumulateU(double **u, int m);
double** acumulateV(double **v, int n);
double** readMatrix(char *file, int &m, int &n);
double** readMatrix(char *file, int &m, int &n, int offset, int k);
bool writeMatrix(char *file, double **a, int m, int n);

int main(int argc, char *argv[]){
    
    int m = 0, n = 0;
    int mn[2];
    double *s;

    ifstream fs("s.bin", ios::in|ios::binary);
    fs.read(reinterpret_cast<char *>(&mn), 2*sizeof(int));

    s = new double[mn[0]];
    for (int i = 0; i < mn[0]; i++) {
        s[i] = 0.0;
    }
    fs.read(reinterpret_cast<char *>(&s[0]), mn[0]*sizeof(double));

    for (int i = 0; i < mn[0]; i++) {
        cout << s[i] << " ";
    }

    cout << endl;
    cin.get();

    //bidiag(a, u, v, m, n);

    //U = acumulateU(u, m);
    //V = acumulateV(v, n);

    //cout << m << " " << n << endl;

    //for (int i = 0; i < m; i++) {
    //    for (int j = 0; j < n; j++) {
    //        cout << a[i][j] << " ";   
    //    }
    //    cout << endl;
    //}

    //for (int i = 0; i < m; i++) {
    //    for (int j = 0; j < n; j++) {
    //        cout << U[i][j] << " ";   
    //    }
    //    cout << endl;
    //}

    //writeMatrix("B.bin", a, m, n);
    ////writeMatrix("u.bin", u, m, n);
    //writeMatrix("U.bin", U, m, m);
    //writeMatrix("V.bin", V, n, n);


    return 0;
}

void svdB(double *a, double *b, int n, int &it) {

    double cs, sn, oldcs, oldsn;
    double h = 0.0, r = 0.0;
    double ref = 0.0;
    double tol = 1e-7;
    double norm = 9999;

    for (int i = 0; i < n; i++) {
        if (ref < abs(a[i])) {
            ref = abs(a[i]);
        }
    }

    tol *= ref;

    while (norm >= tol) {
        cs = 1.0; sn = 0.0;
        oldcs = 1.0; oldsn = 1.0;

        for (int i = 0; i < (n-1); i++) {      
            givens(cs, sn, r, (a[i]*cs), b[i]);          
            if (i != 0) 
                b[i-1] = oldsn*r;        
            givens(oldcs, oldsn, a[i], (oldcs*r), (a[i+1]*sn));
        }

        h        = a[n-1]*cs;
        b[n-2]   = h*oldsn;
        a[n-1]   = h*oldcs;

        norm = 0.0;
        for (int i = 0; i < n; i++) {
            if ((a[i] + a[i-1] + b[i-1]) == (a[i] + a[i-1])) {
                b[i-1] = 0.0;
            }
        }
        for (int i = 0; i < (n-1); i++) {
            norm += b[i]*b[i];
        }
        norm = sqrt(norm);
        it++;
    }
}

void givens(double &cs, double &sn, double &r, double f, double g) {
    double t = 0.0, tt = 0.0;

    if (f==0) {
        cs = 0.0;
        sn = 1.0;
        r = g;
    } else if (abs(f) > abs(g)) {
        t  = g / f;
        tt = sqrt(1+t*t);
        cs = 1 / tt;
        sn = t * cs;
        r  = f * tt;
    } else {
        t  = f / g;
        tt = sqrt(1+t*t);
        sn = 1 / tt;
        cs = t * sn;
        r  = g * tt;
    }
}


void dqd(double *a, double *b, int len) {
	
	double *aa = new double[len];
	double *bb = new double[len-1];
	
	for (int i = 0; i<(len-1); i++) {
		aa[i] = a[i]*a[i];
		bb[i] = b[i]*b[i];
	}
	aa[len-1] = a[len-1]*a[len-1];

	double t = 0.0;
	double tf = 0.0;
	double norm = 0.0;
    double eps  = 1e-8;
    int it = 0;
    
    for(int i = 0; i < len-1; i++) {
    	if(norm < bb[i]) {
    		norm = bb[i];
    	}
    }
        
    while ((norm >= eps) && (it < 1000)) {
    	t = aa[0];
    	for(int i = 0; i < len-1; i++) {
    		aa[i] = t + bb[i];
    		tf = aa[i+1] / aa[i];
    		bb[i] = bb[i] * tf;
    		t = t * tf;
    	}
    	aa[len-1] = t;
    	it++;
        norm = abs(bb[0]);
        cout << bb[0] << ", ";
    	for(int i = 1; i < len-1; i++) {
	    	if(norm < bb[i]) {
	    		norm = abs(bb[i]);
	    	}
            cout << bb[i] << ", ";
	    }
    	cout << endl;
	}
    
    for(int i = 0; i < len; i++) {
    	cout << sqrt(aa[i]) << ", ";
    }

    delete [] aa;
    delete [] bb;
}

string itoa (int n) {
 
        char * s = new char[17];
        string u;
 
        if (n < 0) { //turns n positive
                n = (-1 * n); 
                u = "-"; //adds '-' on result string
        }
 
        int i=0; //s counter
  
        do {
                s[i++]= n%10 + '0'; //conversion of each digit of n to char
                n -= n%10; //update n value
        } while ((n /= 10) > 0);
 
        for (int j = i-1; j >= 0; j--) { 
                u += s[j]; //building our string number
        }
 
        delete[] s; //free-up the memory!
        return u;
}

void tsUV(double *a, double *b, double *s, double *u, double *v, int n) {
    double *lp = new double[n];
    double *lm = new double[n];
    double *up = new double[2*n-1];
    double *um = new double[2*n-1];
    double *bp = new double[n];
    double *bm = new double[n];
    double *x  = new double[n];

    for (int i = 0; i < n; i++) {
        bp[i] = bm[i] = 0.0;
    }

    for (int k = 0; k < n; k++) {
        tLU(lp, up, bp, lm, um, bm, a, b, s[k], n);
        //tFwrd(x, lp, bp, n);
        for (int i = 0; i < n-1; i++) {
            cout << lp[i] << " | " << up[2*i] << " | " << up[2*i+1] << endl;
        }
        cout << 'x' << " | " << up[2*(n-1)] << " | " << 'x' << endl;
        tBckwd(&v[k*n], up, x, n);
        tFwrd(x, lm, bm, n);
        tBckwd(&u[k*n], um, x, n);
        cin.get();
    }

    normalize(u, n);
    normalize(v, n);

    delete [] lp;
    delete [] lm;
    delete [] up;
    delete [] um;
    delete [] bp;
    delete [] bm;
    delete [] x;
}

void normalize(double *v, int n) {
    double *sum = new double[n];
    for (int i = 0; i < n; i++) {
        sum[i] = 0;
    }
    for (int i = 0; i < n*n; i++) {
        sum[i/n] += v[i]*v[i];
    }
    for (int i = 0; i < n; i++) {
        sum[i] = sqrt(sum[i]);
    }
    for (int i = 0; i < n*n; i++) {
        v[i] /= sum[i/n];
    }
    delete [] sum;
}

void tFwrd(double *x, double *l, double *b, int n) {
    x[0] = b[0];
    for (int i = 1; i < n; i++) {
        x[i] = b[i] - l[i-1]*x[i-1];
    }
}

void tBckwd(double *x, double *u, double *b, int n) {
    x[n-1] = b[n-1]/u[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (b[i] - u[2*i+1]*x[i+1])/u[2*i];
    }
}

void tLU(double *lp, double *up, double *bp, double *lm, double *um, double *bm, double *a, double *b, double s, int n) {
    int bsump = 0, bsumm = 0;
    double eps = 1e-3;   
    up[0] = pow(a[0],2) - s; up[1] = a[0]*b[0];
    um[0] = pow(a[0],2) + pow(b[0],2) - s; um[1] = a[1]*b[0];
    for (int i = 0; i < (n-2); i++) {
        lp[i]         = (a[i]*b[i])/up[2*i];
        lm[i]         = (a[i+1]*b[i])/um[2*i];
        up[2*(i+1)]   = (pow(a[i+1],2) + pow(b[i],2) - s) - lp[i]*up[2*i+1];
        um[2*(i+1)]   = (pow(a[i+1],2) + pow(b[i+1],2) - s) - lm[i]*um[2*i+1];
        up[2*(i+1)+1] = a[i+1]*b[i+1];
        um[2*(i+1)+1] = a[i+2]*b[i+1];
        if (abs(up[2*i]) <= eps) {
            up[2*i] = 1.0;
            bp[i] = 1.0;
            bsump++;
        }
        if (abs(um[2*i]) <= eps) {
            um[2*i] = 1.0;
            bm[i] = 1.0;
            bsumm++;
        }
    }
    
    lp[n-2] = a[n-2]*b[n-2]/up[2*(n-2)];
    lm[n-2] = a[n-1]*b[n-2]/um[2*(n-2)];

    up[2*(n-1)] = (pow(a[n-1],2) + pow(b[n-2],2) - s) - lp[n-2]*up[2*(n-1)-1];
    um[2*(n-1)] = (pow(a[n-1],2) - s) - lm[n-2]*um[2*(n-1)-1];

    if (abs(up[2*(n-1)]) <= eps) {
        up[2*(n-1)] = 1.0;
        bp[n-1] = 1.0;
        bsump++;
    }
    if (abs(um[2*(n-1)]) <= eps) {
        um[2*(n-1)] = 1.0;
        bm[n-1] = 1.0;
        bsump++;
    }

    if (bsump == 0) {
        bp[0] = 1.0;
    }
    if (bsumm == 0) {
        bm[0] = 1.0;
    }
}

void bidiag(double **a, double **u, double **v, int row, int col) {
    double alpha = 0.0;
    double beta  = 0.0;
    double sbeta = 0.0;
    double gamma = 0.0;
    double norm  = 0.0;
    double *vk   = new double[row];

    for (int k = 0; k < col; k++) {
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
        sbeta  = sqrt(beta);
        vk[k] += alpha;
        for (int j = 0; j < row; j++) {
            u[k][j] = vk[j]/sbeta;
        }
        for (int j = k; j < col; j++) {
            gamma = 0.0;
            for (int i = 0; i < col; i++) {
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
        sbeta    = sqrt(beta);
        vk[k+1] += alpha;
        for (int j = 0; j < col; j++) {
            v[k][j] = vk[j]/sbeta;
        }
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

double** readMatrix(char *file, int &m, int &n, int offset, int k) {

    double **a;
    fstream f(file, ios::in|ios::binary);
    int size[2];

    if (!f.is_open()) {
        return null;
    }

    f.read(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));

    m = size[0]/k; n = size[1];
    a = new double*[m];

    f.seekg(offset*n*sizeof(double), ios::cur);
    
    for (int i = 0; i < m; i++) {
        a[i] = new double[n];
        for (int j = 0; j < n; j++)
            a[i][j] = 0.0;
        f.read(reinterpret_cast<char *>(a[i]), n*sizeof(double));
        f.seekg((k-1)*n*sizeof(double), ios::cur);
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
