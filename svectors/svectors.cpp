#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

#define sign(a) (a>0) ? 1 : -1
#define null 0
//#define debug

void tsUV(double *a, double *b, double *s, double *u, double *v, int n, int ns);
void normalize(double *v, int m, int n);
void tFwrd(double *x, double *l, double *b, int n);
void tBckwd(double *x, double *u, double *b, int n);
void tLU(double *lp, double *up, double *bp, double *lm, double *um, double *bm, double *a, double *b, double s, int n);
bool writeMatrix(const char *file, double *a, int m, int n);


int main(int argc, char *argv[]) {

    char *bMat, *sNums;
    int m = 0, n = 0, ns = 0;
    int size[2];
    double *s   = NULL;
    double *a   = NULL, *b   = NULL;
    double *u   = NULL, *v   = NULL;
    long start = 0, stop = 0;

    if (argc != 3) {
        cerr << "Bad number of input parameters.";
        return 1;
    } else {
        sNums = argv[1];
        bMat  = argv[2];
    }

    // nacteni dat ze souboru, konretne nacteme singularni cisla a bidiagonalu
    ifstream fs(sNums, ios::in|ios::binary);
    fs.read(reinterpret_cast<char *>(size), 2*sizeof(int));

    ns = size[0];
    s = new double[ns];
    for (int i = 0; i < ns; i++) {
        s[i] = 0.0;
    }
    fs.read(reinterpret_cast<char *>(s), ns*sizeof(double));

    fs.close();
    fs.clear();

    fs.open(bMat, ios::in|ios::binary);
    fs.read(reinterpret_cast<char *>(size), 2*sizeof(int));

    n = size[0];
    a = new double[n];
    b = new double[n-1];

    fs.read(reinterpret_cast<char *>(a), n*sizeof(double));
    fs.read(reinterpret_cast<char *>(b), (n-1)*sizeof(double));

#ifdef debug
    for (int i = 0; i < n-1; i++) {
        cout << b[i] << " ";
    }
    cout << endl;

    cin.get();
#endif

    // data jsou nactena bude se pocitat  v teto chvili uz kazdy proces ma sva data nactena, kazdy ma nacteno vsechno, bo 
    // ta "trocha dat" :)

    u = new double[ns*n];
    v = new double[ns*n];

    start = clock();
    tsUV(a, b, s, u, v, n, ns);
    stop = clock();
    cout << "Left/Right singular vectors of B computed in " << ((double)(stop - start))/CLOCKS_PER_SEC << " s" << endl;

    // a ulozeni dat zpatky do souboru :)
    string o = "u_s.bin";
    writeMatrix(o.data(), u, ns, n);
    o = "v_s.bin";
    writeMatrix(o.data(), v, ns, n);

    delete [] u; delete [] v;
    delete [] s;
    delete [] a; delete [] b;

    return 0;
}

void tsUV(double *a, double *b, double *s, double *u, double *v, int n, int ns) {
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

    /*tLU(lp, up, bp, lm, um, bm, a, b, s[0]*s[0], n);
    for (int i = 0; i < n; i++) {
        cout << lp[i] << ", ";
    } 
    cout << endl;
    for (int i = 0; i < (2*n-1); i++) { 
        cout << up[i] << ", ";
    } 
    cout << endl; cin.get();*/

    for (int k = 0; k < ns; k++) {
        for (int i = 0; i < n; i++) {
            bp[i] = 0.0;
            bm[i] = 0.0;
        }
        tLU(lp, up, bp, lm, um, bm, a, b, s[k]*s[k], n);
        tFwrd(x, lp, bp, n);
        /*for (int i = 0; i < n; i++) {
            cout << bp[i] << " ";
        } cout << endl;*/
        tBckwd(&v[k*n], up, x, n);
        tFwrd(x, lm, bm, n);
        tBckwd(&u[k*n], um, x, n);
    }

    normalize(u, ns, n);
    normalize(v, ns, n);

    delete [] lp;
    delete [] lm;
    delete [] up;
    delete [] um;
    delete [] bp;
    delete [] bm;
    delete [] x;
}

void normalize(double *v, int m, int n) {
    double *sum = new double[n];
    for (int i = 0; i < n; i++) {
        sum[i] = 0;
    }
    for (int i = 0; i < m*n; i++) {
        sum[i/n] += v[i]*v[i];
    }
    for (int i = 0; i < n; i++) {
        sum[i] = sqrt(sum[i]);
    }
    for (int i = 0; i < m*n; i++) {
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
    double eps = 1e-7;   
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
        bsumm++;
    }

    if (bsump == 0) {
        bp[0] = 1.0;
    }
    if (bsumm == 0) {
        bm[0] = 1.0;
    }
}

bool writeMatrix(const char *file, double *a, int m, int n) {

    ofstream f(file, ios::out|ios::binary);
    int size[2] = {m, n};

    if (!f.is_open()) {
        return false;
    }

    f.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));
    f.write(reinterpret_cast<char *>(a), m*n*sizeof(double));

    f.close();

    return true;
}
