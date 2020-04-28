#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <ctime>
#include "mkl.h"
using namespace std;

#pragma warning(disable : 4996)
#define eps 1e-7
//#define debug


void svdB(double *a, double *b, int n, int &it);
void givens(double &cs, double &sn, double &r, double f, double g);
void dqd(double *a, double *b, int len);

int main(int argc, char *argv[]) {
    
    char *bMat;
    int n = 0, it = 0, scount = 0;
    double *s   = NULL;
    double *a   = NULL, *b   = NULL;
    long start = 0, stop = 0;
    int size[2];

	

    if (argc != 2) {
        cerr << "Bad number of input parameters.";
        return 1;
    } else {
        bMat  = argv[1];
    }

    // nacteni bidiagonaly ze souboru
    ifstream fs(bMat, ios::in|ios::binary);
    fs.read(reinterpret_cast<char *>(size), 2*sizeof(int));

    n = size[0];
    a = new double[n];
    b = new double[n-1];

    fs.read(reinterpret_cast<char *>(a), n*sizeof(double));
    fs.read(reinterpret_cast<char *>(b), (n-1)*sizeof(double));

    fs.close();
#ifdef debug
    for (int i = 0; i < n; i++)
        cout << a[i] << " ";
    cout << endl;
#endif
    // a zahajeni rotaci
    start = clock();
    svdB(a, b, n, it);
    stop = clock();
    cout << "Singular values computed in " << ((double)(stop - start))/CLOCKS_PER_SEC << " s" << endl;
    cout << "Iteration count is " << it << endl;

#ifdef debug
    for (int i = 0; i < n; i++)
        cout << a[i] << " ";
    cout << endl;
#endif

    for (int i = 0; i < n; i++) {
        a[i] = abs(a[i]);
        if (a[i] > eps) 
            scount++;
    }

    // zapis spocitanych cisel
    char * file = "s.bin";
    ofstream ofs(file, ios::out|ios::binary);
    size[0] = scount; size[1] = 1;
    cout << "otevreny soubor" << endl;
    ofs.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));
    ofs.write(reinterpret_cast<char *>(a), scount*sizeof(double));
    cout << "zapasno" << endl;
    ofs.close();
    cout << "soubor zavren" << endl;
    delete [] a; delete [] b; delete [] s;
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
        for (int i = 1; i < n; i++) {
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
