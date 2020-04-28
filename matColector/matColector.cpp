#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

string itoa (int n);

int main(int argc, char *argv[]) {
    
    int n = 0, pr = 0, buffLen = 0;
    double *a, *b, *buff;
    int size[2];
    string f, tf;
    fstream fs;


    if (argc != 3) {
        f = argv[1];
        n = atoi(argv[2]);
        pr = atoi(argv[3]);
    } else {
        cerr << "bad nput args" << endl;
        return 1;
    }
    
    // *******************************************************************************
    //  Sestaveni matice B
    // *******************************************************************************
    buffLen = n/pr;

    a = new double[n];
    b = new double[n];
    buff = new double[buffLen];

    cout << "Rading series of B_ matrices" << endl;
    cout << "Input size of matrices is " << buffLen << "x2" << endl;

    for (int k = 0; k < pr; k++) {
        tf = f + itoa(k) + ".bin";
        fs.clear();
        fs.open(tf.data(), ios::in|ios::binary);
        if (!fs.is_open()) {
            return 1;
        }
        fs.read(reinterpret_cast<char *>(size), 2*sizeof(int));
        fs.read(reinterpret_cast<char *>(buff), buffLen*sizeof(double));

        for (int i = 0; i < buffLen; i++) {
            a[(i*pr) + k] = buff[i];
        }
        fs.read(reinterpret_cast<char *>(buff), buffLen*sizeof(double));
        for (int i = 0; i < buffLen; i++) {
            b[(i*pr) + k] = buff[i];
        }

        fs.close();
        fs.flush();
        cout << "Matrix B_" << k << " done" << endl;
    }

    fs.clear();
    tf = f + "A.bin";
    fs.open(tf.data(), ios::out|ios::binary);

    size[0] = n; size[1] = 2;

    fs.write(reinterpret_cast<char *>(size), 2*sizeof(int));
    fs.write(reinterpret_cast<char *>(a), n*sizeof(double));
    fs.write(reinterpret_cast<char *>(b), n*sizeof(double));

    fs.close();

    cout << "All done, size of output matrix is " << n << "x2" << endl;
    cin.get();
    /*for (int i = 0; i < n; i++) {
        cout << a[i] << " , " << b[i] << endl;
    }
    cin.get();*/
    delete [] a;
    delete [] b;

    // *******************************************************************************
    //  Sestaveni matice U
    // *******************************************************************************


    
    return 0;
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
        }
 
        while ((n /= 10) > 0);
 
        for (int j = i-1; j >= 0; j--) { 
                u += s[j]; //building our string number
        }
 
        delete[] s; //free-up the memory!
        return u;
}
