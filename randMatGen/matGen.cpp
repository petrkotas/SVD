#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {

    int m = 0, n = 0;
    string file;

    if (argc > 1 && argc < 5) {
        file = argv[1];
        m    = atoi(argv[2]);
        n    = atoi(argv[3]);
    } else {
        cerr << "bad nput args" << endl;
        return 1;
    }

   
    cout << "Generating random matrix of size : " << endl;
    cout << "m : " << m << ", n : " << n << endl;
    cout << "Output file is : " << file << endl;

    ofstream f(file.data(), ios::out|ios::binary);

    if (!f.is_open()) {
        cerr << "error openning file" << endl;
        return 2;
    }

    int size[2] = {m,n};

    f.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));

    time_t sec;
    time(&sec);
    srand((unsigned int) sec);

    double *r = new double[n];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            r[j] = ((double)rand()/(double)RAND_MAX);
        }
        f.write(reinterpret_cast<char *>(r), n*sizeof(double));
    }

    f.close();
    cout << "Done genrating matrix";

    delete [] r;
    return 0;
}
