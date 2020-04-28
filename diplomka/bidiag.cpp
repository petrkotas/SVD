#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <mpi.h>

using namespace std;

#pragma warning(disable : 4996)

//--------------------
// aux definitions
#define tag 101
#define sign(a)  (a>0) ? 1 : -1
#define max(a,b) (a>b) ? a : b
#define min(a,b) (a>b) ? b : a
//#define debug
#define eps 1e-7
//--------------------    

int rank, pr; 
int mrow;
double t0, t1;
MPI_Status status;

string itoa (int n);
void acumulateU(double **u, double *vk, int m, int n, int rank);
void acumulateV(double **v, double *vk, int m, int n, int rank);
bool writeMatrix(const char *file, double **a, int m, int n);

int main(int argc, char *argv[]) {
	/*Deklarace potrebnych promenych*/
	int row = 0, col = 0;
    int minMN = 0, maxMN = 0;
    char* file;
    t0 = 0.0; 
    double fileTime      = 0.0;
    double bidiagTime    = 0.0, writeTime     = 0.0;
    double allgatherTime = 0.0, allreduceTime = 0.0;
    double bcastTime     = 0.0;

    if (argc == 2) {
        file = argv[1];
    } else {
        cerr << "No input file given." << endl;
        return 1;
    }

	double **myRows;                      // moje radky matice
    double **myUCols, **myVRows;          // moje cast slopcu/radku pro vektoru U/V
    double *vkU, *vkV;                    // householderovy vektory
	double alfa     = 0.0, beta   = 0.0;
    double sbeta    = 0.0;
    int    myRow    = 0,   start  = 0;
    int    myUnCols = 0, myVnRows = 0;
	int size[2];

	/*Inicializacni rutina pro zavedeni MPI*/
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &pr);


    t0 = MPI_Wtime();
    fstream fs(file, ios::in|ios::binary);
    if (!fs.is_open()) {
        cerr << "File inexistent!" << endl;
        return 1;
    }
    
    int t = 0;
    
    fs.read(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));

    row = size[0]; col = size[1];
    mrow     = row/pr; 
    myUnCols = row/pr;
    myVnRows = col/pr;

    fs.seekg(rank*col*sizeof(double), ios::cur);
    myRows = new double*[mrow];
    for (int i = 0; i < mrow; i++) {
        myRows[i] = new double[col];
        fs.read(reinterpret_cast<char *>(myRows[i]), col*sizeof(double));
        fs.seekg((pr-1)*col*sizeof(double), ios::cur);
    }
    fs.close();
    fileTime = (MPI_Wtime() - t0);

#ifdef UV
    vkU = new double[row];
    vkV = new double[col];

    t = rank*myUnCols;
    myUCols = new double*[row];
    for (int i = 0; i < row; i++) {
        myUCols[i] = new double[myUnCols];
        for (int j = 0; j < myUnCols; j++) {
            myUCols[i][j] = (i == (t+j)) ? 1.0 : 0.0;
        }
    }
   
    t = rank*myVnRows;
    myVRows = new double*[myVnRows];
    for (int i = 0; i < myVnRows; i++) {
        myVRows[i] = new double[col];
        for (int j = 0; j < col; j++) {
            myVRows[i][j] = ((t+i) == j) ? 1.0 : 0.0;
        }
    }
#endif
	
#ifdef debug
    cout << "--------------------------" << endl;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < myUnCols; j++) {
			cout << myUCols[i][j] << "|";
		}
		cout << endl;
	}
	cout << endl << "lokalni kousek matice U" << endl;
    cout << "--------------------------" << endl;
	for (int i = 0; i < myVnRows; i++) {
		for (int j = 0; j < col; j++) {
			cout << myVRows[i][j] << "|";
		}
		cout << endl;
	}
	cout << endl << "lokalni kousek matice V" << endl;

    cin.get();
#endif


#ifdef debug
    cout << "--------------------------" << endl;
	for (int i = 0; i < mrow; i++) {
		for (int j = 0; j < col; j++) {
			cout << myRows[i][j] << "|";
		}
		cout << endl;
	}
	cout << endl << "lokalni kousek matice prideleny tomuto procesu" << endl;
    cin.get();
#endif

    
    double sum2   = 0.0, gamma = 0.0;
    /*double *pCol  = (double *)calloc(mrow * sizeof(double));
    double *tv    = (double *)calloc(mrow * pr * sizeof(double));
    double *v     = (double *)calloc(mrow * pr * sizeof(double));
    double *gama  = (double *)calloc(mrow * pr * sizeof(double));
    double *sgama = (double *)calloc(mrow * pr * sizeof(double));*/

    maxMN = max(row, col);
    minMN = min(row, col);
    
    double *pCol  = new double[mrow];
    double *tv    = new double[maxMN];
    double *v     = new double[maxMN];
    double *gama  = new double[row];
    double *sgama = new double[row];

    t0 = MPI_Wtime();
	// Bidiagonalizace
	for (int k = 0; k < minMN; k++) {
   
        // *********** leva strana ************

        myRow = k/pr;  // relativni umisteni radku na jednotlivych procesorech, predpoklada se cyklicke deleni
        start = (rank < (k%pr)) ? (int) myRow + 1 : (int) myRow;

#ifdef debug
        cout << "=========================" << endl;
        cout << "rank  : " << rank  << endl;
        cout << "k     : " << k     << endl;
        cout << "k/pr  : " << k/pr  << endl;
        cout << "start : " << start << endl;
        cout << "-------------------------" << endl;
#endif

        for (int i = 0; i < mrow; i++) {
            pCol[i] = myRows[i][k];
        }
        
        t1 = MPI_Wtime();
        // sesbira veskere sloupce matice takze z nich lze napocitat potrebny householder vektor
        MPI_Allgather(pCol, mrow, MPI_DOUBLE, tv, mrow, MPI_DOUBLE, MPI_COMM_WORLD);
        allgatherTime += (MPI_Wtime() - t1);

        //free((double *)pCol); // promyslet recyklaci v cyklech, nebylo by lepsi uklidit pak?
        // nezapomenout smazat pCol a v

#ifdef debug
        cout << "--------------" << endl;
        for (int j = 0; j < (mrow*pr); j++) {
            cout << tv[j] << "|";
        }
        cout << endl << "docasny vektor v sesbirany ze vsech procesu" << endl;
#endif   

        sum2  = 0.0;
        
        // setrizeni podle cyklicke permutace plus inicializace 
        for (int i = 0; i < row; i++) {
            if (i < k) {
                v[i] = 0;
            } else {
                v[i]  = tv[i/pr + (i%pr)*mrow];
                sum2 += v[i]*v[i];
            }
        }
        
#ifdef debug
        cout << "--------------" << endl;
        for (int j = 0; j < (mrow*pr); j++) {
            cout << v[j] << "|";
        }
        cout << endl << "setrizeny vektor v jeste bez upravy" << endl;
#endif  
        
        alfa  = (double)(sign(v[k]));
        alfa *= sqrt(sum2);
        beta  = 2*(sum2 + v[k]*alfa); // v'*v
        sbeta  = (beta > eps) ? sqrt(beta) : 1.0;
        v[k] += alfa;

#ifdef UV
        for (int i = 0; i < row; i++) {
            vkU[i] = v[i]/sbeta;
        }
		
        // akumulace Householderovych matic provedena na lokalnim kousku pridelenem procesu(rank)
        // deleni je po sloupcich bo se tak nebude komunikovat takze pohoda
        acumulateU(myUCols, vkU, row, myUnCols, rank);
#endif
		
        // nemame radi deleni nulou tak na to jebem
        if (beta < eps)
            continue;

#ifdef debug
        cout << "--------------" << endl;
        cout << "beta : " << beta << endl;
        cout << "alfa : " << alfa << endl;
        cout << "--------------" << endl;
        for (int j = 0; j < (mrow*pr); j++) {
            cout << v[j] << "|";
        }
        cout << endl << "setrizeny vektor v" << endl;
#endif  

        // beta by nemela byt 0 :)
        // paralelni skalarni soucin pro vypocet rady koeficiantu
        for (int j = 0; j < row; j++) {
            gama[j] = 0.0;
            for (int i = 0; i < mrow; i++) {
                gama[j] += myRows[i][j]*v[(i*pr)+rank];
            }
        }
        
#ifdef debug
        cout << "--------" << endl;
        for (int j = 0; j < (mrow*pr); j++) {
            cout << gama[j] << "|";
        }
        cout << endl << "lokalni sumy skalarniho soucinu v*A(:,k)" << endl;
#endif   
        
        t1 = MPI_Wtime();
        // secteni vsech dilcich vysledku skalarniho soucinu do vysledneho vektoru sgama
        MPI_Allreduce(gama, sgama, (mrow*pr), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        allreduceTime += (MPI_Wtime() - t1);

#ifdef debug
        cout << "--------" << endl;
        for (int j = 0; j < (mrow*pr); j++) {
            cout << sgama[j] << "|";
        }
        cout << endl << "globalni suma skalarniho soucinu v*A(:,k)" << endl;
#endif

        // samotna eliminace patricneho sloupce
        for (int j = k; j < col; j++) {
            for (int i = 0; i < mrow; i++) {
                myRows[i][j] = myRows[i][j] - (2*sgama[j]/beta)*v[(i*pr)+rank];
            }
        }

#ifdef debug
        cout << endl << "rank : " << rank << endl;
        cout << "--------" << endl;
        for (int i = 0; i < mrow; i++) {
            for (int j = 0; j < mrow*pr; j++) {
                cout << myRows[i][j] << "|";
            }
            cout << endl;
        }
        //cin.get();
#endif

        // *********** prava strana **************
        
        if (k >= (col-2))
            continue;

        t1 = MPI_Wtime();
        if ((k%pr) == rank) {
            // radek je na tomto procesu -> spocist a rozeslat vektor v
            for (int j = 0; j < col; j++) {
                if (j < (k+1)) {
                    v[j] = 0.0;
                } else {
                    v[j] = myRows[myRow][j];
                }
            }
            MPI_Bcast(v, col, MPI_DOUBLE, rank, MPI_COMM_WORLD);
        } else {
            // radek je jinde -> prijmout vektor v
            MPI_Bcast(v, col, MPI_DOUBLE, (k%pr), MPI_COMM_WORLD);
        }
        bcastTime += (MPI_Wtime() - t1);

        sum2 = 0.0;
        for (int j = (k+1); j < col; j++) {
            sum2 += v[j]*v[j];
        }

        alfa  = (double)(sign(v[k+1]));
        alfa *= sqrt(sum2);
        beta  = 2*(sum2 + v[k+1]*alfa); // v'*v
        sbeta  = (beta > eps) ? sqrt(beta) : 1.0;
        v[k+1] += alfa;
        
#ifdef UV
		for (int j = 0; j < col; j++) {
            vkV[j] = v[j]/sbeta;
        }
		
        // akumulace householderovych matic zprava na lokalnim kousku matice, deleni je po radich kvuli nulove komunikaci
        acumulateV(myVRows, vkV, myVnRows, col, rank);
#endif

        // deleni nulou je zlo
        if (beta < eps)
            continue;
        // samotna eliminace patricneho radku
        for (int i = start; i < mrow; i++) {
            gamma = 0.0;
            for (int j = (k+1); j < col; j++) {
                gamma += myRows[i][j] * v[j];
            }
            for (int j = (k+1); j < col; j++) {
                myRows[i][j] = myRows[i][j] - (2*gamma/beta)*v[j];
            }
        }

#ifdef debug
        cout << endl << "rank : " << rank << endl;
        cout << "--------" << endl;
        for (int i = 0; i < mrow; i++) {
            for (int j = 0; j < mrow*pr; j++) {
                cout << myRows[i][j] << "|";
            }
            cout << endl;
        }
        //cin.get();
#endif
         
	}// end for (i=0;...) bidiagonalizace
    bidiagTime = (MPI_Wtime() - t0);

    // treba po sobe uklidit bordel a smazat vsecnu alokovanou pamet

    delete[] pCol;
    delete[] tv;
    delete[] v;
    delete[] gama;
    delete[] sgama;
    
// ------------------------------------------------------------

    t0 = MPI_Wtime();
    // kazdy proces si ulozi svou cast dat do souboru na splecny diskovy prostor
    // postprocessing data sesbira a zpracuje
    string myFile = "B_" + itoa(rank) + ".bin";
    double *a = new double[mrow];
    double *b = new double[mrow];
    ofstream f(myFile.data(), ios::out|ios::binary);
    size[0] = mrow; size[1] = 2;
    
    f.write(reinterpret_cast<char *>(&size[0]), 2*sizeof(int));
    for (int i = 0; i < mrow; i++) {
        a[i] = myRows[i][(i*pr)+rank];
        if ((i*pr)+rank < mrow*pr-1)
            b[i] = myRows[i][(i*pr)+rank+1];
        
        // usetrime par cyklu, uklidime po sobe pro jednom, mozna ze bude zlobit tak se to pak presune
        delete[] myRows[i];
    }
    b[mrow-1] = 0.0;
    f.write(reinterpret_cast<char *>(a), mrow*sizeof(double));
    f.write(reinterpret_cast<char *>(b), mrow*sizeof(double));
    f.close();
    delete [] myRows;
    delete [] a; delete [] b;

#ifdef UV
    myFile = "bU_" + itoa(rank) + ".bin";
    writeMatrix(myFile.data(), myUCols, row, myUnCols);
    myFile = "bV_" + itoa(rank) + ".bin";
    writeMatrix(myFile.data(), myVRows, myVnRows, col);
    writeTime = (MPI_Wtime() - t0);
#endif

    cout << "Rank : " << rank << ",readT : " << fileTime << ",bidiagT : " << bidiagTime;
    cout << ",writeT : " << writeTime << endl;
    cout << "Rank : " << rank << ",allGatherT : " << allgatherTime << ",allReduceT : " << allreduceTime;
    cout << ",bcastT : " << bcastTime << endl;

	MPI_Finalize();
	return 0;  
}

string itoa (int n) {
    char * s = new char[17];
    string u;
 
    if (n < 0) { //turns n positive
        n = (-1 * n); 
        u = "-"; //adds '-' on result string
    }
 
    int i = 0; //s counter
  
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

void acumulateU(double **u, double *vk, int m, int n, int rank) {
    
    double t = 0.0;
    double *tsum = new double[m];

    for (int j = 0; j < n; j++) {

        for (int i = 0; i < m; i++)
            tsum[i] = 0.0;

        for (int i = 0; i < m; i++) {
            
            for (int l = 0; l < m; l++) {
                t = ((i==l) ? 1.0 : 0.0) - 2*vk[i]*vk[l];
                tsum[i] += t*u[l][j];
            }

        }

        for (int i = 0; i < m; i++)
            u[i][j] = tsum[i];

    }

}

void acumulateV(double **v, double *vk, int m, int n, int rank) {
    
    double t     = 0.0;
    double *tsum = new double[n];
    int st       = rank*m;

    for (int i = 0; i < m; i++) {
    
        for (int j = 0; j < n; j++)
            tsum[j] = 0.0;

        for (int j = 0; j < n; j++) {
        
            for (int l = 0; l < n; l++) {
                t = ((j==l) ? 1.0 : 0.0) - 2*vk[l]*vk[j];
                tsum[j] += v[i][l]*t;
            }

        }

        for (int j = 0; j < n; j++)
            v[i][j] = tsum[j];

    }


}

bool writeMatrix(const char *file, double **a, int m, int n) {

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
