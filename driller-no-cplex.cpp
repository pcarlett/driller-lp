/**
 * @file giornali.cpp
 * @brief 
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
// #include "cpxmacro.h"

using namespace std;

// error status and messagge buffer
// int status;
// char errmsg[BUF_SIZE];

/*************************/
/* Variabili del sistema */
/*************************/

int n; // # di nodi del grafo
std::vector<double> C; // costi degli archi

// funzione per la lettura del file contenente i dati
void read(const char* filename) {
	std::ifstream in(filename);
	// read data from file
	in >> n; // numero di nodi del grafo
	std::cout << "n = " << n << std::endl;
	for (int i = 0; i < n; i++) { //per ogni misura
		for (int j = 0; j < n; j++) {
			double c;
			in >> c;	// costo
			C.push_back(c);
			std::cout << c << "\t";
		}
		cout << endl;
	}
	cout << endl;
	in.close();
}

// funzione per la creazione dell'ambiente del LP
void setupLP( ) {
	
	// aggiungo le variabili y da minimizzare con i costi come da matrice
	//
	//    status =      CPXnewcols (env, lp, ccnt, obj, lb, ub, xctype, colname);

	//	min sum(c_ij * y_ij)	
	cout << "1) Funzione da minimizzare:" << endl;
	cout << "min ";
	for (int i = 0; i < n; i++)	{
		for (int j = 0; j < n; j++)	{
			if(C[i*5+j]) {
				cout << C[i*5+j] << "+ y_" << i+1 << j+1 << " ";
			}
		}
	}
	cout << endl << endl;

	// aggiungo i vincoli del sistema
    //
    //    status = CPXaddrows (env, lp, colcnt, rowcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval , newcolname, newrowname);
	
	//	2) x_0j = |N|
	cout << "2) Vincolo:	x_0j = |N|" << endl;
	for(int i = 0; i < 1; i++) {
		for(int j = 0; j < n; j++) {
			if(C[i*5+j]) {
				cout << "+ x_" << i+1 << j+1 << " ";
			}

		}
		cout << " = " << n << ";" << endl;
	}
	cout << endl;

	//	3) sum(x_ik) - sum(x_kj) = 1	forall k in N
	cout << "3) Vincolo:	sum(x_ik) - sum(x_kj) = 1" << endl;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(C[i*5+j]) {
				cout << "- x_" << i+1 << j+1 << " ";
			}
			if(C[i+j*5]) {
				cout << "+ x_" << j+1 << i+1 << " ";
			}

		}
		cout << " = 1" << endl;
	}
	cout << endl;

	//	4) sum (y_ij) = 1	forall i in N
	cout << "4) Vincolo:	sum(y_ij) = 1	(for i in N)" << endl;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(C[i*5+j]) {
				cout << "+ y_" << i+1 << j+1 << " ";
			}

		}
		cout << " = 1" << endl;
	}
	cout << endl;

	// 	5) sum (y_ij) = 1	forall j in N
	cout << "5) Vincolo:	sum(y_ij) = 1	(for j in N)" << endl;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(C[i+j*5]) {
				cout << "+ y_" << j+1 << i+1 << " ";
			}

		}
		cout << " = 1" << endl;
	}
	cout << endl;

	//	6) x_ij <= |N| * y_ij
	cout << "6) Vincolo:	x_ij <= |N| * y_ij	(for (i,j) in A)" << endl;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(C[i*5+j]) {
				cout << "x_" << i+1 << j+1 << " <= " << n << " * y_" << i+1 << j+1 << endl;
			}

		}
	}
	cout << endl;
}

// programma principale
int main (int argc, char const *argv[]) {
	try {
		
		// legge le variabili dal file di supporto
		read(argv[1]);
		
		// setup LP
		int numVars;
		setupLP();

	} catch(std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}

