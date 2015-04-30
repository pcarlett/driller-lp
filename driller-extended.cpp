/**
 * @file giornali.cpp
 * @brief 
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include "cpxmacro.h"

using namespace std;

// error status and messagge buffer
int status;
char errmsg[BUF_SIZE];

/*************************/
/* Variabili del sistema */
/*************************/

int n; // # di nodi del grafo
std::vector<double> C; // costi degli archi

const int NAME_SIZE = 512;
char name[NAME_SIZE];

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
void setupLP(CEnv env, Prob lp, int & numVars) {
	
	try {
		// aggiungo le variabili y da minimizzare con i costi come da matrice
		//
		//    status =      CPXnewcols (env, lp, ccnt, obj, lb, ub, xctype, colname);

		//	min sum(c_ij * y_ij)	
	
		cout << "1) Funzione da minimizzare:" << endl;
		cout << "min ";
		// adding y vars
		for (int i = 0; i < n; i++)	{
			for (int j = 0; j < n; j++)	{
			
				char xtype = 'B';
				double lb = 0.0;
				double ub = 1.0;
				int l = i+1;
				int r = j+1;
				snprintf(name, NAME_SIZE, "y_%i_%i", l, r);
				char* xname = (char*)(&name[0]);
			
				if(C[i*n+j]) {
					CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &C[i*n+j], &lb, &ub, &xtype, &xname);
					cout << "+ " << C[i*n+j] << " y_" << i+1 << j+1 << " ";
				}
			}
		}
		cout << endl << endl;

		// adding x vars
		for (int i = 0; i < n; i++)	{
			for (int j = 0; j < n; j++)	{
			
				char xtype = 'I';
				double obj = 0.0; // variabile non presente nella funzione obiettivo
				double lb = 0.0;
				double ub = CPX_INFBOUND;
				int l = i+1;
				int r = j+1;
				snprintf(name, NAME_SIZE, "x_%i_%i", l, r);
				char* xname = (char*)(&name[0]);
			
				if(C[i*n+j]) {
					CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &xtype, &xname);
					cout << "+ x_" << i+1 << j+1 << " ";
				}
			}
		}	
		cout << endl << endl;
		
		numVars = CPXgetnumcols(env, lp);
		
		cout << "Num Vars: " << numVars << endl;
	
		// aggiungo i vincoli del sistema
	    //
	    //    status = CPXaddrows (env, lp, colcnt, rowcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval , newcolname, newrowname);
		//
		//		CHECKED_CPX_CALL( CPXaddrows, env, lp, colcnt, rowcnt, nzcnt, &rhs[0], &sense[0], &rmatbeg[0], &rmatind[0], &rmatval[0], newcolnames , rownames );
		//
		// int colcnt = 0; // nessuna nuova colonna
		// int rowcnt = 0; // iteratore per il conteggio delle colonne
		// int nzcnt = 0; // iteratore per il conteggio dei valori diversi da 0
		double rhs; // right-hand-side
		char sense; // constraint type 'L' or 'E' or 'G' ...
		int matbeg = 0; // one element for each constraint (row), reporting the index where each row of the coefficient matrix starts
		vector<int> rmatind; // one element for each zon-zero ceofficient, reporting its column index
		vector<double> rmatval; // the linearized vector of non-zero coefficients
		// char ** newcolnames = NULL; // no names
		// char ** rownames = NULL; // no names
	
		//	2) x_0j = |N|
		cout << "2) Vincolo:	x_0j = |N|" << endl;
		for(int i = 0; i < 1; i++) {
			for(int j = 0; j < n; j++) {
				if(C[i*n+j]) {
					rmatind.push_back(n*n+i*n+j);
					rmatval.push_back(1);
					cout << "+ x_" << i+1 << j+1 << " ";
				}
			}
			matbeg = 0;
			rhs = n;
			sense = 'E';
			cout << " = " << n << ";" << endl;
			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , rmatind.size(), &rhs,  &sense, &matbeg, &rmatind[0], &rmatval[0], NULL      , NULL      );
		    /// status =      CPXaddrows (env, lp, colcnt, rowcnt, nzcnt     	 , rhs  , sense , rmatbeg, rmatind, 	rmatval , 	 newcolname, newrowname);

		}
		cout << endl;
	
		//	3) sum(x_ik) - sum(x_kj) = 1	forall k in N
		cout << "3) Vincolo:	sum(x_ik) - sum(x_kj) = 1" << endl;
		for(int i = 0; i < n; i++) {
		
			// Clear matrix vectors
			rmatind.clear();
			rmatval.clear();
		
			for(int j = 0; j < n; j++) {
				if(C[i*n+j]) {
					rmatind.push_back(n*n+i*n+j);
					rmatval.push_back(-1);
					cout << "- x_" << i+1 << j+1 << " ";
				}
				if(C[i+j*n]) {
					rmatind.push_back(n*n+i*n+j);
					rmatval.push_back(1);
					cout << "+ x_" << j+1 << i+1 << " ";
				}

			}
			matbeg = 0;
			rhs = 1;
			sense = 'E';
			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , rmatind.size(), &rhs,  &sense, &matbeg, &rmatind[0], &rmatval[0], NULL      , NULL      );
		    /// status =      CPXaddrows (env, lp, colcnt, rowcnt, nzcnt     	 , rhs  , sense , rmatbeg, rmatind, 	rmatval , 	 newcolname, newrowname);
			cout << " = 1" << endl;
		}
		cout << endl;

		//	4) sum (y_ij) = 1	forall i in N
		cout << "4) Vincolo:	sum(y_ij) = 1	(for i in N)" << endl;
		for(int i = 0; i < n; i++) {
		
			// Clear matrix vectors
			rmatind.clear();
			rmatval.clear();
		
			for(int j = 0; j < n; j++) {
				if(C[i*n+j]) {
					rmatind.push_back(i*n+j);
					rmatval.push_back(1);
					cout << "+ y_" << i+1 << j+1 << " ";
				}

			}
			matbeg = 0;
			rhs = 1;
			sense = 'E';
			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , rmatind.size(), &rhs,  &sense, &matbeg, &rmatind[0], &rmatval[0], NULL      , NULL      );
		    /// status =      CPXaddrows (env, lp, colcnt, rowcnt, nzcnt     	 , rhs  , sense , rmatbeg, rmatind, 	rmatval , 	 newcolname, newrowname);
			cout << " = 1" << endl;
		}
		cout << endl;

		// 	5) sum (y_ij) = 1	forall j in N
		cout << "5) Vincolo:	sum(y_ij) = 1	(for j in N)" << endl;
		for(int i = 0; i < n; i++) {
		
			// Clear matrix vectors
			rmatind.clear();
			rmatval.clear();
		
			for(int j = 0; j < n; j++) {
				if(C[i+j*n]) {
					rmatind.push_back(j*n+i);
					rmatval.push_back(1);
					cout << "+ y_" << j+1 << i+1 << " ";
				}

			}
			matbeg = 0;
			rhs = 1;
			sense = 'E';
			CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , rmatind.size(), &rhs,  &sense, &matbeg, &rmatind[0], &rmatval[0], NULL      , NULL      );
		    /// status =      CPXaddrows (env, lp, colcnt, rowcnt, nzcnt     	 , rhs  , sense , rmatbeg, rmatind, 	rmatval , 	 newcolname, newrowname);
			cout << " = 1" << endl;
		}
		cout << endl;

		//	6) x_ij <= |N| * y_ij
		cout << "6) Vincolo:	x_ij <= |N| * y_ij	(for (i,j) in A)" << endl;
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
			
				// Clear matrix vectors
				rmatind.clear();
				rmatval.clear();
			
				if(C[i*n+j]) {
					// variabile y
					rmatind.push_back(j*n+i);
					rmatval.push_back(-n);
					// variabile x
					rmatind.push_back(n*n+j*n+i);
					rmatval.push_back(1);
				
					matbeg = 0;
					rhs = 0;
					sense = 'L';
					CHECKED_CPX_CALL( CPXaddrows, env, lp, 0     , 1     , rmatind.size(), &rhs,  &sense, &matbeg, &rmatind[0], &rmatval[0], NULL      , NULL      );
				    /// status =      CPXaddrows (env, lp, colcnt, rowcnt, nzcnt     	 , rhs  , sense , rmatbeg, rmatind, 	rmatval , 	 newcolname, newrowname);
			
					cout << "x_" << i+1 << j+1 << " <= " << n << " * y_" << i+1 << j+1 << endl;
				}
			}
		}
		cout << endl;
	
		// print (debug)
		CHECKED_CPX_CALL( CPXwriteprob, env, lp, "driller.lp", NULL );
		/// status =      CPXwriteprob (env, lp, "myprob"    , filetype_str);

	} catch(std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
}

// programma principale
int main (int argc, char const *argv[]) {
	try {
		///////////////////////// init
		DECL_ENV( env );
		DECL_PROB( env, lp );
		
		// legge le variabili dal file di supporto
		read(argv[1]);
		
		// setup LP
		int numVars;
		setupLP(env, lp, numVars);
		
		// optimize
		CHECKED_CPX_CALL( CPXmipopt, env, lp );
		
		// print
		double objval;
		CHECKED_CPX_CALL( CPXgetobjval, env, lp, &objval );
		std::cout << "******************************" << endl << "Valore Obiettivo: " << objval << std::endl;
		/*int n = CPXgetnumcols(env, lp);
		if (n != 2*I*J+1) { 
			throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + "different number of variables"); 
		}*/
		std::vector<double> varVals;
		varVals.resize(n*n);
		CHECKED_CPX_CALL( CPXgetx, env, lp, &varVals[0], 0, n*n - 1 );
		/// status =      CPXgetx (env, lp, x          , 0, CPXgetnumcols(env, lp)-1);
  		for ( int i = 0 ; i < n ; ++i ) {
	  		for ( int j = 0 ; j < n ; ++j ) {
				if(varVals[i*n+j]) {
					std::cout << "Variabile y_" << i+1 << "_" << j+1 << ": " << varVals[i*n+j] << std::endl;
					/// per leggere i nomi , usare la funzione CPXgetcolname (un po' complicata)
					/// status = CPXgetcolname (env, lp, cur_colname, cur_colnamestore, cur_storespace, &surplus, 0, cur_numcols-1);
				}
			}
	  	}
		CHECKED_CPX_CALL( CPXsolwrite, env, lp, "driller.sol" );
		
		// free
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
	} catch(std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}

