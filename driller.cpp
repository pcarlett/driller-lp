/**
 * @file giornali.cpp
 * @brief 
 */

#include <cstdio>
#include <iostream>
#include <vector>
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

// funzione per la lettura del file contenente i dati
void read(const char* filename) {
	std::ifstream in(filename);
	// read data from file
	in >> n; // numero di nodi del grafo
	std::cout << "n = " << n << std::endl;
	for (int i = 0; i < n*n; i++) { //per ogni misura
		double c;
		in >> c;	// costo
		C.push_back(c);
		std::cout << "c: " << c << std::endl;
	}
	in.close();
}

// funzione per la creazione dell'ambiente del LP
void setupLP(CEnv env, Prob lp, int & numVars ) {
	
	// aggiungo le variabili y da minimizzare con i costi come da matrice
	//
	//    status =      CPXnewcols (env, lp, ccnt, obj, lb, ub, xctype, colname);
	//
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			char xtype = 'I';
			double lb = 0.0;
			double ub = CPX_INFBOUND;
			snprintf(name, NAME_SIZE, "y_%c_%c", i, j);
			char* xname = (char*)(&name[0]);
			CHECKED_CPX_CALL( CPXnewcols, env, lp, 1   , &C[i+j], &lb, &ub, &xtype, &xname );
			/// status =      CPXnewcols (env, lp, ccnt, obj      , lb  , ub, xctype, colname);
		}
	}
	
	// aggiungo i vincoli del sistema
    //
    //    status = CPXaddrows (env, lp, colcnt, rowcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval , newcolname, newrowname);
    //
	
	
	
	for (int j = 0; j < J; ++j)	{
	  std::vector<int> idx;
	  std::vector<double> coef;
		for (int i = 0; i < I; i++)	{
		  if ( S[i*J+j] < T ) continue;
		  idx.push_back(i);
		  coef.push_back(1.0);		  
		}
		idx.push_back(I+j);
		coef.push_back(-1.0);
		char sense = 'G';
		double rhs = 0.0;
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	// add tolerance constraints
	for (int j = 0; j < J; ++j) {
	  std::vector<int> idx;
	  std::vector<double> coef;
		for (int i = 0; i < I; i++)	{
		  if ( S[i*J+j] < T ) continue;
		  idx.push_back(i);
		  coef.push_back(1.0);		  
		}
		idx.push_back(I+j);
		coef.push_back(M[j]);
		char sense = 'L';
		double rhs = N + M[j];
		int matbeg = 0;
		CHECKED_CPX_CALL( CPXaddrows, env, lp, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0], &coef[0], NULL, NULL );
	}
	
	
	// print (debug);
	CHECKED_CPX_CALL( CPXwriteprob, env, lp, "driller.lp", 0 );
}

// programma principale
int main (int argc, char const *argv[]) {
	try {
		// init
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
		std::cout << "Objval: " << objval << std::endl;
		int n = CPXgetnumcols(env, lp);
		cout << "check num var: " << n << " == " << numVars << endl;
		if (n != numVars) { 
			throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + "different number of variables"); 
		}
		std::vector<double> varVals;
		varVals.resize(n);
		CHECKED_CPX_CALL( CPXgetx, env, lp, &varVals[0], 0, n - 1 );
		for ( int i = 0 ; i < n ; ++i ) {
			std::cout << "var in position " << i << " : " << varVals[i] << std::endl;
			CHECKED_CPX_CALL( CPXsolwrite, env, lp, "driller.sol" );
		}
		
		// free
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);
	} catch(std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}
