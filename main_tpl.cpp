/**
 * @file main.cpp
 * @brief 
 */

#include <iostream>
#include "cpxmacro.h"

using namespace std;

// error status and messagge buffer
int status;
char errmsg[BUF_SIZE];

int main (int argc, char const *argv[])
{
	try
	{
	  cout << "Se non compare alcun meggaggio di errore, l'ambiente e` correttamente configurato per l'uso delle API di Cplex" << endl;
		// TODO
	}
	catch(std::exception& e)
	{
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}
	return 0;
}
