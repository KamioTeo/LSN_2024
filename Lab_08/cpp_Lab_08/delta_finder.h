#ifndef __DELTA__
#define __DELTA__

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
//Random numbers
#include "../RandomNumberGenerator/random.h"

using namespace std;

// Generazione del random number generator
Random rnd("../RandomNumberGenerator/Primes", "../RandomNumberGenerator/seed.in");

// simulations
int neq, M, N; // step per equilibrare, step dopo l'equilibrazione, numero di blocchi

double muMin, muMax, muStep; // range for mu
double sigmaMin, sigmaMax, sigmaStep; // range for sigma

int iterations = 0; // tentativi fatti per trovare delta
double acceptance = 0.; // accepted/attempted
int wd=14; // width, larghezza campo output

//functions
double psi_trial(double, double, double); // F. d'onda di prova
double potential(double); // Potenziale confinante
double eigen_val(double, double, double); // calcolo di (H*psi)/psi

vector<double> Metropolis(double, double, double, int, double&);
double Metropolis_delta(double, double, int, double&, double&, int&);

void Input(void);

#endif
