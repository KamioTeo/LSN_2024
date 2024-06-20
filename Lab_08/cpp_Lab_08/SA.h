#ifndef __SA__
#define __SA__

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

double TMin, TMax, TStep; // Temperature range
int Nsteps; // step della temperatura
double mu, sigma; // valori iniziali di mu e sigma

double deltaSA; // delta per campionamento |psi|^2 e Simulated Annealing
int wd=14; // width, larghezza campo output

bool SA; // se runnare l'algortimo SA o no

//functions
double psi_trial(double, double, double); // F. d'onda di prova
double potential(double); // Potenziale confinante
double eigen_val(double, double, double); // calcolo di (H*psi)/psi

double Metropolis_delta(double, double, int, double&);
vector<double> Metropolis(double, double, double, int, double&);

void Input(void);

#endif
