/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

#include <vector>
#include <string>
//Random numbers
#include "random.h"

using namespace std;

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000; // Numero massimo di variabili fisiche misurabili
int n_props, iu, ic, im, ix, ig, iu2, im2;
double nbins;
double walker[m_props]; // array che contiene i valori misurati ad ogni step delle grandezze fisiche

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_u, stima_c, stima_m, stima_x, stima_g, stima_u2, stima_m2;
double err_u, err_c, err_m, err_x, err_g;

//configuration
const int m_spin=50; // numero massimo di spin caricabili
double s[m_spin];

// thermodynamical state
int nspin; // numero di spin effettivamente caricati
double beta,temp,J,h; // proprietà fisiche del sistema caricate

// simulation
int nstep, nblk, metro, restart, Tsequence, Equilibration;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void LastAverages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);

// Funzione per stampare un array su file esterno
void Print_File(const vector<double> &, string );
void EnergyAverage(int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
