/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <tuple> // per la funzione std::tie
#include <utility> // per la funzione std::make_pair

using namespace std;

//#define M_PI       3.14159265358979323846   // pi

// Calcolo della deviazione standard
double error(const vector<double>& , const vector<double>& , int);

// Funzione che restituisce la media progressiva sui blocchi e relativo errore
void Ave_Block(const vector<double>& , int, string);

// Funzione per stampare un array su file esterno
void Print_File(const vector<double> &, string );

// funzione che calcola la media delle medie sui blocchi
double Ave_Block(const vector<double> &, int);

// funzione che calcola la media al quadrato delle medie sui blocchi
double Ave2_Block(const vector<double> &, int); 

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random(std::string primes, std::string seed);
  Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to save the seed to a file
  void SaveSeed();
  // Method to generate a random number in the range [0,1)
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);

  // Generatore numeri casuali su distribuzione esponenziale
  double Exponential(double lamda);
  // Generatore numeri casuali su distribuzione lorentiziana
  double Lorentz(double gamma, double mu);
  // Generatore di theta casuale, angolo di un versore su S2
  double ThetaVersor();

  // "constructor" of RNG in MPI case
  void RandomMPI(string primes, string seed, int rank);
  // Method to save the seed to a file in MPI case
  void SaveSeedMPI(int rank);

};


#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
