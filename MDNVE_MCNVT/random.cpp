/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include <vector>

using namespace std;


// funzione per calcolare la deviazione standard a partire dalle somme cumulate
double error(const vector<double> &sum_prog, const vector<double> &su2_prog,
             int i) {
  if (i == 0)
    return 0;
  double var = su2_prog[i] - pow(sum_prog[i], 2);
  if (var < 0) {
    cerr << "Errore: varianza negativa" << endl;
    return 0;
  }
  return sqrt(var / double(i));
}

// funzione che restituisce la media progressiva sui blocchi e relativo errore, in ingresso gli do la lista delle simulazioni e il numero di blocchi in cui
// dividere
void Ave_Block(const vector<double> &sim_arr, int N, string File_name) {
  int M = sim_arr.size(); // numero di simulazioni
  int L = M / N;          // numero di lanci per blocco

  vector<double> ave(N); // lista che conterrà la media cumulativa delle simulazioni al variare del numero N di blocchi utilizzato
  vector<double> av2(N); // analogo ma contiene la media cumulativa al quadrato
  vector<double> sum_prog(N); // analogo ma è la somma cumulativa dei numeri casuali, per calcolo della dev std
  vector<double> su2_prog(N); // somma cumulativa dei quadrati dei numeri casuali, per calcolo della dev std
  vector<double> err_prog(N); // dev std cumulativa

  for (int i = 0; i < N; i++) { // per ogni blocco calcolo la media sul blocco
    double sum_f = 0;
    for (int j = 0; j < L; j++) { // per ogni numero casuale nel blocco
      int k = j + i * L; // k-esimo numero casuale dallo 0-esimo nello 0-esimo blocco ecc
      sum_f += sim_arr[k];
    }
    ave[i] = sum_f / L;      // valor medio nel blocco i
    av2[i] = pow(ave[i], 2); // (r_i)^2 quadrato della media
  }

  // Genero un file esterno in cui verranno salvate le medie cumulative e relativo errore di numeri pseudocasuali generati con distribuzione uniforme. Questi sono divisi in N blocchi
  ofstream out;
  out.open(File_name);

  // Controllo la corretta apertura del file
  if (!out.is_open()) {
    cerr << "Error: unable to open " + File_name << endl;
    exit(1);
  }

  // Prima riga del file esterno
  out << "Media" << " Errore" << endl;

  // ora calcolo la media cumulativa, sommando in progressione il contributo di ogni blocco
  for (int i = 0; i < N; i++) { // per ogni blocco
    for (int j = 0; j <= i; j++) {
      sum_prog[i] += ave[j]; 	// SUM_{j=0,i} pi_j , sommo tutte le medie stimate fino al blocco i+1-esimo
      su2_prog[i] += av2[j]; 	// SUM_{j=0,i} (pi_j)^2, sommo tutte le medie al quadrato stimate fino al blocco i+1-esimo
    }
    sum_prog[i] /= (i + 1); 	// calcolo la media cumulativa
    su2_prog[i] /= (i + 1); 	// media dei quadrati cumulativa
    err_prog[i] = error(sum_prog, su2_prog,i); // Deviazione standard calcolata fino al blocco i-esimo

    out << sum_prog[i] << " " << err_prog[i] << endl; // stampo i risultati
  }

  out.close();
}

//============================================================
// funzione che stampa su file esterno un array di dati
void Print_File(const vector<double> & arr, string File_name) {
  ofstream out;
  out.open(File_name);

  // Controllo la corretta apertura del file
  if (!out.is_open()) {
    cerr << "Error: unable to open " + File_name << endl;
    exit(1);
  }

	int dim = arr.size();
	for(int i=0; i < dim; i++){
		out << arr[i] << endl;
	}
	
  out.close();
	
}

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){ // Scrive il seed usato su file esterno
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){ // setta dei parametri interni della classe random per generare i numeri casuali
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
