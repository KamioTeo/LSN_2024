/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

// funzione per calcolare la deviazione standard a partire dalle somme cumulate
double error(const vector<double> &sum_prog, const vector<double> &su2_prog, int i) {
  if (i == 0)
    return 0;
  double var = su2_prog[i] - pow(sum_prog[i], 2);
  if (var < 0) {
    cerr << "Errore: varianza negativa" << endl;
    return 0;
  }
  return sqrt(var / double(i));
}

//================================================//

// funzione che stampa un file con le medie sui blocchi e relativo errore. In ingresso gli do la lista delle simulazioni e il numero di blocchi in cui dividere
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

//================================================//
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
//================================================//
// funzione che restituisce la media a N blocchi
double Ave_Block(const vector<double> &sim_arr, int N) {  
	int M = sim_arr.size(); // numero di simulazioni
	int L = M / N;          // numero di lanci per blocco

  	vector<double> ave_block(N); // vettore che contiene le medie su ogni blocco
	double average = 0;  // media delle medie sui blocchi

	for (int i = 0; i < N; i++) { // per ogni blocco calcolo la media sul blocco
    	double sum_f = 0;
    	for (int j = 0; j < L; j++) {  // per ogni numero nel blocco
      		int k = j + i * L;  // k-esimo numero casuale dallo 0-esimo nello 0-esimo blocco ecc
      		sum_f += sim_arr[k];
    	}
    	ave_block[i] = sum_f / L;  // valor medio nel blocco i
  }

	// ora calcolo la media su tutti i blocchi
	for (int i = 0; i < N; i++) {
		average += ave_block[i];
  	}
	return average/N;
}

//================================================//
// funzione che restituisce la media dei quadrati a N blocchi
double Ave2_Block(const vector<double> &sim_arr, int N) {  
	int M = sim_arr.size(); // numero di simulazioni
	int L = M / N;          // numero di lanci per blocco

  	vector<double> ave2_block(N); // vettore che contiene le medie al quadrato su ogni blocco
	double average2 = 0;  // media al qudrato delle medie sui blocchi

	for (int i = 0; i < N; i++) { // per ogni blocco calcolo la media al quadrato sul blocco
    	double sum_f = 0;
    	for (int j = 0; j < L; j++) {  // per ogni numero nel blocco
      		int k = j + i * L;  // k-esimo numero casuale dallo 0-esimo nello 0-esimo blocco ecc
      		sum_f += sim_arr[k];
    	}
    	ave2_block[i] = pow(sum_f / L, 2);  // valor medio al quadrato nel blocco i
  }

	// ora calcolo la media su tutti i blocchi
	for (int i = 0; i < N; i++) {
		average2 += ave2_block[i];
  	}
	return average2/N;
}

//================================================//

Random ::Random() {}
// Default constructor, does not perform any action

Random ::Random(string file_primes, string file_seed) {

  int seed[4]; // Array che contiene i 4 seed di partenza
  int n1, n2;
  ifstream primes(file_primes);
  if (primes.is_open()) {
    primes >> n1 >> n2;
  } else {
    cerr << "Error: unable to open " << file_primes << endl;
    exit(1);
  }
  primes.close();

  ifstream start(file_seed);
  string property;
  if (start.is_open()) {
    while (!start.eof()) {
      start >> property;
      if (property == "RANDOMSEED") {
        start >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        SetRandom(seed, n1, n2); // setta dei parametri interni della classe random per generare i numeri casuali
      }
    }
    start.close();
  } else {
    cerr << " Error: unable to open " << file_seed << endl;
    exit(1);
  }
}

Random ::~Random() {}
// Default destructor, does not perform any action

void Random ::SaveSeed() {
  // This function saves the current state of the random number generator to a
  // file "seed.out"
  ofstream WriteSeed;
  WriteSeed.open("seed.out");
  if (WriteSeed.is_open()) {
    WriteSeed << "RANDOMSEED " << l1 << " " << l2 << " " << l3 << " " << l4
              << endl;
    ;
  } else
    cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random ::Gauss(double mean, double sigma) {
  // This function generates a random number from a Gaussian distribution with
  // given mean and sigma
  double s = Rannyu();
  double t = Rannyu();
  double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
  return mean + x * sigma;
}

double Random ::Rannyu(double min, double max) {
  // This function generates a random number in the range [min, max]
  return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void) {
  // This function generates a random number in the range [0,1)
  const double twom12 = 0.000244140625;
  int i1, i2, i3, i4;
  double r;

  i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
  i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
  i3 = l3 * m4 + l4 * m3 + n3;
  i4 = l4 * m4 + n4;
  l4 = i4 % 4096;
  i3 = i3 + i4 / 4096;
  l3 = i3 % 4096;
  i2 = i2 + i3 / 4096;
  l2 = i2 % 4096;
  l1 = (i1 + i2 / 4096) % 4096;
  r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

  return r;
}

void Random ::SetRandom(int *s, int p1, int p2) {
  // This function sets the seed and parameters of the random number generator
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

double Random ::Exponential(double lambda) {
  // This function generates a random number in the range [0, +inf)
  return -log(1 - Rannyu()) / lambda;
}

double Random :: Lorentz(double gamma, double mu) {
  // This function generates a random number in the range (-inf, +inf)
  return gamma * tan(M_PI * (Rannyu() - 0.5)) + mu;
}

double Random :: ThetaVersor() {
	// This function generates random theta angle of a versor on the unit sphere (theta in [0,PI])
	// phi è uniforme tra 0 e 2PI : Rannyu(0, 2*M_PI)
	return acos(1-2*Rannyu());
}

//===========================================
// per ogni rank carico una coppia diversa di numeri primi
void Random ::RandomMPI(string file_primes, string file_seed, int rank) {

  int seed[4]; // Array che contiene i 4 seed di partenza
  int n1, n2;
  ifstream primes(file_primes);
  if (primes.is_open()) {
	// leggo a vuoto i numeri primi per un multiplo di 10 del rank. Carico la riga 1, 11, 21 ecc
	int void1, void2;
      for(int i=0; i < 10*rank; i++) {
         primes >> void1 >> void2;
	  }
	// carico i numeri primi
    primes >> n1 >> n2;
  } else {
    cerr << "Error: unable to open " << file_primes << endl;
    exit(1);
  }
  primes.close();

  ifstream start(file_seed);
  string property;
  if (start.is_open()) {
    while (!start.eof()) {
      start >> property;
      if (property == "RANDOMSEED") {
        start >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        SetRandom(seed, n1, n2); // setta dei parametri interni della classe random per generare i numeri casuali
      }
    }
    start.close();
  } else {
    cerr << " Error: unable to open " << file_seed << endl;
    exit(1);
  }
}


void Random ::SaveSeedMPI(int rank) {
  // This function saves the current state of the random number generator to a
  // file "seed.out"
  ofstream WriteSeed;
  string rank_str = to_string(rank);
  WriteSeed.open("seed_"+rank_str+".out");
  if (WriteSeed.is_open()) {
    WriteSeed << "RANDOMSEED " << l1 << " " << l2 << " " << l3 << " " << l4
              << endl;
    ;
  } else
    cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}
//===========================================

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
