// ESERCIZIO 5: Hydrogen Wave function Metropolis sampling
// 5.1 Transizione Uniforme

#include "../RandomNumberGenerator/random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

// Scrivo tutto in unità di raggio di Bohr (a_0=1)
// Pdf Ground State, modulo quadro della funzione d'onda
double pdf_gs(double x, double y, double z) {
  return exp(-2*sqrt(x*x+y*y+z*z))/M_PI;
}

// Pdf stato eccitato 2p, modulo quadro della funzione d'onda
double pdf_2p(double x, double y, double z) {
  return  z*z*exp(-sqrt(x*x+y*y+z*z))/(32.*M_PI);
} // OSS: cos(theta) = z / |r|

// Funzione che esegue uno step dell'algoritmo di Metropolis
void Move_Metropolis(double (*pdf)(double, double, double), double &, double &, double &, double, int &, Random &);

// Funzione che restituisce un array di posizioni che segue la pdf data in ingresso, utilizzando l'algoritmo di metropolis a partire dalla posizione iniziale voluta 
std::vector<double> Metropolis(double (*pdf)(double, double, double), double &, double &, double &, int, double, int &, Random &);

int main(int argc, char *argv[]) {
	
  // Generazione del random number generator
  Random rnd("../RandomNumberGenerator/Primes",
			 "../RandomNumberGenerator/seed.in");

  // Definizioni variabili
  bool print_eq = true; // se vero stampa su file esterno i valori durante l'equilibrazione, altrimenti equilibra e basta
  const int M_eq = 1E4;   // Numero di step random walk per equilibrare
  
  // Definisco i parametri delta per cui l'accettazione è circa il 50%. La probabilità di transizione T tra due vettori posizione è uniforme in [x_n-del;x_n+del], [y_n-del;y_n+del] e [z_n-del;z_n+del]
  const double del_gs = 1.223; // delta per il Ground State
  const double del_2p = 2.98;  // delta per il caso eccitato
  
  const int M = 1E6; // Numero di step dopo l'equilibrazione
  const int N = 100; // Numero di blocchi

  // definisco il contatore di mosse accettate
  int count_gs = 0; // caso gs
  int count_2p = 0; // caso 2p
	
  // definizione coordinate iniziali per GS (punto più probabile)
  double x_gs = 0.;
  double y_gs = 0.;
  double z_gs = 0.;
  // definizione coordinate iniziali per 2p (punto più probabile)
  double x_2p = 0.;
  double y_2p = 0.;
  double z_2p = 2.;

  // Equilibrazione, salvo le posizioni che seguono le pdf; in questa fase il valore di delta e della frazione di accettazione non sono fondamentali
  // Caso Ground State
  cout << "--- Equilibrating GS... ---" << endl;
  std::vector<double> Positions_eq_gs = Metropolis(pdf_gs, x_gs, y_gs, z_gs, M_eq, del_gs, count_gs, rnd);
  // azzero il contatore perché conto l'accettazione solo dopo l'equilibrazione
  count_gs = 0;
  cout << "--- Equilibrated GS! ---" << endl << endl;
  
  // Caso eccitato 2p
  cout << "--- Equilibrating 2P... ---" << endl;
  std::vector<double> Positions_eq_2p = Metropolis(pdf_2p, x_2p, y_2p, z_2p, M_eq, del_2p, count_2p, rnd);
  // azzero il contatore perché conto l'accettazione solo dopo l'equilibrazione
  count_2p = 0;
  cout << "--- Equilibrated 2P! ---" << endl << endl;

  // se la condizione è verificata, stampo le posizioni su file esterno dive in blocchi
  if(print_eq) {
    Ave_Block(Positions_eq_gs, N, "dati_unif/Posizioni_eq_gs_unif.txt");
    Ave_Block(Positions_eq_2p, N, "dati_unif/Posizioni_eq_2p_unif.txt");
    cout << "--- Printed equilibration data ---" << endl << endl;
  } 

 // Usando la posizione finale raggiunta dopo l'equilibrazione, procedo con le M simulazioni. Qui delta deve dare accettazione di circa il 50%

  cout << endl << "GS last position: (" << x_gs << ";" << y_gs << ";" << z_gs << ")" << endl << "2P last position: (" << x_2p << ";" << y_2p << ";" << z_2p << ")" << endl << endl;
  
  // Caso Ground State
  cout << "--- Metropoling GS... ---" << endl;
  std::vector<double> Positions_gs = Metropolis(pdf_gs, x_gs, y_gs, z_gs, M, del_gs, count_gs, rnd);
  cout << "--- Motropoled GS! ---" << endl << endl;

  // Caso eccitato 2p
  cout << "--- Metropoling 2P... ---" << endl;
  std::vector<double> Positions_2p = Metropolis(pdf_2p, x_2p, y_2p, z_2p, M, del_2p, count_2p, rnd);
  cout << "--- Metropoled 2P! ---" << endl << endl;

  // Stampo i risultati su file esterno, dividendo per blocchi
  Ave_Block(Positions_gs, N, "dati_unif/Posizioni_gs_unif.txt");
  Ave_Block(Positions_2p, N, "dati_unif/Posizioni_2p_unif.txt");

  cout << "--- Printed data ---" << endl;
  cout << "--- End of code ---" << endl;
  
  rnd.SaveSeed();
  return 0;
}

//========================================================================
// DEFINIZIONE FUNZIONI
//========================================================================

// Funzione Metropolis: restituisce un array con posizioni che seguono la distribuzione pdf in ingresso, a partire da 3 coordinate iniziali
std::vector<double> Metropolis(double (*pdf)(double, double, double), double &x, double &y, double &z, int Nsim, double del, int & count, Random & rnd) {

  vector<double> Position(Nsim);
	
  for (int i = 0; i < Nsim; i++) { // Per ogni simulazione
	  Position[i] = sqrt(x*x+y*y+z*z);
    Move_Metropolis(pdf, x, y, z, del, count, rnd);
  } // fine ciclo / step del random walk

  cout << "Metropolis called with " << Nsim << " values" << endl << "acceptance ratio: " << fixed << setprecision(3) << double(count)/Nsim << endl;
  return Position;
}

//========================================================================
// funzione che esegue uno step di algoritmo di Metropolis, cambiando per riferimento i valori di posizione in ingresso, in modo da poter salvare l'ultima raggiunta dopo equilibrazione
void Move_Metropolis(double (*pdf)(double, double, double), double & x, double & y, double & z, double del, int & count, Random & rnd) {

      double x_new = rnd.Rannyu(x-del, x+del); // posizione al passo successivo, usando la probabilità di Transizione
      double y_new = rnd.Rannyu(y-del, y+del);
      double z_new = rnd.Rannyu(z-del, z+del);

  // dato che T è simmetrica, la probabilità di accettazione  è data solo dal rapporto tra la probabilità da campionare valutata in x_new e x_old
      double A = pdf(x_new, y_new, z_new)/pdf(x,y,z);

      // Prendo il minimo tra A e 1
      if( A>1 ) { A=1; }

      double r = rnd.Rannyu();
      if( r <= A ) { // se è vero accetto la nuova posizione, quindi x diventa x_new
       x = x_new;
       y = y_new;
       z = z_new;
       count += 1; // Conto il fatto che ho accettato
      }
      // Se non accetto, la nuova posizione resta quella vecchia, quindi sovrascrivo x_new con un altro tentativo a partire sempre da x (lo fa in automatico all'inizio del nuovo ciclo)
}
