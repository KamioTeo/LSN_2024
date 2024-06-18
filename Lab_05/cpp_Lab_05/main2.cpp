// ESERCIZIO 5: Hydrogen Wave function Metropolis sampling
// 5.2 Transizione Gaussiana

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
void Move_Metropolis(double (*pdf)(double, double, double), double &, double &, double &, double, int &, Random &, std::vector<float> &, bool);

// Funzione che restituisce un array di posizioni che segue la pdf data in ingresso, utilizzando l'algoritmo di metropolis a partire dalla posizione iniziale voluta. Salva anche le coordinate estratte
std::vector<double> Metropolis(double (*pdf)(double, double, double), double &, double &, double &, int, double, int &, Random &, std::vector<std::vector<float>> &, bool, int);

int main(int argc, char *argv[]) {
	
  // Generazione del random number generator
  Random rnd("../RandomNumberGenerator/Primes",
			 "../RandomNumberGenerator/seed.in");

  // Definizioni variabili
  bool print_eq = true; // se vero stampa su file esterno i valori durante l'equilibrazione, altrimenti equilibra e basta
  const int M_eq = 1E4;   // Numero di step random walk per equilibrare
  
  // Definisco i parametri delta per cui l'accettazione è circa il 50%. La probabilità di transizione T tra due vettori posizione è ora gaussiana
  const double del_gs = 0.873; // delta per il Ground State
  const double del_2p = 1.371;  // delta per il caso eccitato
  
  const int M = 1E6; // Numero di step dopo l'equilibrazione
  const int N = 100; // Numero di blocchi

  const int N_coord = 50000; // Numero di coordinate da salvare
  
  std::vector<std::vector<float>> coord_2p(N_coord, std::vector<float>(3)); // Dove salverò le coordinate x,y,z delle posizioni
  
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
  std::vector<double> Positions_eq_gs = Metropolis(pdf_gs, x_gs, y_gs, z_gs, M_eq, del_gs, count_gs, rnd, coord_2p, false, N_coord);
  // azzero il contatore perché conto l'accettazione solo dopo l'equilibrazione. Qui non salvo le coordinate (false in input)
  count_gs = 0;
  cout << "--- Equilibrated GS! ---" << endl << endl;
  
  // Caso eccitato 2p
  cout << "--- Equilibrating 2P... ---" << endl;
  std::vector<double> Positions_eq_2p = Metropolis(pdf_2p, x_2p, y_2p, z_2p, M_eq, del_2p, count_2p, rnd, coord_2p, false, N_coord);
  // azzero il contatore perché conto l'accettazione solo dopo l'equilibrazione
  count_2p = 0;
  cout << "--- Equilibrated 2P! ---" << endl << endl;

  // se la condizione è verificata, stampo le posizioni su file esterno dive in blocchi
  if(print_eq) {
    Ave_Block(Positions_eq_gs, N, "dati_gauss/Posizioni_eq_gs_gauss.txt");
    Ave_Block(Positions_eq_2p, N, "dati_gauss/Posizioni_eq_2p_gauss.txt");
    cout << "--- Printed equilibration data ---" << endl << endl;
  } 

 // Usando la posizione finale raggiunta dopo l'equilibrazione, procedo con le M simulazioni. Qui delta deve dare accettazione di circa il 50%

  cout << endl << "GS last position: (" << x_gs << ";" << y_gs << ";" << z_gs << ")" << endl << "2P last position: (" << x_2p << ";" << y_2p << ";" << z_2p << ")" << endl << endl;
  
  // Caso Ground State
  cout << "--- Metropoling GS... ---" << endl;
  std::vector<double> Positions_gs = Metropolis(pdf_gs, x_gs, y_gs, z_gs, M, del_gs, count_gs, rnd, coord_2p, false, N_coord);
  cout << "--- Motropoled GS! ---" << endl << endl;

  // Caso eccitato 2p
  cout << "--- Metropoling 2P... ---" << endl;
  std::vector<double> Positions_2p = Metropolis(pdf_2p, x_2p, y_2p, z_2p, M, del_2p, count_2p, rnd, coord_2p, true, N_coord); // salvo le coordinate
  cout << "--- Metropoled 2P! ---" << endl << endl;

  // stampo il file con le posizioni
  cout << "--- Printing coordinates... ---" << endl;
  ofstream out;
  out.open("dati_gauss/Coordinate_2p.txt");

  // Controllo la corretta apertura del file
  if (!out.is_open()) {
    cerr << "Error: unable to open Coordinate_2p.txt" << endl;
    exit(1);
  }

  for(int i=0; i < N_coord; i++) {
    out << coord_2p[i][0] << " " << coord_2p[i][1] << " " << coord_2p[i][2] << endl;
  }

  out.close();
  cout << "--- Saved coordinates! ---" << endl;
  
  // Stampo i risultati su file esterno, dividendo per blocchi
  Ave_Block(Positions_gs, N, "dati_gauss/Posizioni_gs_gauss.txt");
  Ave_Block(Positions_2p, N, "dati_gauss/Posizioni_2p_gauss.txt");

  cout << "--- Printed data ---" << endl;
  cout << "--- End of code ---" << endl;
  
  rnd.SaveSeed();
  return 0;
}

//========================================================================
// DEFINIZIONE FUNZIONI
//========================================================================

// Funzione Metropolis: restituisce un array con posizioni che seguono la distribuzione pdf in ingresso, a partire da 3 coordinate iniziali
// In ingresso passo anche il vettore con le tre coordinate da salvare, facente parte delle matrice di coordinate da salvare
std::vector<double> Metropolis(double (*pdf)(double, double, double), double &x, double &y, double &z, int Nsim, double del, int & count, Random & rnd, std::vector<std::vector<float>> & Coord, bool save_coord, int Ncoord) {

  vector<double> Position(Nsim);
  int k = 0; // indice dell'array di coordinate
  bool save = false; // questo serve alla funzione Move_Metropolis per salvare le coordinate in un vettore SOLO alla fine
  
  for (int i = 0; i < Nsim; i++) { // Per ogni simulazione
	  Position[i] = sqrt(x*x+y*y+z*z);

    if (save_coord) { // se voglio salvare le coordinate
      if (i > (Nsim-Ncoord)) { // appena l'indice arriva alle ultime Ncoord coordinate attivo il bool di Move_Metropolis
        k++; // incremento l'indice
        save = true;
      }
    }
      
    Move_Metropolis(pdf, x, y, z, del, count, rnd, Coord[k], save);
  } // fine ciclo / step del random walk

  cout << "Metropolis called with " << Nsim << " values" << endl << "acceptance ratio: " << fixed << setprecision(3) << double(count)/Nsim << endl;
  return Position;
}

//========================================================================
// funzione che esegue uno step di algoritmo di Metropolis, cambiando per riferimento i valori di posizione in ingresso, in modo da poter salvare l'ultima raggiunta dopo equilibrazione
// In ingresso passo anche il vettore con le tre coordinate da salvare, facente parte delle matrice di coordinate da salvare
void Move_Metropolis(double (*pdf)(double, double, double), double & x, double & y, double & z, double del, int & count, Random & rnd, std::vector<float> & Coord_k, bool save_coord) {
  
      double x_new = rnd.Gauss(x, del*del); // posizione al passo successivo, usando la probabilità di Transizione gaussiana
      double y_new = rnd.Gauss(y, del*del);
      double z_new = rnd.Gauss(z, del*del);
  
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
      // salvo le coordinate estratte
      if (save_coord) {
        Coord_k[0] = x_new;
        Coord_k[1] = y_new;
        Coord_k[2] = z_new;
      }
}
