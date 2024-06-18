// ESERCIZIO 3 Versione 2: Stima di pigreco simulando l'esperimento di Buffon ma senza utilizzare pigreco stesso nel calcolo (nella versione precedente veniva utilizzato per l'estrazione uniforme dell'angolo)

// Genero la posizione casuale di un estremo dell'ago
// Genero altre due coordinate (punto del secondo estremo dell'ago) che vivono in un quadrato di lato pari alla lunghezza dell'ago
// Se questo punto è interno al cerchio tracciato dall'ago, allora lo tengo, altrimenti lo scarto
// La proiezione sull'asse verticale del secondo estremo è pari alla coordinata y del secondo punto estratto per la lunghezza dell'ago diviso la distanza del punto dall'origine

#include "../RandomNumberGenerator/random.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {

  // Generazione del random number generator
  Random rnd("../RandomNumberGenerator/Primes",
			 "../RandomNumberGenerator/seed.in");

  // Definisco le variabili che mi servono per l'analisi
  const int M = 1E7;   // Numero di estrazioni
  const int N = 100;   // Numero di blocchi
  const int L = M / N; // numero di lanci per blocco

  const int d = 1;               // Distanza tra le linee orizzontali
  const double length = 0.85 * d; // Lunghezza ago come percentuale della distanza tra le linee

  std::vector<double> y(M); // Vettore che conterrà le altezze casuali dell'ago
  std::vector<double> p_x(M); // vettore che conterrà la coordinata x del punto del secondo estremo dell'ago
  std::vector<double> p_y(M); // vettore che conterrà la coordinata y del punto del secondo estremo dell'ago
 
  for (int i=0; i<M; i++) { // Per ogni simulazione
	y[i] = rnd.Rannyu(0,d); // L'altezza casuale è tra 0 e d e segue una distribuzione uniforme
	p_x[i] = rnd.Rannyu(0,length); // coordinata x del secondo estremo
  p_y[i] = rnd.Rannyu(0,length); // coordinata y del secondo estremo
  }

  vector<double> ave(N); // lista che conterrà la media cumulativa delle simulazioni al variare del numero N di blocchi utilizzato
  vector<double> av2(N); // analogo ma contiene la media cumulativa al quadrato
  vector<double> sum_prog(N); // analogo ma è la somma cumulativa dei numeri casuali, per calcolo della dev std
  vector<double> su2_prog(N); // somma cumulativa dei quadrati dei numeri casuali, per calcolo della dev std
  vector<double> err_prog(N); // dev std cumulativa

  for (int i=0; i<N; i++) { // per ogni blocco calcolo la media sul blocco della stima di pigreco
	double hit=0;
	for (int j=0; j<L; j++) { // per ogni numero casuale nel blocco
	  int k = j + i * L; // k-esimo numero casuale dallo 0-esimo nello 0-esimo blocco ecc
    double r = sqrt(pow(p_x[k],2)+pow(p_y[k],2)); // distanza del secondo estremo dal centro dell'ago
    while(r>length) { // se il secondo estremo è esterno al cerchio formato dall'ago allora lo scarto
      p_x[k] = rnd.Rannyu(0,length);
      p_y[k] = rnd.Rannyu(0,length);
      r=sqrt(pow(p_x[k],2)+pow(p_y[k],2));
    } 
    
	  double h = y[k] + length * p_y[k] / r; // calcolo l'altezza dell'ago
	  if (h <= 0 or h >= d) { // controllo se l'ago interseca la linea sotto o quella sopra
		hit++;     // conto l'intersezione
	  }
	}

	ave[i] = (2 * length * L) / (d * hit); // valor medio di pigreco nel blocco i
	av2[i] = pow(ave[i], 2);          // (r_i)^2 quadrato della media
  }

  // Genero un file esterno in cui verranno salvate le medie cumulative e
  // relativo errore di pigreco
  ofstream out;
  out.open("dati/Dati_es1-3v2.txt");

  // Controllo la corretta apertura del file
  if (!out.is_open()) {
	cerr << "Error: unable to open Dati_es1-3.txt" << endl;
	exit(1);
  }

  // Prima riga del file esterno
  out << "Media" << " Errore" << endl;

  // ora calcolo la media cumulativa, sommando in progressione il contributo di ogni blocco
  for (int i = 0; i < N; i++) { // per ogni blocco
	for (int j = 0; j <= i; j++) {
	  sum_prog[i] += ave[j]; // SUM_{j=0,i} pi_j , sommo tutte le medie stimate fino al blocco i+1-esimo
	  su2_prog[i] += av2[j]; // SUM_{j=0,i} (pi_j)^2, sommo tutte le medie al quadrato stimate fino al blocco i+1-esimo
	}
	sum_prog[i] /= (i + 1); // calcolo la media cumulativa
	su2_prog[i] /= (i + 1); // media dei quadrati cumulativa
	err_prog[i] = error(sum_prog, su2_prog, i); // Deviazione standard calcolata fino al blocco i-esimo

	out << sum_prog[i] << " " << err_prog[i] << endl; // stampo i risultati
  }

  out.close();
  rnd.SaveSeed();

  return 0;
}
