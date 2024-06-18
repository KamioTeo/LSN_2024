// ESERCIZIO 3: Stima di pigreco simulando l'esperimento di Buffon
// ATTENZIONE: questa è la prima versione del codice e utilizza pigreco stesso per l'astrazione casuale uniforme dell'angolo di orientazione dell'ago
// La versione corretta è la seconda

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
  std::vector<double> angle(
      M); // Vettore che conterrà le inclinazioni casuali dell'ago

  for (int i = 0; i < M; i++) { // Per ogni simulazione
    y[i] = rnd.Rannyu(
        0,
        d); // L'altezza casuale è tra 0 e d e segue una distribuzione uniforme
    angle[i] = rnd.Rannyu(0, 2 * M_PI); // L'angolo casuale è tra 0 e 2PI
  }

  vector<double> ave(
      N); // lista che conterrà la media cumulativa delle simulazioni al variare
          // del numero N di blocchi utilizzato
  vector<double> av2(N); // analogo ma contiene la media cumulativa al quadrato
  vector<double> sum_prog(N); // analogo ma è la somma cumulativa dei numeri
                              // casuali, per calcolo della dev std
  vector<double> su2_prog(N); // somma cumulativa dei quadrati dei numeri
                              // casuali, per calcolo della dev std
  vector<double> err_prog(N); // dev std cumulativa

  for (int i = 0; i < N; i++) { // per ogni blocco calcolo la media sul blocco
    double hit = 0;
    for (int j = 0; j < L; j++) { // per ogni numero casuale nel blocco
      int k = j + i * L; // k-esimo numero casuale dallo 0-esimo nello 0-esimo
                         // blocco ecc
      double h = y[k] + length * sin(angle[k]); // calcolo la proiezione dell'altezza
                                         // di un estremo dell'ago sull'asse y
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
  out.open("dati/Dati_es1-3.txt");

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
