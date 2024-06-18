// ESERCIZIO 1.1 - Plain Vanilla Option pricing
// L'obiettivo è stimare tramite simulazione Monte Carlo il prezzo di un'opzione Europea tipo call e put al tempo presente $t=0$, a partire dalla previsione del prezzo dell'asset sottostante al tempo t=T (tempo di scadenza del contratto) S(T), utilizzando due tecniche di campionamento diverso

// Campionamento diretto
#include "../RandomNumberGenerator/random.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {

  // Generazione del random number generator
  Random rnd("../RandomNumberGenerator/Primes",
             "../RandomNumberGenerator/seed.in");

  // Definizioni variabili per il Numerical Option Pricing
  const double S0 = 100;  // Prezzo dell'asset a t=0
  const double T = 1;     // Tempo scadenza contratto (delivery time)
  const double K = 100;   // Prezzo accordato (strike price)
  const double r = 0.1;   // Tasso d'interesse, costante e senza rischio
  const double sigma = 0.25; // volatilità

  const int M = 100000; // Numero di estrazioni

  // Calcolo il prezzo della Call e del Put per ogni simulazione, utilizzando un campionamento diretto di S(T) e riportando il valore all'istante t=0 tramite moltiplicazione per l'esponenziale del tasso d'interesse

  double z = 0;      // Variabile normale per il campionamento diretto
  double S = 0;      // Conterrà il prezzo dell'asset al tempo T
  double Profit = 0; // Conterrà il profitto al tempo T

  std::vector<double> Call(M); // Vettore che conterrà tutti i prezzi delle Call riportati al tempo t=0
  std::vector<double> Put(M); // Vettore che conterrà tutti i prezzi Put riportati al tempo t=0

  for (int i = 0; i < M; i++) { // Per ogni simulazione
    z = rnd.Gauss(0., T); // z è distribuito normalmente tra 0 e T, serve al GBM
    S = S0 * exp((r - 0.5*sigma*sigma)*T - sigma*z*sqrt(T)); // prezzo dell'asset a t=T come GBM (Geometric Brownian Motion)

    // Ora calcolo il profitto come max(0, S(T) - K), ovvero lo considero solo se S(T) > K
    Profit = S - K;
    if (Profit > 0) {
      Call[i] = Profit * exp(-r * T); // lo riporto al tempo t=0 tramite l'esponenziale
    } else {
      Put[i] = -Profit * exp(-r * T); // il profitto delle opzioni put è esattamente il contrario
    }
  }

  const int N = 100; // Numero di blocchi

  Ave_Block(Call, N, "prove/Direct_Call.txt");
  Ave_Block(Put, N, "prove/Direct_Put.txt");

  // Campionamento discreto
  // Divido l'intervallo [0,T] in 10^2 intervalli

  int time_divisions = 100;
  std::vector<double> time(time_divisions + 1); // vettore che conterrà t_0=0, t_1, ... , t_100

  for (int i = 0; i <= time_divisions; i++) {
    time[i] = i * T / time_divisions;
  }

  std::vector<double> Call2(M); // Vettore che conterrà tutti i prezzi delle Call riportati al tempo t=0
  std::vector<double> Put2(
      M); // Vettore che conterrà tutti i prezzi Put riportati al tempo t=0

  for (int i = 0; i < M; i++) { // Per ogni simulazione
    S = S0;                     // definisco S al tempo iniziale
    for (int j = 1; j <= time_divisions; j++) {
      z = rnd.Gauss(0., 1); // z è distribuito normalmente
      S *= exp(
          (r - 0.5 * sigma * sigma) * (time[j] - time[j - 1]) -
          sigma * z * sqrt(time[j] - time[j - 1])); // Formula per il prezzo dell'asset come percorso di una GBM, ogni volta incremento il tempo fino ad arrivare a t=T
    }

    // Ora calcolo il profitto come max(0, S(T) - K), ovvero lo considero solo se positivo
    Profit = S - K;
    if (Profit > 0) {
      Call2[i] = Profit * exp(-r * T); // lo riporto al tempo t=0 tramite l'esponenziale
    } else {
      Put2[i] = -Profit * exp(-r * T); // il profitto delle opzioni put è esattamente il contrario
    }
  }

  Ave_Block(Call2, N, "prove/Discrete_Call.txt");
  Ave_Block(Put2, N, "prove/Discrete_Put.txt");

  rnd.SaveSeed();
  return 0;
}
