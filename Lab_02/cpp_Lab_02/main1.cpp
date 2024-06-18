// ESERCIZIO 1: Integrazione Monte Carlo di M_PI/2*cos(M_PI*x/2) da 0 a 1, utilizzando un'estrazione uniforme nel primo punto e l'importance sampling nel secondo punto

#include <iostream>
#include <fstream>
#include <string>
#include "../RandomNumberGenerator/random.h"

using namespace std;

// Funzione di cui calcolare l'integrale, f(x)
double func(double x) {
    return (M_PI/2) * cos(M_PI*x/2);
}

// Funzione della pdf da usare nel secondo punto, p(x)
double pdf(double x) {
    return 2. * (1-x);
}

// Generatore distribuzione pdf del tipo 2(1-x), F^-1(y)
double cumulative_inverse(double y){
	  return 1-sqrt(1-y);
  }

// g(x)=f(x)/p(x), con x che segue la distribuzione p(x), è la variabile di cui calcolare la media per ottenere la stima dell'integrale in questo secondo metodo
double g(double x){
	return func(cumulative_inverse(x))/pdf(cumulative_inverse(x));
}

int main (int argc, char *argv[]){

// Generazione del random number generator 
	Random rnd("../RandomNumberGenerator/Primes", "../RandomNumberGenerator/seed.in");

// Definisco le variabili che mi servono per l'analisi
	const int M = 1E6; // Numero di estrazioni
	const int N = 100; // Numero di blocchi
	
	const int a = 0; // Primo estremo d'integrazione
	const int b = 1; // Secondo estremo d'integrazione

// Punto 1: estrazione uniforme della variabile casuale
	std::vector<double> Integral_unif(M); // Vettore che conterrà gli addendi per il calcolo della media dell'integrale

	for (int i = 0; i < M; i++) { // Per ogni simulazione
		Integral_unif[i] = (b-a)*func(rnd.Rannyu(a, b)); // aggiungo al vettore l'addendo che serve per il calcolo della media dell'integrale sul blocco. I numeri uniformi vanno estratti nell'intervallo d'integrazione
	}

// Uso la funzione Ave_Block per salvare su file esterno le medie cumulative e i loro errori, al variare degli N blocchi
// queste medie per ogni blocco sono i valori stimati dell'integrale, questo perché l'approssimazione dell'integrale sui dati di un blocco è I=(a-b)/L*sum(f(x_i)), e quindi facendo la media a blocchi di Integral_unif[i] divido già per L=numero di dati per blocco
	Ave_Block(Integral_unif, N, "dati/Dati_es2-1-1.txt");

// PUNTO 2: Analogo ma con estrazione da distribuzione non uniforme
// In questo caso i numeri sono estratti seguendo una distribuzione come la funzione pdf definita sopra [2(1-x)]. Per farlo si utilizza il metodo dell'inversa della cumulativa, quindi il set di variabili di cui lacolare la media è
// In questo caso non impongo i limiti d'integrazione perché, per come è ricavata l'inversa, vale solo tra 0 e 1

	std::vector<double> Integral_non_unif(M); // Vettore che conterrà gli addendi per il calcolo della media dell'integrale

	for (int i = 0; i < M; i++) { // Per ogni simulazione
		Integral_non_unif[i] = g(rnd.Rannyu()); // Creo l'addendo per il calcolo della media dell'integrale con il secondo metodo
	}

// Analogo a prima, uso la funzione Ave_Block ma sul nuovo vettore
	Ave_Block(Integral_non_unif, N, "dati/Dati_es2-1-2.txt");

   rnd.SaveSeed();
   return 0;
}
