// ESERCIZIO 1
/*
- Stampa su file esterno la media cumulativa (e relativo errore) fatta su 100 blocchi da 10^4 numeri casuali uniformi (per verificare che tende a 1/2)
- Stampa su file esterno la deviazione standard cumulativa (e relativo errore) di quei numeri (per verificare che tende a 1/12)
- Stampa su file esterno il chi quadro associato al dataset estratto
- Genera un file esterno con 3 distribuzioni di numeri casuali (uniforme, esponenziale e lorentziana)
*/

#include <iostream>
#include <fstream>
#include <string>
#include "../RandomNumberGenerator/random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Generazione del random number generator 
	Random rnd("../RandomNumberGenerator/Primes", "../RandomNumberGenerator/seed.in");

// Definisco le variabili che mi servono per l'analisi
	const int M = 1E6; // Numero di estrazioni
	const int N = 100; // Numero di blocchi

// Punto 1: test della media di numeri pseudocasuali con distribuzione uniforme

// Vettore che conterrà i numeri psudocasuali
	std::vector<double> Unif(M);

	for (int i = 0; i < M; i++) { // Per ogni simulazione
		Unif[i] = rnd.Rannyu(); // aggiungo al vettore un numero pseudocasuale con distribuzione uniforme
	}

// Uso la funzione Ave_Block per salvare su file esterno le medie cumulative e i loro errori, al variare degli N blocchi
	Ave_Block(Unif, N, "dati/Dati_es1-1-1.txt");
// ovvero la j-esima riga del file è il risultato della media sui primi j blocchi, in il valore numerico di ogni blocco è la media su L numeri casuali.
	
// PUNTO 2: Test della deviazione standard di numeri pseudocasuali con distribuzione uniforme
// Analogo a prima, salvo su file esterno le medie cumulative e relativo errore delle varianze dei numeri pseudocasuali generati in precedenza.
	
	std::vector<double> Unif2(M); // Vettore che conterrà le varianze dei numeri psudocasuali trovati in precedenza

	for (int i = 0; i < M; i++) { // Per ogni simulazione
		Unif2[i] = pow(Unif[i]-0.5, 2); // calcolo la varianza usando i numeri pseudocasuali trovati in precedenza
	}

// Analogo a prima, uso la funzione Ave_Block ma sul nuovo vettore
	Ave_Block(Unif2, N, "dati/Dati_es1-1-2.txt");

	
// PUNTO 3: Calcolo del chi quadro
// divido l'intervallo [0,1] in n_bins sottointervalli
	int n_bins = 100;
	int L = M/N; // numeri per blocco
	int exp = L/n_bins; // valore atteso di numeri casuali in ogni sottointervallo
	std::vector<double> chi(N); // conterrà i valori di chi quadro calcolati per blocco

// per ogni blocco calcolo quanti numeri sono effettivamente finiti in ogni intervallino
	for (int i = 0; i < N; i++) { // per ogni blocco
        std::vector<double> n_i(n_bins, 0.0); // vettore che contiene il numero totale di valori osservato in ogni intervallino, lo inizializzo a 0

        for (int j = i * L; j < (i * L) + L; j++) { // per ogni numero casuale nel blocco
            int bin_index = static_cast<int>(Unif[j] * n_bins); // calcolo l'indice del bin corrispondente al numero casuale.
// I numeri casuali sono in [0,1) quindi bin_index sarà un numero in [0,n_bins-1). Castando il risultato del prodotto ottengo un intero che uso per incrementare il relativo elemento nell'array n_i. Infatti, se ad esempio avessi Unif[0]=0.356 lo moltiplico per n_bins = 100 e diventa 35 a causa del cast (come la funzione floor approssima in difetto), quindi vuol dire che Unif[0] appartiene al 34_esimo intervallo, così come tutti gli eventuali numeri casuali che stanno in [0.35, 0.3599]
            n_i[bin_index]++; // incremento il conteggio del bin corrispondente
        }

        double chi_squared = 0.0; // definisco la variabile chi quadro del singolo blocco

        for (const auto& n : n_i) { // per ogni numero di numeri casuali di ogni intervallino
            chi_squared += std::pow((n - exp), 2) / exp; // incremento la sommatoria per il calcolo del chi quadro in questo blocco
        }

        chi[i]=chi_squared; // aggiungo il valore del chi quadro alla lista
    }

	Print_File(chi, "dati/Dati_es1-1-3.txt");
	
// ESERCIZIO 2:
// Genero un file esterno in cui salvo tre colonne di dati casuali, generati da distribuzione uniforme, esponenziale e di Lorentz
	ofstream out;
   out.open("dati/Dati_es1-2.txt");
   out << "Uniform" << " Exponential" << " Lorentz" << endl;

// Controllo sulla corretta apertura del file
	if(!out.is_open()){
		cerr << "Error: unable to open Dati_es1-2.out.txt" << endl;
		exit(1);
	}

   for(int i=0; i<1E6; i++){
	out << rnd.Rannyu() << " " << rnd.Exponential(1) << " " << rnd.Lorentz(1,0) << endl;
   }
	
   out.close();

   rnd.SaveSeed();
   return 0;
}
