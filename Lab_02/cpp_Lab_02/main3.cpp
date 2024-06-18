// ESERCIZIO 2-2: Random Walks continui
// Calcolare la distanza media (come radice della media dei quadrati delle posizioni finali) di una serie di random walks che procedono nello spazio in direzioni continue, al variare del numero di passi effettuati

#include <iostream>
#include <fstream>
#include <string>
#include "../RandomNumberGenerator/random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Generazione del random number generator 
	Random rnd("../RandomNumberGenerator/Primes", "../RandomNumberGenerator/seed.in");

// Definisco le variabili che mi servono per l'analisi
	const int Ndati = 1E6; 		// Numero di estrazioni
	const int P = 100; 			// Numero di passi totali
	const int M = Ndati / P;	// Numero di random walks lunghi P passi
	const int N = 100; 			// Numero di blocchi
	const int L = M / N; 		// Dati per blocco

	const double a = 1.0; 		// Lunghezza di un passo

	std::vector<double> phi(Ndati); // Conterrà gli angoli phi casuali per la direzione spaziale continua
	std::vector<double> theta(Ndati); // Conterrà gli angoli theta casuali per la direzione spaziale continua

	for(int i=0; i< Ndati; i++){
		phi[i] = rnd.Rannyu(0, 2*M_PI); // phi è uniforme e va da 0 a 2PI
		theta[i] = rnd.ThetaVersor(); // theta non è uniforme perché deve bilanciare il contributo dello Jacobiano in coordinate sferiche
	}

/* Logica del codice:
 Divido le Ndati=10e6 direzioni casuali in M=10e4 gruppi da P=100 passi.
 Metto tutto in una matrice MxP, in modo che ogni colonna sia un passo fissato (da 1 a P) e ogni riga sia una simulazione di Random Walk (da 1 a M).
 Fissato il passo (e quindi la colonna) eseguo una media a blocchi (N=100 blocchi) delle posizioni quadro raggiunte da tutti i Random Walk.
 Il risultao per ogni passo è dato dalla radice di questo valore, e l'errore associato si ottiene tramite propagazione degli errori.
*/
	
	std::vector<std::vector<double>> ave_bloc(N, std::vector<double>(P)); // lista che conterrà la media su un blocco della posizione quadro raggiunta dopo un certo numero di passi
	std::vector<std::vector<double>> ave2_bloc(N, std::vector<double>(P)); // analogo ma ogni termine è al quadrato, serve per l'errore

// Vector per la divisione a blocchi
	std::vector<double> av(P); // Per ogni passo, contiene la media sui blocchi della posizione al quadrato raggiunta
	std::vector<double> av2(P); // analogo ma al quadrato, serve per l'errore

	std::vector<double> err_av(P); // Conterrà gli errori per ogni passo della media della posizione quadro nei blocchi
	
	std::vector<std::vector<double>> s(M, std::vector<double>(P)); // matrice che contiene M random Walk di P passi, ovvero a passo fissato contiene la distanza quadro raggiunta, ad esempio s[k] è una lista di P distanze quadro riguardanti il k-esimo Random Walk

// prima calcolo le posizioni quadro finali raggiunte da ogni random walk a passo fissato
	for (int m = 0; m < M; m++) { // Per ogni random walk
	// Ogni simulazione parte dall'origine
	    double dir_x = 0.0; 
	    double dir_y = 0.0;
	    double dir_z = 0.0;
	    for (int p = 0; p < P; p++) { // ad ogni passo calcolo la distanza quadro raggiunta
	        int k = p + m * P;
		// la posizione al passo successivo si basa su quella raggiunta al passo precedente
	        dir_x += a * sin(theta[k]) * cos(phi[k]);
	        dir_y += a * sin(theta[k]) * sin(phi[k]);
	        dir_z += a * cos(theta[k]);
	        s[m][p] = dir_x * dir_x + dir_y * dir_y + dir_z * dir_z; 
	    }
	}

// Ora per ogni passo faccio la media a blocchi dei risultati ottenuti
	for (int p = 0; p < P; p++) {
	    for (int i = 0; i < N; i++) {
	        double sum = 0.0;
	        for (int j = 0; j < L; j++) {
	            int k = j + i * L;
	            sum += s[k][p];
	        }
	        ave_bloc[i][p] = sum / L; // media sui blocchi
	        ave2_bloc[i][p] = ave_bloc[i][p] * ave_bloc[i][p];
	    }
	}

  // Genero un file esterno in cui verranno salvate le medie cumulative e relativo errore della distanza raggiunta
	  ofstream out;
	  out.open("dati/Dati_es2-2-2.txt");

  // Controllo la corretta apertura del file
	  if (!out.is_open()) {
	    cerr << "Error: unable to open Dati_es2-2-2.txt" << endl;
	    exit(1);
	  }
	
	// Prima riga del file esterno
	out << "Media" << " Errore" << endl;

	// ora calcolo per ogni passo la media delle medie sui blocchi (delle distanze al quadrato)
	// ovvero è come fare la media cumulativa e prendere il risultato ottenuto cumulando tutti i blocchi

	for (int p = 0; p < P; p++) {
	    double sum = 0.0;
	    double sum2 = 0.0;
	    for (int n = 0; n < N; n++) { // per ogni blocco
	        sum += ave_bloc[n][p];
	        sum2 += ave2_bloc[n][p];
	    }
	    av[p] = sum / N; // media sui blocchi delle distanze quadro
	    av2[p] = sum2 / N;
	    err_av[p] = error(av, av2, p); // errore sulla media delle distanze quadro
		
		out << sqrt(av[p]) << " " << err_av[p] / (2. * sqrt(av[p])) << endl; // Stampo su file esterno la distanza, considero la propagazione degli errori
	}
	
   out.close();
   rnd.SaveSeed();
   return 0;
}
