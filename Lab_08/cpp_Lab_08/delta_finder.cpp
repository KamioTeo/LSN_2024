// Codice per l'analisi del delta ottimale da utilizzare e per la
// ricerca dell'intervallo di mu e sigma in cui eseguire l'algoritmo SA

// Variational Monte Carlo, Metropolis sampling
// uso unità di hbar=1 e m=1

#include "delta_finder.h"
using namespace std;

int main(int argc, char *argv[]) {

	Input();

	vector<double> Energy(M); // valori di aspettazione di H (H*psi/psi) DOPO l'equilibrazione
	
	ofstream out;
	out.open("./dati/best_delta.dat");
	out << "mu" << setw(wd) << "sigma" << setw(wd) << "energy"<< setw(wd) << "err" << setw(wd) << "delta" << setw(wd) << "acc/att" << setw(wd) << "iterations" << endl;

	const double tolerance = 1e-9; // Soglia di tolleranza per evitare che le approssimazioni taglino l'ultimo valore del ciclo for
	for (double mu = muMin; mu <= muMax + tolerance; mu += muStep) { // per ogni valore di mu e sigma da testare
	    for (double sigma = sigmaMin; sigma <= sigmaMax + tolerance; sigma += sigmaStep) {
			double x_start = 0.; // valore iniziale dell'algoritmo

			// estraggo il valore di delta ottimale grazie a questa funzione, modifica per riferimento anche la posizione iniziale, il valore di accettazione e il numero di iterazioni
			double delta = Metropolis_delta(mu, sigma, neq, x_start, acceptance, iterations);
			// salvo il vettore x che contiene la distribuzione di psi^2 campionata con Metropolis (M valori), utilizzando la specifica coppia di valori mu e sigma
			vector<double> x = Metropolis(mu, sigma, delta, M, x_start);
			// calcolo la distribuzione dell'energia, è Hpsi/psi mediata sulla distribuzione |psi|^2
			for (int i = 0; i < M; i++) { // Per ogni step campionato dopo equilibrazione
				Energy[i] = eigen_val(x[i], mu, sigma); // calcolo il valore dell'energia
			}
			// calcolo media e dev std dell'energia
			double e_mean = Ave_Block(Energy, N);
			double e_mean2 = Ave2_Block(Energy, N);
			double e_sigma = sqrt( (e_mean2-pow(e_mean, 2))/N );
			
			out << mu << setw(wd) << sigma << setw(wd) << e_mean << setw(wd) << e_sigma << setw(wd) << delta << setw(wd) << acceptance << setw(wd) << iterations << endl;
		}
	}
	rnd.SaveSeed();
	return 0;
}


//=================================================================
// Funzione d'onda di prova
double psi_trial(double x, double mu, double sigma) {
	double exp1 = exp(-pow((x-mu)/sigma, 2)/2.);
	double exp2 = exp(-pow((x+mu)/sigma, 2)/2.);
	return exp1 + exp2;
}

// definizione del potenziale confinante
double potential(double x) {
	return pow(x, 4) - (5./2.)*pow(x, 2);
}

// calcolo di (Hpsi)/psi
double eigen_val(double x, double mu, double sigma) {
	// calcolo la derivata seconda della funzione d'onda di prova (serve per il calcolo dell'Hamiltoniana applicata a psi)
	double exp1 = exp(-pow((x-mu)/sigma, 2)/2.);
	double exp2 = exp(-pow((x+mu)/sigma, 2)/2.);
	double fac1 = pow((x-mu)/sigma, 2);
	double fac2 = pow((x+mu)/sigma, 2);
	return 0.5*pow(sigma, -2)*( 1-(fac1*exp1+fac2*exp2) / psi_trial(x, mu, sigma) ) + potential(x);
}

//===================================================================
// funzione che restituisce un vettore di nstep valori di x che seguono la distribuzione di probabilità |psi|^2, prendendo in ingresso mu, sigma, delta, il valore di step e il valore iniziale di x ottenuto da dopo l'equilibrazione
vector<double> Metropolis(double mu, double sigma, double delta, int nstep, double& x) {
	double accepted = 0; // Conteggio delle volte in cui accetto la mossa proposta
	// double x = 0.; // valore iniziale, non lo definisco qui perché lo prendo dalla funzione sopra che equilibra il sistema, in questo modo posso scegliere la posizione di partenza e non è obbligatoriamente zero
	vector<double> x_metro(nstep);  // vettore che contiene tutti i valori di x che seguono la distribuzione di |psi|^2
	
	for (int i = 0; i < nstep; i++) {
		double x_new = rnd.Rannyu(x-delta/2, x+delta/2); // estraggo il nuovo valore
		
		// dato che T è simmetrica, il rapporto q è dato solo dal rapporto tra la probabilità da campionare valutata in x_new e x_old
		double A = pow(psi_trial(x_new, mu, sigma)/psi_trial(x, mu, sigma), 2);
		
		if( A > 1 ) { A=1; }
		
		double r = rnd.Rannyu();
		if( r <= A ) { // se è vero accetto la nuova posizione, quindi x diveta x_new
			x = x_new;
			accepted += 1.; // Conto il fatto che ho accettato
		}
		// se non accetto, il valore resta quello vecchio
		x_metro[i] = x;
	}

	cout << "Acc/att: " <<  accepted/nstep << " delta: " << delta << " mu: " << mu << " sigma: " << sigma << endl;
	
	return x_metro;
}

//================================================
// funzione che restituisce il valore ottimale di delta (che garantisce un'accettazione del 50%) dati mu, sigma e gli step di equilibrazione. Modifica per riferimento la variabile che conta il numero di tentativi fatti "iterations" e l'accettazione "acceptance", in modo da salvarli poi su file esterno insieme ai valori di accettazione analizzati a blocchi (serve ad analizzarne la velocità di convergenza in termini di tentativi). Inoltre, modifica per riferimento il valore iniziale dell'algoritmo di Metropolis, in modo da utilizzarlo dopo l'equilibrazione
double Metropolis_delta(double mu, double sigma, int neq, double& x, double& acceptance, int& iterations) {
	iterations = 0; // conteggio iterazioni per la ricerca di delta
	double delta = 5;  // valore iniziale proposto di delta
	acceptance = 0.; // accepted/attempted
	
	while(acceptance > 0.51 || acceptance <0.49) {

		if (acceptance > 0.49) { delta += 0.1; }
		else if (acceptance < 0.51) { delta -= 0.1; }
		// per ogni delta testato azzero i conteggi
		double accepted = 0; // Conteggio delle volte in cui accetto la mossa proposta

		for (int i = 0; i < neq; i++) {
			double x_new = rnd.Rannyu(x-delta/2, x+delta/2);
			
			// dato che T è simmetrica, calcolo solo il rapporto tra la probabilità da campionare valutata in x_new e x_old
			double A = pow(psi_trial(x_new, mu, sigma)/psi_trial(x, mu, sigma), 2);
			
			if( A > 1 ) { A = 1; }
			double r = rnd.Rannyu();
			if( r <= A ) { // se è vero accetto la nuova posizione, quindi x diveta x_new
				x = x_new;
				accepted += 1.; // Conto il fatto che ho accettato
			} 
			// Se non accettato il valore resta quello vecchio, quindi sovrascrivo x_new con un altro tentativo a partire sempre da x
		}
		
		acceptance = accepted/neq; // calcolo l'accettazione
		iterations += 1; // conto il fatto che ho fatto un tentativo
	}
		
	return delta;
}


//===========================================================

void Input(void) {
	ifstream ReadInput;
	
	cout << "_______________VMC_______________" << endl;
	cout << "Code that finds the best delta value for Metropolis algorithm" << endl << endl;
	
	//Read input informations
	ReadInput.open("delta_finder_input.in");
	
	ReadInput >> neq;
	ReadInput >> M;
	ReadInput >> N;
	
	cout << "equilibration steps = " << neq << endl;
	cout << "Number of steps after equilibration = " << M << endl;
	cout << "Number of blocks = " << N << endl;
	
	ReadInput >> muMin;
	ReadInput >> muMax;
	ReadInput >> muStep;

	ReadInput >> sigmaMin;
	ReadInput >> sigmaMax;
	ReadInput >> sigmaStep;
	
	cout << "mu range = [" << muMin << ";" << muMax << "]" << " step = " << muStep << endl;
	cout << "sigma range = [" << sigmaMin << ";" << sigmaMax << "]" << " step = " << sigmaStep <<  endl << endl;

	ReadInput.close();
	
} // Fine Input()
