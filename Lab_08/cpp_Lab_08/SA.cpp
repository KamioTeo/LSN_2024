// Implementazione dell'algoritmo SA ricercando mu e sigma nell'intervallo
// trovato con il codice del file delta_finder.cpp

// Variational Monte Carlo, Metropolis sampling
// uso unità di hbar=1 e m=1

#include "SA.h"
using namespace std; 

int main(int argc, char *argv[]) {
	
	Input();

	if (SA == 1) { // se voglio utilizzare l'algoritmo SA
		ofstream out;
		out.open("./dati/SA_results.dat");
		out << "Temp" << setw(wd) << "mu" << setw(wd) << "sigma" << setw(wd) << "energy"<< setw(wd) << "err" << endl;
		
		double accepted = 0; // Conteggio delle volte in cui accetto la mossa proposta nel Simulated Annealing
		
		const double tolerance = 1e-9; // Soglia di tolleranza per evitare che le approssimazioni taglino l'ultimo valore del ciclo for

		for (double temp = TMax; temp >= TMin - tolerance; temp -= TStep) {
	
			// Energia nelle configurazioni attuale e proposta
			vector<double> Energy(M);
			vector<double> Energy_new(M);
			
			double x_start = 0.; // valore iniziale campionamento |psi|^2
			double x_start_new = 0.;
			
			// estraggo i valori proposti di mu e sigma
			double mu_new = rnd.Rannyu(mu-deltaSA/2., mu+deltaSA/2.);
			double sigma_new = rnd.Rannyu(sigma-deltaSA/2., sigma+deltaSA/2.);
	
			// Equilibrio il sistema ed estraggo delta ottimale per l'energia vecchia e quella nuova proposta
			double delta = Metropolis_delta(mu, sigma, neq, x_start);
			double delta_new = Metropolis_delta(mu_new, sigma_new, neq, x_start_new);
			
			//cout << "------OLD-----" << endl;
			// Campiono la configurazione attuale
			vector<double> x = Metropolis(mu, sigma, delta, M, x_start);
			
			//cout << "------NEW-----" << endl;
			// Campiono la configurazione nuova
			vector<double> x_new = Metropolis(mu_new, sigma_new, delta_new, M, x_start_new);
			
			// Calcolo i valori di aspettazione dell'Hamiltoniana nelle due configurazioni
			for (int i = 0; i < M; i++) {
				Energy[i] = eigen_val(x[i], mu, sigma);
				Energy_new[i] = eigen_val(x_new[i], mu_new, sigma_new);	
			}
			
			// calcolo le due energie a blocchi, in base a quale scelgo calcolo anche l'errore
			double E = Ave_Block(Energy, N);
			double E_new = Ave_Block(Energy_new, N);
			double E2 = 0.;
			double E_sigma = 0.;
		
			// Eseguo Metropolis con probabilità di Boltzmann
			double A = exp(-(E_new-E)/temp);

			if( A>1 ) { A = 1.; } // ovvero se la nuova energia è minore di quella vecchia allora l'esponente è positivo e quindi l'esponenziale è maggiore di 1, ovvero accetto sempre
			double r = rnd.Rannyu();
			// dato che r è minore di 1, in realtà non serve contrllare che A sia più grande di 1, basta vedere se è più grande di r
			if( r <= A ) { // se è vero accetto la nuova posizione, quindi mu e sigma diventano mu_new e sigma_new
				mu = mu_new;
				sigma = sigma_new;
				// sovrascrivo l'energia nuova così poi la stampo
				E = E_new;
				// Calcolo l'errore dell'energia della nuova configurazione
				E2 = Ave2_Block(Energy_new, N);
				E_sigma = sqrt( (E2-pow(E_new, 2))/N );
				
				accepted += 1.; // Conto il fatto che ho accettato
			}
			else { // Non accettato, il valore resta quello vecchio
				// calcolo l'errore della vecchia configurazione
				E2 = Ave2_Block(Energy, N);
				E_sigma = sqrt( (E2-pow(E, 2))/N );
			}
			
			out << temp << setw(wd) << mu << setw(wd) << sigma << setw(wd) << E << setw(wd) << E_sigma << endl;
		}

		cout << "Acc/att: " <<  accepted/Nsteps << endl;
	}
	else { // altrimenti non faccio l'algoritmo SA ma calcolo l'energia a blocchi relativa alla funzione d'onda di test con i parametri mu e sigma estratti
		// Energia
		vector<double> Energy(M);
		// valore iniziale campionamento |psi|^2
		double x_start = 0.;
		// Equilibrio il sistema ed estraggo delta ottimale
		double delta = Metropolis_delta(mu, sigma, neq, x_start);
		// Campiono la configurazione attuale
		vector<double> x = Metropolis(mu, sigma, delta, M, x_start);
		// Calcolo i valori di aspettazione dell'Hamiltoniana
		for (int i = 0; i < M; i++) {
			Energy[i] = eigen_val(x[i], mu, sigma);
		}
		// Stampo l'energia a blocchi
		Ave_Block(Energy, N, "./dati/Min_Energy.dat");
		// Stampo il campionamento di |psi|^2
		Print_File(x, "./dati/Psi2_sampling.dat");
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

// calcolo di (H*psi)/psi
double eigen_val(double x, double mu, double sigma) {
	// calcolo derivata seconda della funzione d'onda di prova
	double exp1 = exp(-pow((x-mu)/sigma, 2)/2.);
	double exp2 = exp(-pow((x+mu)/sigma, 2)/2.);
	double fac1 = pow((x-mu)/sigma, 2);
	double fac2 = pow((x+mu)/sigma, 2);
	return 0.5*pow(sigma, -2)*( 1-(fac1*exp1+fac2*exp2) / psi_trial(x, mu, sigma) ) + potential(x);
}

//================================================
// funzione che restituisce il valore ottimale di delta (che garantisce un'accettazione del 50%) dati mu, sigma e gli step di equilibrazione
double Metropolis_delta(double mu, double sigma, int neq, double& x) {
	double iterations = 0; // conteggio iterazioni per la ricerca di delta
	double delta = 4.26;  // valore iniziale proposto di delta, è la media dei delta ottenuti negli intervalli ristretti di mu e sigma
	double acceptance = 0.; // accepted/attempted
	
	while(acceptance > 0.51 || acceptance <0.49) {

		if (acceptance > 0.51) { delta += 0.1; }
		else if (acceptance < 0.49) { delta -= 0.1; }
		
		double accepted = 0; // Conteggio delle volte in cui accetto la mossa proposta

		for (int i = 0; i < neq; i++) {
			double x_new = rnd.Rannyu(x-delta/2, x+delta/2);
			
			// dato che T è simmetrica, il rapporto q è dato solo dal rapporto tra la probabilità da campionare valutata in x_new e x_old
			double A = pow(psi_trial(x_new, mu, sigma)/psi_trial(x, mu, sigma), 2);
			
			if( A>1 ) { A = 1.; }
			double r = rnd.Rannyu();
			if( r <= A ) { // se è vero accetto la nuova posizione, quindi x diveta x_new
				x = x_new;
				accepted += 1.; // Conto il fatto che ho accettato
			}
			// Se non accetto, il valore resta quello vecchio, quindi sovrascrivo x_new con un altro tentativo a partire sempre da x
		}
		
		acceptance = accepted/neq;
		iterations += 1;
	}
		
	return delta;
}

//===================================================================
// funzione che restituisce un vettore di nstep valori di x che seguono la distribuzione di probabilità |psi|^2, prendendo in ingresso mu, sigma, delta, il valore di step e il valore iniziale di x ottenuto da dopo l'equilibrazione
vector<double> Metropolis(double mu, double sigma, double delta, int nstep, double& x) {
	double accepted = 0; // Conteggio delle volte in cui accetto la mossa proposta
	vector<double> x_metro(nstep);  // vettore che contiene tutti i valori di x che seguono la distribuzione di |psi|^2
	
	for (int i = 0; i < nstep; i++) {
		double x_new = rnd.Rannyu(x-delta/2, x+delta/2); // genero il nuovo punto estratto
		
		// dato che T è simmetrica, il rapporto q è dato solo dal rapporto tra la probabilità da campionare valutata in x_new e x_old
		double A = pow(psi_trial(x_new, mu, sigma)/psi_trial(x, mu, sigma), 2);
		
		if( A>1 ) { A=1.; }
		double r = rnd.Rannyu();
		if( r <= A ) { // se è vero accetto la nuova posizione, quindi x diveta x_new
			x = x_new;
			accepted += 1.; // Conto il fatto che ho accettato
		}
		// Se non accetto, il valore resta quello vecchio
		x_metro[i] = x;
	}

	//cout << "Acc/att: " <<  accepted/nstep << " delta: " << delta << " mu: " << mu << " sigma: " << sigma << endl;
	
	return x_metro;
}

//===========================================================
void Input(void) {
	ifstream ReadInput;
	
	cout << "_______________VMC_______________" << endl;
	cout << "Simulated Annealing code" << endl << endl;
	
	//Read input informations
	ReadInput.open("SA_input.in");
	
	ReadInput >> neq;
	ReadInput >> M;
	ReadInput >> N;
	
	cout << "equilibration steps = " << neq << endl;
	cout << "Number of steps after equilibration = " << M << endl;
	cout << "Number of blocks = " << N << endl;

	ReadInput >> TMax;
	ReadInput >> TMin;
	ReadInput >> TStep;
	
	cout << "Temperature range = [" << TMin << ";" << TMax << "]" << " step = " << TStep << endl;

	Nsteps = static_cast<int>((TMax - TMin) / TStep) + 1; // numero totale di step eseguiti qui di seguito
	cout << "Temperature steps: " << Nsteps << endl;

	ReadInput >> deltaSA;

	cout << "Delta for Simulated Annealing = " << deltaSA << endl;
	
	ReadInput >> mu;
	ReadInput >> sigma;

	cout << "Starting mu = " << mu << endl;
	cout << "Starting sigma = "<< sigma <<  endl << endl;

	ReadInput >> SA;

	if (SA == 0) {
		cout << "Sampling the wave function with the given parameters and printig the energy results" << endl;
	}
	
	ReadInput.close();
	
} // Fine Input()
