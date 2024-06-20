/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Implementazione del modello di Ising 1D: Montecarlo e Gibbs

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
  Input(); //Inizialization, setta il Random Number generator e i parametri di Ising, carico il metodo da utilizzare (Metropolis o Gibbs) e carico se calcolare tutti i parametri fisici ad una specifica temperatura oppure se calcolarli in un range

  if (Equilibration == 1) { // il codice stampa i valori dell'energia per spin durante l'equilibrazione
	int equil_steps = 5000;
	vector<double> Energy(equil_steps); // vettore che conterrà l'energia calcolata ad ogni step
	Reset(1);
	// porto il sistema all'equilibrio e salvo l'energia al variare dei passi
	for(int istep=1; istep <= equil_steps; istep++) { // per ogni incremento temporale (step)
        Move(metro); // aggiorno gli spin
        Measure(); // Misuro tutti i parametri nella nuova configurazione
        Energy[istep-1] = walker[iu]/nspin; // Salvo l'energia della nuova configurazione
    }
	// salvo i valori su file esterno
	if (metro == 1)
		Print_File(Energy, "dati/Metropolis/Equilibration.dat");
	else
		Print_File(Energy, "dati/Gibbs/Equilibration.dat");
	return 0;
  }
	
  if (Tsequence == 1) { // se voglio calcolare le grandezze nel range di temperatura
    for(int i = 0; i<=15; ++i) {
      temp = 0.5 +0.1*i; // sovrascrivo la temperatura caricata, vado da 0.5 a 2 a passi di 0.1
      beta = 1./temp;

      // for( int i=0; i< 100000; ++i)   // Equilibration
        // Move(metro);
    
      for(int iblk=1; iblk <= nblk; ++iblk) { // Simulation
        Reset(iblk);   // Reset block averages
        for(int istep=1; istep <= nstep; ++istep) {
          Move(metro);
          Measure();
          Accumulate(); // Update block averages
        }
	  LastAverages(iblk); // stampo l'ultimo valore delle medie a blocchi di ogni grandezza fisica
	  }
      // ConfFinal(); // Write final configuration

	} // fine ciclo temperatura
  } else if (Tsequence == 0) { // stampo le grandezze termodinamiche alla temperatura di input
	// for( int i=0; i< 100000; ++i)   // Equilibration, non serve, converge subito
 //        Move(metro);
	  
  	for(int iblk=1; iblk <= nblk; ++iblk) {  // Simulation
	  Reset(iblk);   // per ogni blocco resetto i valori ed eseguo nstep
	  for(int istep=1; istep <= nstep; ++istep) { // per ogni step evolvo il sistema e salvo le grandezze termodinamiche
	    Move(metro); // aggiorno il sistema con metodo di Metropolis o Gibbs
	    Measure(); // Misuro le proprietà termodinamiche
	    Accumulate(); //Update block averages
	  }
	  Averages(iblk);   // Una volta finiti gli step, stampo le medie dei parametri di quel blocco
    }
	ConfFinal(); //Write final configuration dopo tutti i blocchi e salva il seed
  } else { // calcolo l'errore sull'energia al variare del numero di blocchi
	  int nsim = 1000000; // numero totale di simulazioni
	  for(int i = 1; i<=10; ++i) { // vario il numero di blocchi
        nblk = 10*i; // sovrascrivo il numero di blocchi da utilizzare
        nstep = int(nsim/nblk); // numero di simulazioni (step) per blocco
        // for( int i=0; i< 100000; ++i)   // Equilibration
        // Move(metro);
        
        for(int iblk=1; iblk <= nblk; ++iblk) { // Per ogni simulazione nel blocco
          Reset(iblk);   // Reset block averages
          for(int istep=1; istep <= nstep; ++istep) {
            Move(metro);
            Measure();
            Accumulate(); // Update block averages
          }
	      EnergyAverage(iblk); // stampo l'ultimo valore delle medie a blocchi dell'energia
	    }
	  }
    } // fine else
  return 0;
}

// Inizializzazione, setta il Random Number generator e carica i parametri in input
void Input(void) {
  ifstream ReadInput, Readconf, input;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs
  if(metro==1)
	  cout << "The program perform Metropolis moves" << endl;
  else
	  cout << "The program perform Gibbs moves" << endl;
	
  ReadInput >> nblk;
  cout << "Number of blocks = " << nblk << endl;
	
  ReadInput >> nstep;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> restart; // vale 0 se no, 1 se sì (ovvero se parto da una configurazione esistente, caricata da file esterno)
  if(restart==1)
	  cout << "The program restart from the last moves" << endl;
  else
	  cout << "The program restart again" << endl;
	
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
	
   if(restart)
	   input.open("seed.out"); // se devo partire da una configurazione esistente devo anche continuare le simulazioni usando lo stesso seed
   else
	  input.open("seed.in"); // uso il seed nuovo appena generato
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
	
  ReadInput >> Tsequence; // se 1 allora calcola i dati termodinamici al variare della temperatura

  if (Tsequence == 1)
	  cout << "Calculating values with different temperatures" << endl << endl;
  else
	  cout << "Calculating values using the temperature given" << endl << endl;

  ReadInput >> Equilibration; // se 1 allora salva un file con il valore dell'energia del sistema ad ogni passo eseguito
	
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
	
  iu2 = 4;
  im2 = 5;
  n_props = 6; //Number of observables

// Inizializzazione del reticolo, può essere casuale o dipendere da una specifica configurazione precedentemente salvata

//initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart) {
    cout<<"Configurazione iniziale caricata da file esterno"<<endl;
	Readconf.open("config.final");
      for (int i=0; i<nspin; ++i)
        Readconf >> s[i];   
  } else {
    cout<<"Configurazione iniziale casuale"<<endl;
    for (int i=0; i<nspin; ++i) {
      if(rnd.Rannyu() >= 0.5)
		  s[i] = 1;
      else
		  s[i] = -1;
    }
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
} // Fine Input()


void Move(int metro) { // funzione che aggiorna il reticolo in base all'algoritmo utilizzato, Metropolis o Gibbs
  int o; // coordinata spin estratto casualmente
  double A; // probabilità di accettazione
  // double energy_up, energy_down;

  for(int i=0; i<nspin; ++i) { // Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) { //Metropolis
    // INCLUDE YOUR CODE HERE
    // Calcolo la probabilità di accettazione, a partire dalla differenza in energia tra la configurazione attuale e quella nuova in cui flippo lo spin (nota segno meno in ingresso): uso la funzione Boltzmann
      A = min(1., exp(-beta * (Boltzmann(-s[o],o) - Boltzmann(s[o],o)) )); // ovvero vale 1 se l'energia nuova è inferiore, quindi visto che dopo lo confronto con un numero random tra 0 e 1 in quel caso verrebbe sempre accettato
      if(rnd.Rannyu() < A) {
        s[o] *= -1; // accetto la nuova configurazione in cui inverto lo spin
        accepted = accepted + 1; 
      }
        attempted = attempted + 1;
    }
    else { //Gibbs sampling
    // INCLUDE YOUR CODE HERE
      A = 1./(1. + exp(-beta*(Boltzmann(-1,o) - Boltzmann(1,o)) )); // probabilità di accettazione, coincide con la probabilità di estrarre spin up
      if(rnd.Rannyu() < A) // condizione di accettazione
        s[o] = 1 ; // se accetto devo fissare spin up perché quella sopra è la prob di estrarre spin up
      else s[o] = -1; // se rifiuto setto spin down
    }
  } // fine ciclo sul reticolo
} // fine funzione Move()


double Boltzmann(int sm, int ip) { // Calcolo dell'energia che UNO spin (di valore sm e in posizione ip) fornisce al valore totale del sistema, presa da sola sto considerando il doppio dell'interazione, quindi serve solo a calcolare la differenza di energia tra due stati in cui flippo un solo spin
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm; // uso le condizioni periodiche
  return ene;
}


void Measure() {
  // int bin;
  double u = 0.0, m = 0.0, u2 = 0.0, m2 = 0.0; // azzero energia e magnetizzazione

//cycle over spins
  for (int i=0; i<nspin; ++i) { // calcolo l'energia di TUTTO il reticolo considerando i primi vicini. Per non contare due volte le interazioni uso solo lo spin successivo a quello dell'indice di somma, non quello precedente
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
// INCLUDE YOUR CODE HERE
  	 m += s[i]; // calcolo la magnetizzazione
  }
  walker[iu] = u; // aggiorno l'array walker con il valore corrente misurato dell'energia
// INCLUDE YOUR CODE HERE
  u2 = u*u; // energia al quadrato per il calcolo del calore specifico
  m2 = m*m; // magnetizzazione al quadrato per il calcolo della suscettibilità magnetica
  walker[im] = m;
  walker[iu2] = u2;
  walker[im2] = m2;
} // fine funzione Measure()


void Reset(int iblk) { // Reset block averages
   if(iblk == 1) { // se è il primo blocco, per ogni variabile fisica annullo l'array rispettivo in cui sono salvate le media (se ho un solo blocco coincide con blk_av)
       for(int i=0; i<n_props; ++i) {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
// per tutti gli altri invece azzero le medie sui blocchi
   for(int i=0; i<n_props; ++i) {
     blk_av[i] = 0; // annullo il valore delle medie sui blocchi
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) { // Update block averages
   for(int i=0; i<n_props; ++i) {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0; // pari al numero di volte che chiamo Accumulate(), serve a calcolare la media sul blocco (è il numero degli step fatti a quel momento)
}


void Averages(int iblk) { // Print results for current block
   ofstream Ene, Heat, Mag, Chi;
   const int wd=15; // width
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl; // serve per capire quando ho raggiunto il 50%

	if (metro == 1) {
    	Ene.open("dati/Metropolis/T_fissata/output.ene.0",ios::app);
	} else 
    	Ene.open("dati/Gibbs/T_fissata/output.ene.0",ios::app);

    stima_u = blk_av[iu]/blk_norm/(double)nspin; // Energy per spin, divido per blk_norm che è il numero di volte in cui è chiamata la funzione Accumulate, ovvero il numero degli step, serve a calcolare la media su quel blocco
    glob_av[iu]  += stima_u; // in questa variabile aggiungo tutti i valori di tutti i blocchi, serve a calcolare la media tra i blocchi (media delle medie)
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close(); // divido per il numero di blocchi per trovare la media su quei blocchi considerati

// INCLUDE YOUR CODE HERE
// Calcolo calore specifico, mi serve la media dell'energia al quadrato
	if (metro == 1) {
    	Heat.open("dati/Metropolis/T_fissata/output.heat.0",ios::app);
	} else 
    	Heat.open("dati/Gibbs/T_fissata/output.heat.0",ios::app);
	
    stima_u2 = blk_av[iu2]/blk_norm; // Media sul blocco dell'energia quadro
    stima_c = beta*beta*(stima_u2  - stima_u*stima_u*nspin*nspin)/(double)nspin; // per utilizzare l'energia media devo moltiplicare per il numero di spin perché stima_u è l'energia media per spin
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

// Calcolo magnetizzazione
	if (metro == 1) {
    	Mag.open("dati/Metropolis/T_fissata/output.mag.0",ios::app);
	} else 
    	Mag.open("dati/Gibbs/T_fissata/output.mag.0",ios::app);

    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

	if (metro == 1) {
    	Chi.open("dati/Metropolis/T_fissata/output.chi.0",ios::app);
	} else 
    	Chi.open("dati/Gibbs/T_fissata/output.chi.0",ios::app);

    // Calcolo suscettibilità magnetica, mi serve la magnetizzazione al quadrato
    stima_m2 = blk_av[im2]/blk_norm/(double)nspin; // magnetizzazione al quadrato media sul blocco
    stima_x = beta*stima_m2;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();
	
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void) { // Stampa sul file config.final la configurazione finale del sistema
  ofstream WriteConf;
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
  rnd.SaveSeed();
}

int Pbc(int i) {  //Algorithm for Periodic Boundary Conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) { // dev std cumulativa
    if(iblk==1)
		return 0.0;
    else
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


void LastAverages(int iblk) { // Stampa il risultato finale, il valore in ingresso serve a capire quando sono arrivato alla fine.
   ofstream Ene, Heat, Mag, Chi;
   const int wd=15; // width
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl; // serve per capire quando ho raggiunto il 50%
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin, divido per blk_norm che è il numero di volte in cui è chiamata la funzione Accumulate, ovvero il numero degli step, serve a calcolare la media su quel blocco
    glob_av[iu]  += stima_u; // in questa variabile aggiungo tutti i valori di tutti i blocchi, serve a calcolare la media tra i blocchi
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);

	if (iblk == nblk) { // Stampo solo l'ultimo
		if (metro == 1) {
	    	Ene.open("dati/Metropolis/T_sequenza/output.ene.0",ios::app);
		} else 
	    	Ene.open("dati/Gibbs/T_sequenza/output.ene.0",ios::app);
	    
		Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << setw(wd) << temp << endl; // Scrivo anche la temperatura utilizzata
	    Ene.close();
	}
	
//===============================		
    stima_u2 = blk_av[iu2]/blk_norm; //Energy^2
    stima_c = beta*beta*(stima_u2  - stima_u*stima_u*nspin*nspin)/(double)nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);

	if (iblk == nblk) { // Stampo solo l'ultimo
		if (metro == 1) {
	    	Heat.open("dati/Metropolis/T_sequenza/output.heat.0",ios::app);
		} else 
	    	Heat.open("dati/Gibbs/T_sequenza/output.heat.0",ios::app);

	    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << setw(wd) << temp << endl;
	    Heat.close();
	}

//========================
    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
	
	if (iblk == nblk) { // Stampo solo l'ultimo
		if (metro == 1) {
	    	Mag.open("dati/Metropolis/T_sequenza/output.mag.0",ios::app);
		} else 
	    	Mag.open("dati/Gibbs/T_sequenza/output.mag.0",ios::app);

	    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << setw(wd) << temp << endl; // Scrivo anche la temperatura utilizzata
	    Mag.close();
	}

//==================
    stima_m2 = blk_av[im2]/blk_norm/(double)nspin; // Mag^2
    stima_x = beta*(stima_m2  - (h!=0.0 ? stima_m*stima_m*nspin : 0.00000));
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);

	if (iblk == nblk) { // Stampo solo l'ultimo
		if (metro == 1) {
	    	Chi.open("dati/Metropolis/T_sequenza/output.chi.0",ios::app);
		} else 
	    	Chi.open("dati/Gibbs/T_sequenza/output.chi.0",ios::app);

	    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << setw(wd) << temp << endl; // Scrivo anche la temperatura utilizzata
	    Chi.close();
	}
	
    cout << "----------------------------" << endl << endl;
}


void EnergyAverage(int iblk) { // Analogo alla funzione sopra ma stampa solo i valori finali dell'energia con relativo numero di blocchi utilizzato
   ofstream Ene;
   const int wd=15; // width
    
    cout << "Block number " << iblk << "/ " << nblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl; // serve per capire quando ho raggiunto il 50%
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin, divido per blk_norm che è il numero di volte in cui è chiamata la funzione Accumulate, ovvero il numero degli step, serve a calcolare la media su quel blocco
    glob_av[iu]  += stima_u; // in questa variabile aggiungo tutti i valori di tutti i blocchi, serve a calcolare la media tra i blocchi
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);

	if (iblk == nblk) { // Stampo solo l'ultimo
		if (metro == 1) {
	    	Ene.open("dati/Metropolis/Energy_Data_Over_Blocks.dat",ios::app);
		} else 
	    	Ene.open("dati/Gibbs/Energy_Data_Over_Blocks.dat",ios::app);
	    
		Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl; // Non serve scrivere il numero di blocchi usati perché basta guardare iblk che sarà pari a nblk perché stampo solo l'ultimo
	    Ene.close();
	}	
}


//================================================//
// funzione che stampa su file esterno un array di dati
void Print_File(const vector<double> & arr, string File_name) {
  ofstream out;
  out.open(File_name);

  // Controllo la corretta apertura del file
  if (!out.is_open()) {
    cerr << "Error: unable to open " + File_name << endl;
    exit(1);
  }

	int dim = arr.size();
	for(int i=0; i < dim; i++){
		out << arr[i] << endl;
	}
	
  out.close();
	
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
