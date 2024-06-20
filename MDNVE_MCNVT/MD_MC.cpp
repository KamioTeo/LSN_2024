/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Esercitazione 4-7 
// Dinamica Molecolare
// Implementazione del codice di Dinamica Molecolare, utilizzando l'algoritmo di Verlet (NVE) e Metropolis Monte Carlo (NVT)

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"
#include <vector>

using namespace std;

int main() { 
  Input("input.solid"); // inizializzo il sistema

int neq = 0; // Definizione variabile del numero di step per equilibrare il sistema (in Verlet 2000 step per solido e liquido e 100k per gas, in Metropolis 200)
if (iNVET == 1) { // se 1 è NVT MC Metropolis altrimenti Verlet NVE
  neq = 200; 
} else { // Se NVE Verlet
  neq = 2000; // per solido e liquido
	// neq = 100000; // per gassoso
}
//===========================================
// se voglio stampare su file esterno i valori di temperatura ed energia istantanei DURANTE l'equilibrazione
  if (PrintEqValues == 1) { 
	vector<double> Temperature(neq); // vettore che conterrà la temperatura calcolata ad ogni step per l'equilibrazione
	vector<double> Energy(neq); // Per contenere U/N nel caso di algoritmo di metropolis

	// Salvo la temperatura al variare dei passi, serve a sapere quando ho raggiunto l'equilibrio e come modificare i parametri per garantire che sia la temperatura desiderata
	for(int istep=1; istep <= neq; istep++) { // per ogni incremento temporale (step)
    	Move(); // Movimento delle molecole, dato da algoritmo di Verlet o Metropolis
    	Measure(); // Misuro nuovamente temperatura, energia potenziale, cinetica e totale a questo step
		Temperature[istep-1] = walker[it]; // Salvo la temperatura ad ogni step
		//if (iNVET == 1) // salvo l'andamento dell'energia potenziale per particella
			Energy[istep-1] = walker[iv]/(double)npart;
	}
	// salvo i valori dell'equilibrazione su file esterno
	Print_File(Temperature, "./dati/Equilibration_T.dat");
	Print_File(Energy, "./dati/Equilibration_U-N.dat");
	
return 0; // fermo qui il programma, questa parte la uso solo per salvare i risultati durante l'equilibrazione per poterli plottare
}

// se voglio stampare i valori dell'energia istantanea DOPO l'equilibrazione
  if (PrintEnergy == 1 and iNVET == 1) {
	cout<<"Equilibration"<<endl;
	for(int i = 0; i < neq; ++i) // equilibrio
		Move();

	vector<double> Energy(nstep); // Per contenere U/N

	for(int istep=1; istep <= nstep; istep++) { // per ogni incremento temporale (step)
    	Move(); // Movimento delle molecole, dato da algoritmo di Verlet o Metropolis
    	Measure(); // Misuro nuovamente temperatura, energia potenziale, cinetica e totale a questo step
		Energy[istep-1] = walker[iv]/(double)npart;
	}
	// salvo i valori DOPO equilibrazione su file esterno. Chiamo il file così perché è il dataset da cui calcolerò l'autocorrelazione, non è esso stesso l'autocorrelazione
	Print_File(Energy, "./Autocorrelation_U-N.dat");
	
return 0; // fermo qui il programma, questa parte la uso solo per salvare i risultati durante l'equilibrazione
}
	
	
//===========================================
// Ora inizio l'analisi: nella parte precedente ho trovato il numero di step, il numero di blocchi e il valore di delta ottimali per l'analisi. Qui equilibrio il sistema con quel numero di step
// è molto più efficiente iniziare la simulazione a partire dai file di configurazioni di posizione e velocità raggiunte in precedenza, perchè non devo rifare tutte le operazioni di equilibrazione
	cout<<"Equilibration"<<endl;
	for(int i = 0; i < neq; ++i) // Per equilibrare uso lo stesso numero di passi controllati prima
		Move();

	// se voglio stampare i valori di accettazione (codice MC NVT) e fermare il programma
	if (PrintAcceptance == 1) {
		cout << "Acceptance after " << neq << " step and delta: " << delta << " is:" << accepted/attempted << endl;
	return 0;
	}
	
// il numero di blocchi da utilizzare è passato da file di input
  int nconf = 1; // setto a 1 questa variabile all'inizio
  for(int iblk=1; iblk <= nblk; iblk++) { // Per ogni blocco (nblk lo legge dal file di input)
    Reset(iblk);   // Resettto le medie sui blocchi
    for(int istep=1; istep <= nstep; istep++) { // per ogni incremento temporale (step)
      Move(); // Movimento delle molecole, dato dall'algoritmo di Verlet o Metropolis
      Measure(); // Misuro nuovamente temperatura, energia potenziale, cinetica e totale a questo step
      Accumulate(); // Update block averages, sommo il valore calcolato a questo step dei 4 parametri a quelli calcolati agli step precendenti
      
		if(istep%10 == 0) { // ogni 10 step
//        ConfXYZ(nconf); // Write actual configuration in XYZ format // Commented to avoid "filesystem full"! Stampa la posizione x, y, z per ogni particella
        nconf += 1;
      }
		
    } // Fine step totali

    Averages(iblk);   // Print results for current block, alla fine degli step
  } // fine blocchi
	
  ConfFinal(); // Write final configuration, posizione finale in config.out e velocità finale in velocity.out
	
  return 0;
}

//===========================================================//
// Inizializza il sistema, caricando i parametri di input, le posizioni e velocità iniziali da file esterno e il seed del Random Number Generator
void Input(string FileToOpen) {
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

// Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ; // Carico i numeri da file
  Primes.close();

//Read input informations
  ReadInput.open(FileToOpen);

// ATTENZIONE: densità, raggio di cutoff, temperatura dipendono dallo stato di aggregazione della sostanza, i valori riportati nei commenti sono quelli nell'input.in file
  ReadInput >> iNVET; // Carico il tipo di analisi che voglio fare, se 0 è Molecular Dynamichs con algoritmo di Verlet (NVE, ovvero sono constai il numero di particelle N, il volume V e l'energia E), se 1 è Monte Carlo con algoritmo di Metropolis (NVT, costanti particelle N, volume V e temperatura T). Le medie temporali sulle traiettorie di quantità termodinamiche di un macrostato coincidono con le medie d'insieme sui microstati
  ReadInput >> restart; // se 1 comincia dalla configurazione caricata da file esterno invece che da quella di cristallo perfetto

  if(restart) Seed.open("seed.out"); // se devo partire da una configurazione esistente devo anche continuare le simulazioni usando lo stesso seed
  else Seed.open("seed.in"); // se sono già partito allora apro il file con i seed già creati
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2); // setta questi valori a dei parametri interni dell'istanza random, serve per generare i numeri casuali
  Seed.close();

  ReadInput >> temp; // Temperatura in unità ridotte (vale 1.1 liquido, 0.8 solido)
  beta = 1.0/temp; // (beta = 1/T)
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart; // Numero di particelle in un volume di riferimento (qui vale 108)
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho; // Densità di particelle per unità di volume (qui vale 0.8 liquido, 1.1 solido)
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho; // il volume di riferimento è semplicemente Nparticelle / densità
  box = pow(vol,1.0/3.0); // è la lunghezza di un lato del volume di riferimento (radice cubica del volume)
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut; // (qui 2.5 liquido, 2.2 solido) è la distanza da una particella fino a dove considero l'interazione con altre particelle, ovvero tutte quelle nella sfera di questo raggio
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta; // (qui 0.0005 per NVE e 0.112 per NVT) per NVT è il dt^* d'integrazione, ovvero l'intervallino di tempo per integrare, in unità ridotte (dt^* = dt*sqrt(epsilon / (mass*sigma^2) ), per NVE è il delta per estrarre la posizione con transizione uniforme
  ReadInput >> nblk; // numero di blocchi (20)
  ReadInput >> nstep; // numero di step da calcolare del movimento delle particelle (2000)
  ReadInput >> PrintAcceptance; // 1 se voglio stampare il file con i valori di accettazione
  ReadInput >> PrintEqValues; // 1 se voglio stampare i risultati di temperatura ed energia durante l'equilibrazione
  ReadInput >> PrintEnergy; // se 1 stampa i valori dell'energia per NVT dopo l'equilibrazione
	
  cout << "The program also perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

// Ora che ho caricato tutti i parametri posso calcolare le correzioni di coda al potenziale di Lennard-Jones e alla pressione; le calcolo subito qui perché sono costanti

  vtail = (8./3.)*M_PI* rho *( 1./(3.*pow(rcut,9)) - 1./pow(rcut,3) ) * npart;
  ptail = (32./3.)*M_PI* rho *( 1./(3.*pow(rcut,9)) - 1./(2.*pow(rcut,3)) ) * rho;

  cout << "Tail correction for the potential energy = " << vtail/npart << endl;
  cout << "Tail correction for the pressure = " << ptail << endl; 

// osservazione sulla correzione alla pressione: il fattore moltiplicativo a destra è dato dal fatto che la correzione è w/3N, ma per trovarla in termini di pressione devo dividere w per 3V, ovvero complessivamente moltiplicare per rho
	
// Prepare arrays for measurements
  iv = 0; // Potential energy 
  it = 1; // Temperature
  ik = 2; // Kinetic energy
  ie = 3; // Total energy
  ip = 4; // Pressure
  ig = 5; // Indice per la distribuzione a coppie
  n_props = 5+nbins; // Number of observables
	
  bin_size = (box/2.) / nbins; // Larghezza bin, devono dividere metà della larghezza del volume
	
// Read initial configuration, se voglio partire da quella finale di una simulazione precedente (salvata su file esterno)
  cout << "Read initial configuration" << endl << endl;
  if(restart) {
    ReadConf.open("config.out"); // Posizioni iniziali
    ReadVelocity.open("velocity.out"); // Velocità iniziali
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i]; // carico le velocità iniziali (nelle 3 direzioni) per ognuna delle npart=108 particelle, la uso per calcolare la posizione al tempo precedente
  }
  else { // altrimenti se non carico da file esterno, creo tutte le velocità iniziali di ogni particella
    ReadConf.open("config.in");  // carico  le posizioni iniziali delle particelle
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
	  
    for (int i=0; i<npart; ++i) { // setto una velocità del centro di massa di ogni particella da distribuzione di Maxwell-Boltzmann (quindi gaussiana con dev std pari alla radice della temperatura)
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
 // Calcolo il vettore velocità totale (velocità di drift) da tutte le particelle. Quindi sommo tutte le velocità di tutte le particelle nelle 3 direzioni
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
	  
// per ogni direzione calcolo la velocità di drift media per particella (velocità totale diviso il numero di particelle)
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  
    double sumv2 = 0.0, fs;
// tolgo la velocità di drift per particella alla rispettiva velocità, ovvero mi metto nel sdr del CM
    for (int i=0; i<npart; ++i) {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
	// somma di tutti i moduli quadri delle singole velocità, serve a calcolare l'energia cinetica media, che infatti dipende dalla somma delle singole velocità delle particelle al quadrato
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
	  
// velocità quadro media per particella
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor, è pari a sqrt(massa/epsilon)
// Dimostrazione: voglio trovare il coefficiente fs t.c. v = fs * v^* (perché finora ho le velocità in unità ridotte e voglio trovarle non ridotte). Questa vale sia per la media che per le singole velocità. fs lo calcolo a partire dalla media delle velocità, più nello specifico a partire dalla media dei quadrati delle velocità perché così uso i legami conl'energia cinetica:

// <E_k> = 0.5 * m * <v^2> = 0.5 * m * <v^*^2> * fs^2
// ma E_k si può scrivere anche come:
// <E_k> = (3/2) * k_b * T * (epsilon/epsilon) = (3/2) * epsilon * T^*
// metto a sistema le due e ricavo fs:
// 0.5 * m * <v^*^2> * fs^2 = (3/2) * epsilon * T^* -> 
// fs^2 = 3 * epsilon * T^* / (m * <v^*^2>)
// ma sono in unità ridotte quindi epsilon = m = 1, quindi:
// fs = sqrt(3 * T^* /  <v^*^2>)
	  
    cout << "velocity scale factor: " << fs << endl << endl;
	  
    for (int i=0; i<npart; ++i) {
// converto da velocità ridotte a velocità normali
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
	  
  } // fine crezione velocità

// carico dal file config.in le posizioni iniziali delle particelle
  for (int i=0; i<npart; ++i) {
    ReadConf >> x[i] >> y[i] >> z[i];
// uso l'algoritmo delle Periodic Boundary Condition per sfruttare la periodicità delle celle del reticolo
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
	
  ReadConf.close();

  for (int i=0; i<npart; ++i) {
    if(iNVET) { // se è Monte Carlo NVT: N particelle, Volume e Temperatura costanti, la posizione nuova proposta va accettata, quindi per ora setto la posizione iniziale come quella letta da file esterno
	    xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else { // Molecular Dynamichs (NVE, Energy costant)
// calcolo posizione precedente (infatti sottraggo) a partire dalla posizione e dalla velocità
      xold[i] = Pbc(x[i] - vx[i] * delta); // posizione x al tempo t-dt
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  } // fine set old position
  
// Evaluate properties of the initial configuration
// Misura di E potenziale, E cinetica, E totale, temperatura all'inizio
  Measure(); // setta i valori dell'array walker stampato di seguito
	
// Recall sugli indici:
  // iv = 0; // Potential energy 
  // it = 1; // Temperature
  // ik = 2; // Kinetic energy
  // ie = 3; // Total energy
  // ip = 4; // Pressure
  // ig = 5; // Indice per la distribuzione a coppie
  // n_props = 5+nbins; // Number of observables

// Print initial values for measured properties per particle
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure			= " << walker[ip] << endl;
  
cout << endl << endl << "----------------------------" << endl << endl;
	
  return;
} // Fine definizione funzione Input()

//===============================================//
// Movimento delle molecole
void Move() {
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

// Monte Carlo (NVT) move
  if(iNVET) {
    for(int i=0; i<npart; ++i) { // per ogni particella
    // Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart); // estraggo una particella a caso

    // Old
      energy_old = Boltzmann(x[o],y[o],z[o],o); // misuro l'energia di questa particella (data da potenziale additivo di coppia di Lennard Jones con tutte le altre particelle nella cella attuale)

    // New, la nuova posizione proposta la estraggo casualmente in un range di +/- delta/2, dove delta va settato in modo da avere accettazione del 50%
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    // Calcolo l'energia nella nuova posizione proposta
      energy_new = Boltzmann(x[o],y[o],z[o],o);

    // Metropolis test, la probabilità da campionare è quella di Boltzmann
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu()) { // accetto la proposta
      // Update
        xold[o] = x[o]; // quindi la posizione proposta diventa ora la posizione vecchia da cui far partire la prossima simulazione
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0; // conto che ho accettato
      } else { // altrimenti rifiuto e la posizione vecchia resta quella da cui partirò
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0; // conto i rifiuti
    }
  } else { // Molecular Dynamics (NVE) move
    double fx[m_part], fy[m_part], fz[m_part]; // Arrays che contengono la forza totale (nelle 3 direzioni) subìta da ogni particella. m_part è il numero di particelle nella cella fondamentale, variabile globale definita nel .h

    for(int i=0; i<npart; ++i) { // Force acting on particle i
      fx[i] = Force(i,0); // 0, 1, 2 sono le direzioni in cui calcolo le forze con la funzione Force
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }
// ora ho calcolato tutte le forze su ogni particella
// OSSERVAZIONE: l'accelerazione di una particella è anche la forza totale su di essa perché sono in unità ridotte
    for(int i=0; i<npart; ++i) { // Verlet integration scheme
// RECALL: r(t+dt) = 2*r(t) - r(t-dt) + dt^2 * a(t), ovvero mi serve la posizione a questo step, la posizione precedente e la forza totale su quella particella
// xnew è r_x(t+dt), x[i] è r_x(t), xold è r_x(t-dt)
      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

// sovrascrivo le velocità con quella nuova ma ricavata dall'algoritmo, ovvero calcolo v(t) grazie alla posizione a t+dt e alla posizione a t-dt. Al primo step coincide con v caricata da input file (infatti ho usato quella velocità per calcolare x a t-dt, ma qui sto riusando xold)
      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
// ora la nuova posizione a t-dt è sovrascritta con quella a t
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
// mentre la nuova posizione a t è sovrascritta con quella a t+dt
      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;
// ovvero globalmente sono andato avanti di uno step temporale
      accepted = accepted + 1.0; // alla fine calcolerò il rapporto tra accepted/attempted
      attempted = attempted + 1.0; // nell'algoritmo di Verlet è sempre 100%, mentre in MC no perché sommo accepted solo in certi casi
    }
  } // fine Move MD
  return;
} // fine Move()

//===============================================//
// serve solo all'algoritmo MC
// Calcolo energia potenziale totale in una sfera di raggio cutoff centrata nelle coordinate xx, yy, zz
double Boltzmann(double xx, double yy, double zz, int ip) {
  double ene = 0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i) {
    if(i != ip) {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr); // modulo distanza dalla particella ip-esima

      if(dr < rcut) { // Potenziale di Lennard Jones
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  } // fine ciclo particelle

  return 4.0*ene + vtail; // aggiungo la tail correction anche qui, anche se in realtà non serve perché poi devo fare la differenza tra due valori di energia e vtail si semplifica
}

//===============================================//
// Calcolo della forza subìta su una particella al centro di una sfera di raggio cutoff
// in ingresso dò l'indice della particella centrale (ip) e la direzione in cui calcolare la forza (idir)
double Force(int ip, int idir) { // Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i) { // calcolo la distanza tra la particella centrale e tutte le altre
    if(i != ip){ // guardo solo le altre particelle, non considero quella centrale
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr); // calcolo modulo distanza

      if(dr < rcut) { // considero solo le particelle nella sfera cutoff
	// derivata direzionale in direzione idir come prodotto scalare <r_dir ; -\nabla V>
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r), calcolo la forza totale, infatti viene il pot di LJ moltiplicato per r^-2, raccogliendo 48 viene un 1/2 davanti al termine in r^-6
      }
    }
  }
 
  return f; // Forza totale in direzione idir sulla particella ip
}

//===============================================//
// Misura di E potenziale, E cinetica, E totale, temperatura, pressione
void Measure() { // Properties measurement
  double v = 0.0, kin=0.0, w=0.0; // viriale
  double vij, wij; // addendi
  double dx, dy, dz, dr;

// Devo azzerare i valori dell'istogramma per il calcolo di g
  for ( int i = ig; i < ig+nbins; i++)
	  walker[i] = 0;

  // cycle over pairs of particles
  for (int i=0; i<npart-1; ++i) {
    for (int j=i+1; j<npart; ++j) {
    // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr); // elemento di distanza infinitesimo
    // controllo che sia dentro il raggio di cutoff, solo in questo caso ne conto il contributo
      if(dr < rcut) {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6); // contributo ij-esimo al Potenziale di Lennard-Jones in unità ridotte (epsilon = sigma = 1) (manca un 4 moltiplicativo)
        v += vij; // Totale sommatoria del V_LJ
		    wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
		    w += wij;
      }

    // update the histogram of g(r)
      idr = ig + int(dr/bin_size); // Questa è la distanza nell'istogramma (quindi il numero di bin) alla quale si trova dr. Infatti, dr/bin_size è il numero di bin che ci vogliono per raggiungere una distanza dr, poi aggiungo l'indice ig=5 di partenza per traslare
      if (idr <= ig+nbins) // conto solo le particelle che stanno nell'istogramma, ovvero che distano tra loro una distanza inferiore a L/2
        walker[idr] += 2; // ogni contributo è doppio perché sto sommando su i e j
		
    }
  } // fine ciclo sulle coppie di particelle

  for (int i=0; i<npart; ++i)
	  kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); // E_kin TOTALE in unità ridotte

/* RECALL:
	iv = 0; // Potential energy 
  	it = 1; // Temperature
  	ik = 2; // Kinetic energy
  	ie = 3; // Total energy
    ip = 4; // Pressure
    ig = 5; // Indice per la distribuzione a coppie
  	n_props = 5+nbins; // Number of observables
*/
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperatura, infatti E_k = (3/2) * N * K_b * T, scritto però in unità ridotte

  walker[iv] = 4.0 * v; // Potential energy, moltiplico il 4 che mancava al valore trovato prima per V_LJ
  walker[ip] = rho*walker[it] + 1./vol*16.*w;
		
  if (iNVET) { // se è Metropolis NVT devo aggiungere le tail corrections
    walker[iv] += vtail;
    walker[ip] += ptail;
  }

  walker[ie] = walker[iv] + kin;  // Total energy; potenziale più cinetica

  return;
}

//===============================================//
//Reset block averages
void Reset(int iblk) {
   if(iblk == 1) {
       for(int i=0; i<n_props; ++i) { // n_props = 5+nbins, numero di osservabili, ma glob_av ha dimensione m_props=1000, quella massima pari al numero massimo di variabili misurabili
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i) {
     blk_av[i] = 0; // azzero le medie sui blocchi, anche questo ha dim m_props=1000
   }
	
   blk_norm = 0; // serve per il calcolo della media sui blocchi
   attempted = 0;
   accepted = 0;
}

//===============================================//
//Update block averages, eseguito per ogni passo, a blocco fissato
void Accumulate(void) {
// per onguno dei 5 parametri da calcolare:
   for(int i=0; i<n_props; ++i) { // n_props = 5+nbins
     blk_av[i] = blk_av[i] + walker[i]; // allo step successivo sommo ognuno dei parametri con il valore precedente
   }
   blk_norm = blk_norm + 1.0; // block_norm parte da 0 all'inizializzazione (vedi funzione subito sopra)
}

//===============================================//
//Print results for current block
void Averages(int iblk) { // funzione chiamata ad ogni blocco successivo, dopo aver fatto tutti i passi
	cout << "Averages called" << endl;
   ofstream Epot, Ekin, Etot, Temp, Press, G_r;
   const int wd=12; // width, larghezza campo output
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("./output_epot.dat",ios::app);
    Ekin.open("./output_ekin.dat",ios::app);
    Temp.open("./output_temp.dat",ios::app);
    Etot.open("./output_etot.dat",ios::app);
    Press.open("./output_pressure.dat",ios::app);
    G_r.open("./output_g_r.dat",ios::app);
	
// RECALL:
  // iv = 0; // Potential energy 
  // it = 1; // Temperature
  // ik = 2; // Kinetic energy
  // ie = 3; // Total energy
  // ip = 4; // Pressure
  // ig = 5; // Indice per la distribuzione a coppie
  // n_props = 5+nbins; // Number of observables

// Visto che Averages è eseguito fuori dal ciclo degli step, blk_av qui è la somma di TUTTI i parametri ottenuti ad ogni step (a fine ciclo), quindi per trovare la media basta normalizzare sul numero di step per blocco (blk_norm) e per il numero di particelle (uso il fatto che la media temporale è circa la media spaziale sui microstati)
    stima_pot = blk_av[iv]/blk_norm/(double)npart; // Potential energy
    glob_av[iv] += stima_pot; // media cumulativa
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm/(double)npart; // Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; // Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; // Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    stima_press = blk_av[ip]/blk_norm; // Pressure
    glob_av[ip] += stima_press;
    glob_av2[ip] += stima_press*stima_press;
    err_press=Error(glob_av[ip],glob_av2[ip],iblk);
  
  // stima di g(r) (esercitazione 6), è una funzione quindi devo calcolarla per ogni bin dell'istogramma
    double r_norm, r; // devo calcolare anche la normalizzazione e la posizione r
	  double rad[nbins];
    for(int i = ig; i < n_props; i++) { // ciclo sui bin dell'istogramma
      rad[i-ig] = (i+1-ig)*bin_size; // la distanza nel bin k-esimo è data da quel numero di bin moltiplicato per la larghezza dei bin. L'indice dell'istogramma non parte da 1 quindi devo traslarlo. Ad esempio il primo i=ig è ad una distanza r=bin_size
      r = rad[i-ig];
		  r_norm = rho*(double)m_part*4./3.*M_PI*( pow( (r +bin_size), 3 ) - pow(r,3) );     // la normalizzazione è data da rho*N*deltaV, con deltaV = 4/3*Pi*[(r+dr)^3-r^3]       

    // err_g è un array, salvo i valori della stima di g nello stesso array in cui salvo le altre misure
	   stima_g = blk_av[i]/((double)blk_norm*r_norm); //(double)iblk;
    //cout<<blk_norm<<endl;
      glob_av[i] += stima_g;
      glob_av2[i] += stima_g*stima_g;
      err_g[i-ig] = Error(glob_av[i],glob_av2[i],iblk); 
    }

// Stampo su file i risultati di ogni blocco, cumulando i blocchi e l'errore cumulato
// Potential energy per particle, uso setwidth per impostare la larghezza del campo di output e garantire un corretto allineamento
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
// Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
// Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
// Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
// Pressure
    Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_press << endl;

// g(r), stampo l'istogramma per ogni blocco
	G_r << "___BLOCK NUMBER___" << iblk << endl;
	for (int i = ig; i < n_props; i++) {
		G_r << setw(wd) << iblk << setw(wd) << rad[i-ig] << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_g[i-ig] << endl;
	}
	G_r << endl;
	
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Press.close();
	G_r.close();
}

//===============================================//
// stampa posizioni e velocità finali
void ConfFinal(void) {
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("./config.out");
  WriteVelocity.open("./velocity.out");
  for (int i=0; i<npart; ++i) {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

//===============================================//

void ConfXYZ(int nconf) { //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//===============================================//
// Algoritmo delle Periodic Boundary Conditions, per un volume cubico di lato L=box
double Pbc(double r) {
    return r - box * rint(r/box);
}

//===============================================//
// Deviazione standard per divisione a blocchi
double Error(double sum, double sum2, int iblk) {
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
