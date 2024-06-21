// Salesman Problem with Genetic Algorithm and MPI
/* ATTENZIONE:
Se dà questo errore:
-> "The value of the MCA parameter "plm_rsh_agent" was set to a path
that could not be found:

  plm_rsh_agent: ssh : rsh

-> da terminale setta questa variabile d'ambiente:
-> export OMPI_MCA_plm_rsh_agent=
*/
#include "functions.h"

// opzioni per il debugging
bool globalPrint = false; // se stampare o no tutte le configurazioni
bool globalPrintMutation = false; // se stampare o no tutte le configurazioni per ogni mutazione

// probabilità di mutazione e di crossover
double p_mut = 0.20;
double p_cross = 0.95;
// probabilità di sopravvivenza, ovvero prob. di lasciare il genitore se il figlio dovesse essere più lungo (al contrario, se è più corto lo tengo sempre)
double p_surv = 0.5;

// Stampa una sequenza
void PrintPos(vector<Posizione> pos) {
	int dim = pos.size();
	for (int i=0; i<dim; i++) {
		cout << pos[i].index << " " << setprecision(3) << pos[i].x << " " << pos[i].y << endl;
	}
	cout << "-----------" << endl;
}

// stampa una popolazione di sequenze
void PrintPop(vector<Individual> pop) {
	int dim = pop.size();
	for (int i=0; i<dim; i++) {
		cout << "---" << i+1 << "---" << endl;
		PrintPos(pop[i].pos);
		cout << "dist:" << endl;
		cout << pop[i].dist << endl;
		cout << "-----------" << endl;
	}
}

int main(int argc, char *argv[]) {
	// Definisco i parametri MPI
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if ( size < 3 || size > 10) {
		if (rank == 0) {
			cout << "The number of cores must be between 3 and 10" << endl;
		}
		MPI_Finalize();
		exit(0);
	} else if (rank == 0) {
		cout << "_______________TSP_______________" << endl;
		cout << "The Salesman Problem code using " << size << " cores." << endl << endl;
	}

	// Per ogni core setto un seed diverso
    for(int i=0; i < size; i++) {
        if(rank == i) {
			rnd.RandomMPI( "../RandomNumberGenerator/primes32001.in", "../RandomNumberGenerator/seed.in", rank);
		}
    }
	
	// leggo il numero di città da usare, il tipo di posizionamento e la dimensione della popolazione
	Input(rank);

	// Per tutti i core creo il vector che ospiterà le posizioni delle città
	vector<Posizione> x(N_cities+1);
	if (rank == 0) {
		if (shape == 2 || shape == 3 || shape == 4) { // carico le città da file esterno
			x = CityLoader(shape);
		} else { // creo le città random
			x = CityPlacer(N_cities, shape);
		}
		CheckBoundOk(x);
	}
	
	// Spedisco il risultato agli altri core dal nodo 0 (altrimenti ogni core avrebbe cità in posizioni diverse), ma prima devo specificare il data type della struct Posizione
	MPI_Datatype MPI_Posizione;
	// così interpreto il datatype 'int' dell'index come float, è tutto pi+ù compatto ma non è proprio corretto per MPI (anche se funziona)
	// MPI_Type_contiguous(3, MPI_FLOAT, &MPI_Posizione);
	int blocklengths[] = {1, 1, 1};
	MPI_Aint displacements[] = {offsetof(Posizione, index), offsetof(Posizione, x), offsetof(Posizione, y)};
	MPI_Datatype types[] = {MPI_INT, MPI_FLOAT, MPI_FLOAT};
	MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_Posizione);
	// salvo il datatype personalizzato
	MPI_Type_commit(&MPI_Posizione);
	// In Bcast devo specificare il puntatore al buffer di dati da spedire, non il vector stesso
	MPI_Bcast(&x.front(), x.size(), MPI_Posizione, 0, MPI_COMM_WORLD);
	// Libera il tipo di dati MPI personalizzato quando hai finito
	// MPI_Type_free(&MPI_Posizione);

	// controllo che le posizioni siano correttamente inviate a tutti i core
	// for (int i=0; i<N_cities; i++) {
	// 	cout << "rank " << rank << x[i].index <<  " "<< fixed << setprecision(5) << x[i].x  << " " << x[i].y << endl;
	// }
	// cout << "RANK " << rank << " END" << endl;
	// MPI_Finalize();
		
	// return 0;

	// stampo le posizioni ottenute, uguali per tutti i core quindi lo faccio solo per il rank 0
	if (rank == 0) {
		PrintPositionsFile(x, shape);
	}
	
	// creo la popolazione (prima generazione), da permutazioni casuali della configurazione iniziale (Check BC svolto nella funzione)
	vector<Individual> population = GeneratePopulation(dim_pop, x);

	// Riordino la popolazione in ordine crescente di distanza (serve per il crossover)
    sort(population.begin(), population.end(), compareByDist);

	if (globalPrint == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		PrintPop(population);
	}
	cout << "--------RANK " << rank << "--------" << endl;

	// Vector che conterrà tutte le generazioni
	vector<vector<Individual>> AllGenerations(N_gen);
	// Salvo la prima generazione
	AllGenerations[0] = population;
	
	// Vector che conterrà la distanza minima per ogni generazione
	vector<Individual> AllminSeq(N_gen);
	// Salvo la distanza minima della prima generazione
	AllminSeq[0] = MinimumDistance(population);
	
	if (globalPrint == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		cout << "0 minimum distance: " << AllminSeq[0].dist << endl;
}

	for (int gen=1; gen<N_gen; gen++) {
		// creo un nuovo vector per ospitare la nuova generazione
		vector<Individual> new_population(dim_pop);
		// creo la nuova generazione a partire da quella vecchia (Check BC svolto nella funzione)
		NewGeneration(new_population, population, rank);
		// salvo la sequenza minima della nuova generazione
		AllminSeq[gen] = MinimumDistance(new_population);
		
		if (globalPrint == true) {
			cout << "==== RANK " << rank << " ====" << endl;
			// stampo la nuova popolazione
			cout << "-----Generation_" << gen+1 << "----" << endl;
			PrintPop(new_population);
			cout << "--------" << endl;
		}

		// Riordino la nuova generazione
    	sort(new_population.begin(), new_population.end(), compareByDist);

		// Ogni N_migr generazioni scambio i primi elementi delle popolazioni (ovvero i migliori individual) tra i core
		if (gen % N_migr == 0) {
		  // per tutti i core creo il vector che ospiterà le coppie di rank che comunicheranno
		  vector<pair<int,int>> coppie;
		  if (rank == 0) {
		    // genero le coppie di rank casuali che comunicheranno, la dimensione dipende da size e dalla sua parità (ma non è size)
			coppie = GeneraCoppie(size);
		  }

		  // devo dire agli altri core (diversi da 0) quanto è la dimesione del vector coppie
		  int numCoppie = coppie.size();
		  // Broadcast del numero di coppie
		  MPI_Bcast(&numCoppie, 1, MPI_INT, 0, MPI_COMM_WORLD);
		  // setto il size del vector per tutti i rank che non sono 0
		  if (rank != 0) {
			coppie.resize(numCoppie);
		  }
	
		  // Broadcast delle coppie di rank
		  MPI_Bcast(&coppie.front(), numCoppie * 2, MPI_INT, 0, MPI_COMM_WORLD);
			
		  if (rank==0) {
		    for (int i=0; i<numCoppie; i++) {
			  cout << "gen " << gen << " coppia " << i+1 << ": " << coppie[i].first << " - " << coppie[i].second << endl;
			}
		  }

		  // ora invece di spedire un Individual che ha una struttura complessa, spedisco un vector di Posizioni e la sua distanza in momenti separati, per poi riunirli
		  vector<Posizione> recvPosizione(N_cities+1);
		  float recvDist;

/*
		  if (rank==1) {
			cout << "rank 1 to SEND individual: " << endl;
		    PrintPos(new_population[0].pos);
		cout << "dist: " << new_population[0].dist << endl;
		  }
		  if (rank==3) {
			cout << "rank 3 to SEND individual: " << endl;
		    PrintPos(new_population[0].pos);
		cout << "dist: " << new_population[0].dist << endl;
		  }
			if (rank==0) {
			cout << "rank 0 to SEND individual: " << endl;
		    PrintPos(new_population[0].pos);
		cout << "dist: " << new_population[0].dist << endl;
		  }
			if (rank==2) {
			cout << "rank 2 to SEND individual: " << endl;
		    PrintPos(new_population[0].pos);
cout << "dist: " << new_population[0].dist << endl;
		  }
*/
			
		  // Per ogni coppia, il primo invia la sua configurazione migliore al secondo
		  for (int i=0; i<numCoppie; i++) {
			// Il primo rank invia al secondo
			if (rank==coppie[i].first) {
				MPI_Send(&new_population[0].pos.front(), N_cities+1, MPI_Posizione, coppie[i].second, 0, MPI_COMM_WORLD);
				MPI_Send(&new_population[0].dist, 1, MPI_FLOAT, coppie[i].second, 0, MPI_COMM_WORLD);
			} else if (rank==coppie[i].second) {
			  // OSS: questo metodo è bloccante
		      MPI_Recv(&recvPosizione.front(), N_cities+1, MPI_Posizione, coppie[i].first, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

MPI_Recv(&recvDist, 1, MPI_FLOAT, coppie[i].first, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		  }

		  // Analogo ma per ogni coppia il secondo invia la sua configurazione migliore al primo
		  for (int i=0; i<numCoppie; i++) {
			// Il primo rank invia al secondo
			if (rank==coppie[i].second) {
				MPI_Send(&new_population[0].pos.front(), N_cities+1, MPI_Posizione, coppie[i].first, 0, MPI_COMM_WORLD);
				MPI_Send(&new_population[0].dist, 1, MPI_FLOAT, coppie[i].first, 0, MPI_COMM_WORLD);
			} else if (rank==coppie[i].first) {
			  // OSS: questo metodo è bloccante
		      MPI_Recv(&recvPosizione.front(), N_cities+1, MPI_Posizione, coppie[i].second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
			  MPI_Recv(&recvDist, 1, MPI_FLOAT, coppie[i].second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		  }

		  // individual vuoto su cui salverò quelli inviati dai vari rank
		  Individual RecvIndividual;
		  RecvIndividual.pos = recvPosizione;
		  RecvIndividual.dist = recvDist;

/*
		  if (rank==2) {
			cout << "rank 2 RECV individual: " << endl;
		    PrintPos(RecvIndividual.pos);
			  cout << "dist: " << RecvIndividual.dist << endl;
		  }
		  if (rank==0) {
			cout << "rank 0 RECV individual: " << endl;
		    PrintPos(RecvIndividual.pos);
			 cout << "dist: " << RecvIndividual.dist << endl; 
		  }
		  if (rank==1) {
			cout << "rank 1 RECV individual: " << endl;
		    PrintPos(RecvIndividual.pos);
			  cout << "dist: " << RecvIndividual.dist << endl;
		  }
		  if (rank==3) {
			cout << "rank 3 RECV individual: " << endl;
		    PrintPos(RecvIndividual.pos);
			  cout << "dist: " << RecvIndividual.dist << endl;
		  }
*/
	
		// ora tutti i rank (e l'ultimo se size è pari) hanno salvato in RecvIndividual la configurazione migliore rivevuta da altri rank
		// lo aggiungo al primo posto della popolazione
    	new_population.insert(new_population.begin(), RecvIndividual);
    	// Elimino l'ultimo individuo (dist più elevata) per riavere la dimensione corretta
    	new_population.pop_back();
		// oss: in generale la nuova popolazione non è più in ordine, ma prediligo il fatto di aver scambiato delle configurazioni tra processi
		}
		
		// Salvo la nuova generazione
		AllGenerations[gen] = new_population;
		// la nuova generazione diventa quella vecchia nel nuovo ciclo
		population = new_population;

		if (globalPrint == true) {
			cout << "==== RANK " << rank << " ====" << endl;
			cout << gen << " minimum distance: " << AllminSeq[gen].dist << endl;
		}
	}

// stampo tutte le distanze per ogni generazione
PrintDistancesFile(AllGenerations, shape, rank);
// stampo la sequenza con la distanza minima
PrintMinSequence(AllminSeq, shape, rank);

// Salvo il seed di ogni rank in un file diverso
rnd.SaveSeedMPI(rank);

cout << "RANK " << rank << " END" << endl;
MPI_Finalize();
	
return 0;

}

//===========================================================
void Input(int rank) {
	ifstream ReadInput;
	//Read input informations
	ReadInput.open("TSP_input.in");

	// numero città, se shape==2 o 3 lo sovrascrive con il numero delle città lette dal file esterno
	ReadInput >> N_cities;
	ReadInput >> shape;

	// se le posizioni le carico da file esterni devo ricontare il numero delle città
	if (shape == 2 || shape == 3 || shape == 4) {
		CityCounter(shape, N_cities);
	}
			
	if (N_cities <=3) { // altrimenti la soluzione è la somma delle distanze
		cerr << "Error: the number of cities must be greater than 3" << endl;
		MPI_Finalize();
		exit(0); // -1
	}
	
	if (shape == 0 && rank == 0) {
		cout << "Randomly placing " << N_cities << " cities inside a square..." << endl;
	} else if (shape == 1 && rank == 0) {
		cout << "Randomly placing " << N_cities << " cities on a circumference..." << endl;
	} else if (shape == 2 && rank == 0) {
		cout << "Placing cities on the given positions..." << endl;
	} else if (shape == 3 && rank == 0) {
		cout << "Using American capitals..." << endl;
	} else if (shape == 4 && rank == 0) {
		cout << "Using Italian towns..." << endl;
	}
	
	if (shape != 0 && shape != 1 && shape != 2 && shape != 3 && shape != 4) {
		cerr << "Error: shape must be 0, 1, 2, 3 or 4 in Input; exit" << endl;
		MPI_Finalize();
		exit(0);
	}

	ReadInput >> dim_pop;
	if (dim_pop <= 1) { // altrimenti non posso fare crossover
		cerr << "Error: the population must be greater than 1" << endl;
		MPI_Finalize();
		exit(0); // -1
	}
	if (rank == 0) {
	cout << "Creating a population of " << dim_pop << " elements..." << endl;
	}
	
	ReadInput >> N_gen;
	if (N_gen <= 0) { // altrimenti non ha senso
		cerr << "Error: the number of generations must be greater than 0" << endl;
		MPI_Finalize();
		exit(0); // -1
}
	if (rank == 0) {
	cout << "Creating " << N_gen << " generations..." << endl << endl;
	}

	ReadInput >> N_migr;
	if (N_migr < 1 || N_migr > N_gen) {
		cerr << "Error: the number of generations after migration must be greater than 1 and lower than the number of generations." << endl;
		MPI_Finalize();
		exit(0); // -1
}
	if (rank == 0) {
		cout << "Sending best individuals every " << N_migr << " generations." << endl;
		ReadInput.close();
	}
} // Fine Input()

//===========================================================
//Funzione che restituisce un vettore di struct con il numero della città e loro posizione, specificando la forma da utilizzare per la generazione casuale
// OSS: una volta chiamata questa funzione, le posizioni restano fissate ma cambiano solo gli indici
vector<Posizione> CityPlacer(int N_cty, int sh) {
	vector<Posizione> positions;

	for (int n = 1; n <= N_cty; n++) {
		Posizione position;
		position.index = n;
		
		if (sh == 0) { // Inside a Square
			position.x = rnd.Rannyu();
			position.y = rnd.Rannyu();
		} else if (sh == 1) { // On a circumference
			double theta = rnd.Rannyu(0, 2.*M_PI);
			position.x = 0.5*cos(theta);
			position.y = 0.5*sin(theta);
		} else {
			cerr << "Error: shape must be 1 or 0 in CityPlacer; exit" << endl;
			MPI_Finalize();
			exit(0);
		}
		positions.push_back(position);
	}

	// Aggiungo la prima posizione anche all'ultimo posto (stessa città di partenza e arrivo)
	Posizione first_pos = positions[0];
	positions.push_back(first_pos);
	
	return positions;
}

void PrintPositionsFile(vector<Posizione> positions, int sh) {
	// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Positions.dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Positions.dat";
	} else if (sh == 2) {
		file_path = "./data/InputShape/Positions.dat";
	} else if (sh == 3) {
		file_path = "./data/AmericanCities/Positions.dat";
	} else if (sh == 4) {
		file_path = "./data/ItalianCities/Positions.dat";
	} else {
		cerr << "Error: shape must be 0, 1, 2, 3 or 4 in PrintPositionsFile; exit" << endl;
		MPI_Finalize();
		exit(0);
	}
	
	const int wd = 15;
	// numero città
	int N = positions.size()-1;

	// Creo lo stream per salvare i dati
    ofstream outFile;
	outFile.open(file_path);
    if (outFile.is_open()) {
		// stampo la posizione di ogni città
		outFile << "index" << setw(wd) << "x" << setw(wd) << "y" << endl;
		for (int i=0; i<N; i++) {
        	outFile <<  positions[i].index <<  setw(wd) << fixed << setprecision(5) << positions[i].x  << setw(wd) << positions[i].y << endl;
		}

        outFile.close();
    } else {
        cerr << "Impossibile aprire il file al percorso " << file_path << endl;
    }
}

//===========================================================
//Funzione che crea la popolazione iniziale di dimPop elementi, a partire dalla configurazione iniziale delle città e anche la distanza associata ad ogni configurazione
vector<Individual> GeneratePopulation(int dimPop, vector<Posizione> startPos) {

	int length = startPos.size();
	vector<Individual> Population(dimPop);
	// aggiungo la prima sequenza con la sua distanza alla popolazione
	Population[0].pos = startPos;
	Population[0].dist = L_1(startPos);
//	Population[0].parent = "M0";
	// genero dimPop popolazioni, facendo uno shuffle della sequenza iniziale (partendo da i=1 perché la prima l'ho già settata)
	for (int i=1; i<dimPop; i++) {
		// copio il set di posizioni iniziali
		vector<Posizione> positions = startPos;

		// non uso la funzione integrata shuffle per avere un controllo diretto sul seed
		// scambio ogni città dalla 2 alla penultima con un'altra città ottenuta con un indice casuale
		for (int j=1; j<length-1; j++) {
			// genero un indice casuale (osserva che arriva fino a dimPop-2, perché l'estremo dx del Rannyu non è incluso)
			int randIndex = int(rnd.Rannyu(1, length-1));
			// scambio i suoi elementi
			swap(positions[j], positions[randIndex]);
		}
		// salvo la nuova successione
		Population[i].pos = positions;
		Population[i].dist = L_1(positions);
	//	Population[i].parent = "M0";
	}
	// controllo BC
	CheckPopBoundOk(Population);
	return Population;
}

//===========================================================
// Funzione che calcola la distanza tra due città
// Se carico le città americane o italiane, queste hanno le coordinate in longitudine e latitudine, quindi la distanza la calcolo diversamente
float Distance(Posizione p1, Posizione p2, int sh) {
	if (shape < 3) {
		return sqrt(pow((p2.x-p1.x), 2) + pow((p2.y-p1.y), 2));
	} else { // Calcolo la distanza tra due città in chilometri
		// Raggio medio della Terra in chilometri
		const double earthRadiusKm = 6371.0;

		// Converte la latitudine e la longitudine in radianti
		double lat1 = p1.y * M_PI / 180.0; // latitudine
		double lon1 = p1.x * M_PI / 180.0;
		double lat2 = p2.y * M_PI / 180.0;
		double lon2 = p2.x * M_PI / 180.0;

		// Calcola la differenza tra le longitudini
		double dlon = lon2 - lon1;
		// Uso formula distanza sulla superficie della Terra
		double distance = earthRadiusKm * acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon));

		return distance;
	}
}

//===========================================================
//Funzione che calcola la distanza totale data una certa configurazione (funzione costo)
float L_1(vector<Posizione> positions) {
	int N_dist = positions.size();
	float TotDist = 0.;

	for (int n=1; n<N_dist; n++) {
		float distance = Distance(positions[n], positions[n-1], shape);
		// cout << distance << endl;
		TotDist += distance;
	}

	return TotDist;
}

//===========================================================
//Funzione che controlla le condizioni al contorno (ogni città si ripete 1! volta, con la prima che compare sempre al primo e ultimo posto) e ferma il codice se non vengono rispettate
void CheckBoundOk(vector<Posizione> positions) {
	bool is_ok = true; // creo una variabile boolena invece di usare return per analizzare tutti gli indici senza uscire in anticipo
	int N = positions.size();
	if (positions[0].index != 1) {
		cerr << "Invalid BC: city nr.1 is not at the first position (" << positions[0].index << ")" << endl;
		is_ok = false;
	} 

	for (int i=1; i<N-1; i++) {			
		for (int j=i+1; j<N; j++) {
			if (positions[i].index == positions[j].index) {
				cerr << "Invalid BC: repeated city (" << positions[i].index << ")" << endl;
				is_ok = false; // Gli indici non sono tutti distinti
			}
		}
	}

	if (positions[N-1].index != 1) {
		cerr << "Invalid BC: city nr.1 is not at the last position (" << positions[N-1].index << ")" << endl;
		is_ok = false;
	}

	if (is_ok == false) { 
		MPI_Finalize();
		exit(0); // -1
	} ;
}

//========================================================================
//Funzione che controlla le BC di una popolazione
void CheckPopBoundOk(vector<Individual> pop) {
	int dim = pop.size();
	for (int i=0; i< dim; i++) {
		CheckBoundOk(pop[i].pos);
	}
}

//========================================================================
//Funzione che seleziona una configurazione casuale tra le popolazioni, con probabilità proporzionale all'inverso della distanza
// Utilizzo l'algortimo della roulette wheel selection
int Selector(vector<Individual> population) {
	int dimPop = population.size();
	// Calcola la somma dell'inverso delle distanze (per a normalizzazione)
    double inverSum = 0.0;
    for (const Individual &config : population) {
        inverSum += 1. / config.dist;
    }

    double r = rnd.Rannyu(); // numero casuale per la selezione
    double cumul_Sum = 0.; // somma cumulativa delle probabilità
    int sel_Index = 0; // indice selezionato
	
    for (int i = 0; i < dimPop; i++) {
		// se il numero casuale è troppo grande continuo a incrementare la soglia di selezione
        cumul_Sum += (1. / population[i].dist) / inverSum;
		
        if (r <= cumul_Sum) {
            sel_Index = i;
            break;
        }
    }
	// ovvero, le configurazioni con distanza bassa corrispondono a un grande incremento della probabilità cumulata, quindi è più probabile raggiungere il livello selezionato da r e quindi che vengano estratte
	
	// estraggo la configurazione selezionata
	return sel_Index;
}

// creare una sigola funzione mutazione?
//===========================================================
//Funzione mutazione 1: Pair Permutation, seleziona casualmente due città consecutive e le scambia
void Mutazione1(Individual & sequence, int rank) {
	// numero città
	int N = sequence.pos.size()-1;
	// seleziono un indice casuale, dalla seconda alla terzultima città
	int m = int(rnd.Rannyu(1, N-1));
	// scambio la coppia consecutiva
	swap(sequence.pos[m], sequence.pos[m+1]);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		cout << "mutation 1, index: " << m << endl;
		PrintPos(sequence.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << sequence.dist << endl;
		cout << "--------" << endl;
	}
}

//===========================================================
//Funzione mutazione 2: shifto i primi m elementi (dalla seconda città) di n posizioni a destra
void Mutazione2(Individual & sequence, int rank) {
	// numero città
	int N = sequence.pos.size()-1;
	// l'indice che indentifica il numero di elementi da spostare (m) va da 2 a N-1, ovvero posso shiftare minimo le prime 2 città oltre la prima e massimo le prime N-1 (più la prima fanno N città)
	int m = int( rnd.Rannyu(2, N) );
	// l'indice casuale di shift (n) va da 1 a m-1 (altrimeni non rimango nella sottosequenza dei primi m elementi o resterebbe invariata)
	int n = int( rnd.Rannyu(1, m) );
	// rotate prende come primo e terzo ingresso l'inizio e la fine (meno 1) del vector a cui applicare la rotazione (indice finale escluso, quindi lo incremento di 1)
	// il secondo ingresso rappresenta il valore che diventerà il primo posto dopo la rotazione (questo è incluso, ma volgio n_min=1 quindi lo incremento di 1)
    rotate(sequence.pos.begin() + 1, sequence.pos.begin() + m + 1 - n, sequence.pos.begin() + m + 1);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		cout << "mutation 2, shift and move index: " << n << " " << m << endl;
		PrintPos(sequence.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << sequence.dist << endl;
		cout << "--------" << endl;
	}
}

//===========================================================
//Funzione mutazione 3: Swap del primo blocco di m posizioni con un secondo blocco delle successive m (m<N_cities/2), a partire dalla seconda città
void Mutazione3(Individual & sequence, int rank) {
	// numero città
	int N = sequence.pos.size()-1;
	// genero l'indice casuale da 1 a N/2 escluso
	// in questo modo se fosse N=7 -> N/2=3.5 viene arrotondato a 4 e quindi l'intero random arriva fino a 3
	// se N=6 -> N/2=3 resta 3 e il random arriva fino a 2
	int m = int( rnd.Rannyu(1, static_cast<int>( round(N/2.) )) );
	// swappo il primo blocco di m elementi (dalla seconda città) con il successivo blocco di m elementi
	for (int i=0; i<m; i++) {
		swap(sequence.pos[1+i], sequence.pos[1+i+m]);
	}
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		cout << "mutation 3, index: " << m << endl;
		PrintPos(sequence.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << sequence.dist << endl;
		cout << "--------" << endl;
	}
}

//===========================================================
//Funzione mutazione 4, inversione della sequenza per i primi m valori oltre la prima (2<=m<=N-2 altrimenti rimane invariata)
void Mutazione4(Individual & sequence, int rank) {
	// numero città
	int N = sequence.pos.size()-1;
	// genero l'indice casuale da 2 a N-2
	int m = int(rnd.Rannyu(2, N-1));
	// inverto la sequenza dalla seconda città all'm-esima (l'ultimo indice non è preso)
	reverse(sequence.pos.begin()+1, sequence.pos.begin()+m+1);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		cout << "mutation 4, index: " << m << endl;
		PrintPos(sequence.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << sequence.dist << endl;
		cout << "--------" << endl;
	}
}

//===========================================================
//Crossover: prende una sequenza e la spezza ad una certa lunghezza casuale. Estrae la prima sequenza dalle popolazioni (per questo la popolazione dev'essere ordinata, altrimenti accoppierei il primo genitore con uno a caso), la spezza alla stessa lunghezza e prova a unire l'inizio della prima con la fine della seconda. Continua ad estrarre una sequenza (diversa dalla sequenza fornita) finché l'unione non rispetta le BC.
//Modifica per riferimento anche l'indice che identifica il secondo genitore nella popolazione precedente. Se provo tutte le configurazioni della popolazione e non riesce a fare crossover lo modifica in -1
Individual Crossover(Individual & sequence, vector<Individual> pop, int & index, int rank) {
	// copia del genitore 1, lo salvo così se i figli sono più lunghi ho modo di restituire i genitori
	Individual parent_1 = sequence;
	// candidato genitore 2, se non riesce il crossover equivale all'ultimo elemento della popolazione precedente
	// lo inizializzo qui perché altrimenti non posso accedervi fuori dal ciclo for
	Individual selected_seq;

	// creo la variabile che ospiterà una copia del secondo genitore
	Individual parent_2;
	
	// indice della sequenza con la quale avverrà il crossover (se resta -1 non è avvenuto)
	index = -1;
	// dimensione popolazione
	int dimPop = pop.size();
	// numero città
	int N = sequence.pos.size()-1;
	// genero l'indice casuale da 2 a N-2 (Escludo 1, N-1 e N altrimenti i figli sarebbero identici), mi dice quante città prendere (contando la prima)
	int m = int(rnd.Rannyu(2, N-1));
	// se è riuscito a fare il crossover oppure no
	bool canCross = false; 
	
	// salvo gli indici della sequenza fornita (genitore 1)
	vector<int> seq1;
	for(Posizione pos: sequence.pos) { seq1.push_back(pos.index); }
	
	// continuo a selezionare una sequenza della popolazione finché non avviene il crossover
	// OSSERVAZIONE: il crossover può avvenire (rispettando le BC) sse i primi m valori, dopo il primo, delle due sequenze sono uguali tra loro
	for (int i=0; i<dimPop; i++) {
		// candidato genitore 2
		selected_seq = pop[i];
		// salvo gli inidici della sequenza selezionata
	    vector<int> seq2;
	    for(Posizione pos: selected_seq.pos) { seq2.push_back(pos.index); }

		// controllo se le due sequenze sono diverse (evito che due sequenze identiche diano una figlia identica; anche nel caso siano una l'opposto dell'altra sono uguali, ma non possono accoppiarsi perché hanno i primi valori sicuramente diversi)
		if (seq1 != seq2) { // se invece sono uguali ne seleziona un'altra
			// controllo se le prime m sequenze sono uguali (anche non nello stesso ordine)
		    canCross = std::is_permutation(seq1.begin(), seq1.begin()+m, seq2.begin());
			// se posso fare il crossover
			if (canCross == true) {
				// salvo una copia del genitore 2
				parent_2 = selected_seq;
				// salvo una copia del genitore 1
				Individual sequence_copy = sequence;
				// eseguo il crossover: scambio la seconda metà del gen 2 con quella del gen 1
	        	copy(selected_seq.pos.begin() + m, selected_seq.pos.end(), sequence.pos.begin() + m);
				// scambio la seconda metà del gen 1 con quella del gen 2
	        	copy(sequence_copy.pos.begin() + m, sequence_copy.pos.end(), selected_seq.pos.begin() + m);
				// salvo l'indice del gen 2
				index = i;
				break;
			}
		}
	}
	// calcolo la nuova distanza per entrambi i figli
	// OSS: se non riesce a fare crossover qui sto settando il secondo genitore come l'ultimo della popolazione, ma tanto l'indice di controllo è negativo quindi me ne accorgo e non lo aggiungo alla nuova generazione
	sequence.dist = L_1(sequence.pos);
	selected_seq.dist = L_1(selected_seq.pos);

	if (globalPrintMutation == true) {
		cout << "==== RANK " << rank << " ====" << endl;
		if (canCross) { cout << "crossover, cut after: " << m << endl; }
		else { cout << "fail crossover, cut after: " << m << endl; }
		cout << "child 1:" << endl;
		PrintPos(sequence.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << sequence.dist << endl;
		cout << "--------" << endl;
		cout << "child 2, from index: " << index << endl;
		PrintPos(selected_seq.pos);
		cout << "--------" << endl;
		cout << "dist:" << endl;
		cout << selected_seq.dist << endl;
		cout << "--------" << endl;
	}

	return selected_seq;
}

//==========================================
//Funzione che crea una nuova generazione a partire da una vecchia, utilizzando le mutazioni e il crossover
void NewGeneration(vector<Individual> & new_pop, vector<Individual> pop, int rank) {
	int dimPop = pop.size();
	int i = 0; // contatore elementi nella popolazione
	while (i < dimPop) {
		// seleziono una configurazione casuale della generazione precedente
		int sel_index = Selector(pop);
		// creo una copia  e poi la passo per riferimento alle mutazioni, così la popolazione precedente non viene modificata e posso pescare ancora le stesse configurazioni
		Individual sel_config = pop[sel_index];
		// creo un'altra copia per poter tornare al genitore se il figlio dovesse avere una distanza maggiore
		Individual sel_config_c = sel_config;
		
		if (globalPrint == true) {
			cout << "==== RANK " << rank << " ====" << endl;
			cout << "Individual " << i << " , selected sequence: " << sel_index + 1 << endl;
			cout << "-------" << endl;
		}

		// creo l'etichetta che identifica il genitore dalla quale ha mutato
	//	string parent_code = "M"+to_string(sel_index);
		
		// chiamo le mutazioni e il crossover con certe probabilità
		if (rnd.Rannyu() < p_mut) {
			Mutazione1(sel_config, rank);
			// se la mutazione non ha ridotto la distanza torno al genitore con una certa prob di sopravvivenza
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "==== RANK " << rank << " ====" << endl;
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} // se invece la distanza viene ridotta la salvo
			else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione2(sel_config, rank);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "==== RANK " << rank << " ====" << endl;
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione3(sel_config, rank);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "==== RANK " << rank << " ====" << endl;
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione4(sel_config, rank);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "==== RANK " << rank << " ====" << endl;
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_cross) {
			// restituisce il secondo figlio Individual da crossover con il sel_config; se non riesce a fare crossover l'indice è negativo e la sequenza resta invariata. Cambia per riferimento l'indice della sequenza con cui l'ha fatto
			int partner_index = 0; // indice del secondo genitore
			Individual child_2_config = Crossover(sel_config, pop, partner_index, rank);
			// solo se il crossover è avvenuto aggiungo il secondo figlio alla generazione
			
			if (partner_index >= 0) {
				// controllo le BC del partner				
				CheckBoundOk(child_2_config.pos);
				// salvo il secondo figlio nella nuova popolazione se è più corto del primo genitore (con una certa probabilità)
				if (child_2_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
					child_2_config = sel_config_c;
					if (globalPrint == true) {
						cout << "==== RANK " << rank << " ====" << endl;
						cout << "Distance unreduced, returning the parent" << endl;
					}
				}
				
				new_pop[i] = child_2_config;
				// incremento l'indice per indicare che ho aggiunto un individual in più alla popolazione, così fuori da questo if aggiungo ancora un altro indice per il sel_config
				i+=1;
			} //else { parent_code = "E"+to_string(sel_index); }
		}
		// Check BC after mutation
		CheckBoundOk(sel_config.pos);
		
		// aggiungo l'indice del genitore
	//	sel_config.parent = parent_code;
		
		// prima di aggiungere controllo se posso farlo (magari l'ultima estrazione è un crossover quindi mi aggiungerebbe due sequenze, ma per quella in più non c'è posto)
		if (i<dimPop) { new_pop[i] = sel_config; }
		// incremento l'indice della popolazione
		i+=1;
	}
}

//==========================================
// Funzione che stampa le distanze di ogni generazione (appendendo i valori)
void PrintDistancesFile(vector<vector<Individual>> generations, int sh, int rank) {
	// cambio percorso in base alla shape e al rank
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Distances_" + to_string(rank) + ".dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Distances_" + to_string(rank) + ".dat";
	} else if (sh == 2) {
		file_path = "./data/InputShape/Distances_" + to_string(rank) + ".dat";
	} else if (sh == 3) {
		file_path = "./data/AmericanCities/Distances_" + to_string(rank) + ".dat";
	} else if (sh == 4) {
		file_path = "./data/ItalianCities/Distances_" + to_string(rank) + ".dat";
	} else {
		cerr << "Error: shape must be 0, 1, 2, 3 or 4 in PrintDistanceFile; exit" << endl;
		MPI_Finalize();
		exit(0);
	}
		
	const int wd = 15;
	int dimGen = generations.size();	// numero generazioni
	int dimPop = generations[0].size();	// numero individui per popolazione

	// Creo lo stream per salvare i dati
    ofstream outFile;
	outFile.open(file_path);
    if (outFile.is_open()) {		
		// aggiungo l'intestazione
		for (int i=1; i<dimPop; i++) {
			outFile << "D_" << i << setw(wd);
		}
		outFile << "D_" << dimPop << endl;

		// Per ogni generazione:
		for (int gen=0; gen<dimGen; gen++) {
			// stampo la distanza di ogni sequenza della sua popolazione
			for (int i=0; i<=dimPop-2; i++) {
	        	outFile << fixed << setprecision(5) << generations[gen][i].dist  << setw(wd);
			}
			// per l'ultima valore non uso la spaziatura
	        outFile << fixed << setprecision(5) << generations[gen][dimPop-1].dist  << endl;
		}
		// chiudo lo stream
        outFile.close();
    } else {
        cerr << "Impossibile aprire il file al percorso " << file_path << endl;
    }
}

//==========================================
// NON PIU' USATA
// Funzione che stampa le sequenze di ogni generazione (appendendo i valori) e per ciascuna l'indice della sequenza da cui è avvenuta la mutazione o la coppia di indici da cui è avvenuto il crossover
void PrintSequencesFile(vector<Individual> pop, int sh, int rank) {
	// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Sequences_" + to_string(rank) + ".dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Sequences_" + to_string(rank) + ".dat";
	} else if (sh == 2) {
		file_path = "./data/InputShape/Sequences_" + to_string(rank) + ".dat";
	} else if (sh == 3) {
		file_path = "./data/AmericanCities/Sequences_" + to_string(rank) + ".dat";
	} else if (sh == 4) {
			file_path = "./data/ItalianCities/Sequences_" + to_string(rank) + ".dat";
	} else {
		cerr << "Error: shape must be 0, 1, 2, 3 or 4 in PrintSequencesFile; exit" << endl;
		MPI_Finalize();
		exit(0);
	}
	
	const int wd = 10;
//	const int wd2 = 3;
	int dimPop = pop.size();
	int dimSeq = pop[0].pos.size();

	// Creo lo stream per salvare i dati
    ofstream outFile;
	
	// verifico l'esistenza del file
    ifstream isFile(file_path);
	// se il file non esiste lo creo e aggiungo l'intestazione
    if (!isFile.is_open()) {		
		// creo il file
		outFile.open(file_path);
		// aggiungo l'intestazione
		for (int i=1; i<dimPop; i++) {
			outFile << "Par_" << i << setw(wd) << "Seq_" << i << setw(wd);
		}
		outFile << "Par_" << dimPop << setw(wd) << "Seq_" << dimPop << endl;
		outFile.close();
    }
	// chiudo lo stream di controllo
	isFile.close();

	// apro il file esistente in modalità append
	outFile.open(file_path, ios::app);
    if (outFile.is_open()) {
		// stampo ogni sequenza della popolazione
		for (int i=0; i<=dimPop-2; i++) {
		//	outFile << pop[i].parent << setw(wd2);
			// stampo gli indici di ogni sequenza
			for (int j=0; j<=dimSeq-2; j++) {
				outFile << pop[i].pos[j].index << "-";
			}
			// dopo l'ultimo indice non metto il trattino
			outFile << pop[i].pos[dimSeq-1].index << setw(wd);
		}
		
	//	outFile << pop[dimPop-1].parent << setw(wd2);
		// stampo l'ultimo senza il setw finale
		for (int j=0; j<=dimSeq-2; j++) {
			outFile << pop[dimPop-1].pos[j].index << "-";
		}
		// dopo l'ultimo indice non metto il trattino
		outFile << pop[dimPop-1].pos[dimSeq-1].index << endl;
		
        outFile.close();
    } else {
        cerr << "Impossibile aprire il file al percorso " << file_path << endl;
    }
}

//===========================================
// Funzione per trovare il minimo della distanza in una popolazione
Individual MinimumDistance(vector<Individual> pop) {

    auto minElement = min_element(pop.begin(), pop.end(), [](const Individual& seq1, const Individual& seq2) {
        return seq1.dist < seq2.dist;
    });

	// Estraggo l'intera struct
    Individual minIndividual = *minElement;

	return minIndividual;
}

//============================================
// Stampa un file con tutte le sequenze minime per ogni generazione
void PrintMinSequence(vector<Individual> seq, int sh, int rank) {
// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/minSequences_" + to_string(rank) + ".dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/minSequences_" + to_string(rank) + ".dat";
	} else if (sh == 2) {
		file_path = "./data/InputShape/minSequences_" + to_string(rank) + ".dat";
	} else if (sh == 3) {
		file_path = "./data/AmericanCities/minSequences_" + to_string(rank) + ".dat";
	} else if (sh == 4) {
			file_path = "./data/ItalianCities/minSequences_" + to_string(rank) + ".dat";
	} else {
		cerr << "Error: shape must be 0, 1, 2, 3 or 4 in PrintMinSequence; exit" << endl;
		MPI_Finalize();
		exit(0);
	}
	
	const int wd = 3;
	int dimGen = seq.size();
	int dimSeq = seq[0].pos.size();

	// Creo lo stream per salvare i dati
    ofstream outFile;
	outFile.open(file_path);
    if (outFile.is_open()) {
		// stampo tutte le distanze
		for (int gen=0; gen<dimGen; gen++) {
			outFile << seq[gen].dist << setw(wd);
			for (int i=0; i<dimSeq-1; i++) {
				outFile << seq[gen].pos[i].index << "-";
			}
			// dopo l'ultimo indice non metto il trattino
			outFile << seq[gen].pos[dimSeq-1].index << endl;
		}
		outFile.close();
	} else {
		cerr << "Impossibile aprire il file al percorso " << file_path << endl;
	}
}

//=====================================================
// funzione che confronta due sequenze in base alla distanza, in modo da riordinare una popolazione in ordine crescente
bool compareByDist(const Individual & seq1, const Individual & seq2) {
    return seq1.dist < seq2.dist;
};

//=====================================================
// Funzione che genera casualmente le coppie di rank che devono scambiarsi i best individual

vector<pair<int,int>> GeneraCoppie(int size) {
	vector<pair<int, int>> coppie;
     // Creo un vettore di rank da 0 a size-1, poi lo mischio per ottenere coppie di rank senza ripetizioni
	vector<int> ranks(size);
	for (int i = 0; i < size; i++) {
		ranks[i] = i; // es: 0,1,2,3
	}

	// scambio ogni rank in modo casuale (è come shuffle ma non uso la funzione integrata per avere controllo sul seed)
	for (int j=0; j<size; j++) {
		// genero un indice casuale (osserva che arriva fino a size-1, perché l'estremo dx del Rannyu non è incluso)
		int randIndex = int(rnd.Rannyu(0, size));
		// scambio i suoi elementi
		swap(ranks[j], ranks[randIndex]);
	}
	// es. ho ottenuto 2,0,1,3
    // Accoppio un rank con il suo successivo
    for (int i = 0; i < size - 1; i += 2) {
        coppie.push_back(std::make_pair(ranks[i], ranks[i	+1]));
	}
	// saltando di due in due se size è dispari l'ultimo rank non lo accoppio
	
	return coppie;
}

//===========================================================
//Funzione che carica da file esterno le posizioni delle città in un vettore di struct, con il loro numero e posizione. Carica anche le posizioni delle città americane se sh=3
// OSS: lo faccio eseguire solo dal rank 0 per poi passare i dati
vector<Posizione> CityLoader(int sh) {
	string file_path;
	
	if (sh == 2) {
		file_path = "./Positions.dat";
	} else if (sh == 3) {
		file_path = "./American_capitals.dat";
	} else if (sh == 4) {
		file_path = "./cap_prov_ita.dat";
	} else {
		cerr << "Error: shape must be 2, 3 or 4 in CityLoader; exit" << endl;
		MPI_Finalize();
		exit(0);
	}

	vector<Posizione> positions;

	ifstream file(file_path);

	if (!file.is_open()) {
		cerr << "Errore nell'apertura del file delle città." << endl;
		return positions;
	}

	string line;

	// dato che le città italiane non hanno l'intestazione nel file devo saltare la prima riga solo nel caso in cui sh=2 o 3
	if (sh==2 || sh==3) {
		// Ignoro la prima riga (intestazione)
		getline(file, line);
	}

	// mi serve solo a scartare i nomi delle città americane
	string name;
	string description;
	int idx = 1; // indice per contare le città

	while (getline(file, line)) {
		Posizione city_pos;
		istringstream iss(line);

		if (sh == 2) {
			iss >> city_pos.index >> city_pos.x >> city_pos.y;
		} else if (sh == 3) { // carico le città americane
			iss >> name >> description >> city_pos.x >> city_pos.y;
			city_pos.index = idx;
			idx++;
		} else { // le città italiane nel file non hanno nome e stato
			iss >> city_pos.x >> city_pos.y;
			city_pos.index = idx;
			idx++;
		}
		
		positions.push_back(city_pos);
	}

	// Aggiungo la prima posizione anche all'ultimo posto (stessa città di partenza e arrivo)
	Posizione first_pos = positions[0];
	positions.push_back(first_pos);

	return positions;
}

//==============================================
// funzione che conta il numero di città nel file
void CityCounter(int sh, int & N_cty) {
	N_cty = 0; // azzero il numero di città della variabile globale
	string file_path;

	if (sh == 2) {
		file_path = "./Positions.dat";
	} else if (sh == 3) {
		file_path = "./American_capitals.dat";
	} else if (sh == 4) {
			file_path = "./cap_prov_ita.dat"; 
	} else {
		cerr << "Error: shape must be 2, 3 or 4 in CityCounter; exit" << endl;
		MPI_Finalize();
		exit(0);
	}

	ifstream file(file_path);

	if (!file.is_open()) {
		cerr << "Errore nell'apertura del file delle città." << endl;
		MPI_Finalize();
		exit(0);
	}

	string line;
	// dato che le città italiane non hanno l'intestazione nel file devo saltare la prima riga solo nel caso in cui sh=2 o 3
	if (sh==2 || sh==3) {
		// Ignoro la prima riga (intestazione)
		getline(file, line);
	}

	// conto le città contando le righe nel file
	while (getline(file, line)) {
		N_cty++;  // incremento il numero di città
	}

	// controllo che il numero di città sia maggiore di 3 (altrimenti la distanza minima è la somma delle distanze)
	if (N_cty < 3) {
		cerr << "Error: Il numero di città è minore di 3." << endl;
		MPI_Finalize();
		exit(0);
	}
}
//==============================================
