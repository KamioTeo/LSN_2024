// Salesman Problem with Genetic Algorithm

#include "functions.h"

// opzioni per il debugging
bool globalPrint = false; // se stampare o no tutte le configurazioni / sequenze
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
	// leggo il numero di città da usare, il tipo di posizionamento e la dimensione della popolazione
	Input();

	// creo la configurazione iniziale delle città (indice e coordinate)
	vector<Posizione> x = CityPlacer(N_cities, shape);
	// controllo che le condizioni al contorno siano rispettate
	CheckBoundOk(x);
	// stampo le posizioni ottenute
	PrintPositionsFile(x, shape);

	// creo la popolazione, da permutazioni casuali della configurazione iniziale (Check BC svolto nella funzione)
	vector<Individual> population = GeneratePopulation(dim_pop, x);

	// Riordino la popolazione in ordine crescente di distanza (serve per il crossover)
    sort(population.begin(), population.end(), compareByDist);

	if (globalPrint == true) { PrintPop(population); }
	cout << "--------" << endl;
	
	// stampo la prima generazione e le sequenze
	PrintDistanceFile(population, shape);
	// PrintSequencesFile(population, shape);

	// salvo tutti gli Individual con distanza minima per ogni generazione
	vector<Individual> AllminSeq(N_gen);

	// della prima generazione
	AllminSeq[0] = MinimumDistance(population);
	if (globalPrint == true) {
		cout << "0 minimum distance: " << AllminSeq[0].dist << endl;
}

	for (int gen=1; gen<N_gen; gen++) {
		// creo un nuovo vector per ospitare la nuova generazione
		vector<Individual> new_population(dim_pop);
		// creo la nuova generazione a partire da quella vecchia (Check BC svolto nella funzione) usano le mutazioni e il crossover
		NewGeneration(new_population, population);
		// salvo la sequenza minima della nuova generazione
		AllminSeq[gen] = MinimumDistance(new_population);
		
		if (globalPrint == true) {
			// stampo la nuova popolazione
			cout << "-----Generation_" << gen+1 << "----" << endl;
			PrintPop(new_population);
			cout << "--------" << endl;
		}

		// riordino la nuova generazione
    	sort(new_population.begin(), new_population.end(), compareByDist);
		// stampo le distanze della nuova generazione
		PrintDistanceFile(new_population, shape);
		// stampo le sequenze (serve a disegnarle)
		// PrintSequencesFile(new_population, shape);

		// la nuova generazione diventa quella vecchia nel nuovo ciclo
		population = new_population;

		if (globalPrint == true) {
			cout << gen << " minimum distance: " << AllminSeq[gen].dist << endl;
		}
	}

// stampo la sequenza con la distanza minima
PrintMinSequence(AllminSeq, shape);
rnd.SaveSeed();

return 0;

}

//===========================================================
void Input(void) {
	ifstream ReadInput;
	
	cout << "_______________TSP_______________" << endl;
	cout << "The Salesman Problem code" << endl << endl;
	
	//Read input informations
	ReadInput.open("TSP_input.in");
	
	ReadInput >> N_cities;
	if (N_cities <=3) { // altrimenti la soluzione è la somma delle distanze
		cerr << "Error: the cities must be greater than 3" << endl;
		exit(-1);
	}
	
	ReadInput >> shape;
	if (shape == 0) {
		cout << "Randomly placing " << N_cities << " cities inside a square..." << endl;
	} else if (shape == 1) {
		cout << "Randomly placing " << N_cities << " cities on a circumference..." << endl;
	} else {
		cerr << "Error: shape must be 1 or 0; exit" << endl;
		exit(0);
	}

	ReadInput >> dim_pop;
	if (dim_pop <= 1) { // altrimenti non posso fare crossover
		cerr << "Error: the population must be greater than 1" << endl;
		exit(-1);
	}
	cout << "Creating a population of " << dim_pop << " elements..." << endl;

	ReadInput >> N_gen;
	if (N_gen <= 0) { // altrimenti non ha senso
		cerr << "Error: the number of generations must be greater than 0" << endl;
		exit(-1);
}
	cout << "Creating " << N_gen << " generations..." << endl << endl;
	
	ReadInput.close();
} // Fine Input()

//===========================================================
// Funzione che restituisce un vettore di struct con il numero della città e loro posizione, specificando la forma da utilizzare per la generazione casuale
// OSS: una volta chiamata questa funzione, le posizioni restano fissate ma cambiano solo gli indici
vector<Posizione> CityPlacer(int N_cty, bool sh) {
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
			cerr << "Error: shape must be 1 or 0 in function 'CityPlacer'; exit" << endl;
			exit(0);
		}
		positions.push_back(position);
	}

	// Aggiungo la prima posizione anche all'ultimo posto (stessa città di partenza e arrivo)
	Posizione first_pos = positions[0];
	positions.push_back(first_pos);
	
	return positions;
}

void PrintPositionsFile(vector<Posizione> positions, bool sh) {
	// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Positions.dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Positions.dat";
	} else {
		cerr << "Error: shape must be 1 or 0; exit" << endl;
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
//Funzione che calcola la distanza tra due città
float Distance(Posizione p1, Posizione p2) {
	return sqrt(pow((p2.x-p1.x), 2) + pow((p2.y-p1.y), 2));
}

//===========================================================
//Funzione che calcola la distanza totale data una certa configurazione (funzione costo)
float L_1(vector<Posizione> positions) {
	int N_dist = positions.size();
	float TotDist = 0.;

	for (int n=1; n<N_dist; n++) {
		float distance = Distance(positions[n], positions[n-1]);
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

	if (is_ok == false) { exit(-1); } ;
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
// Utilizzo l'algoritmo della roulette wheel selection
int Selector(vector<Individual> population) {
	int dimPop = population.size();
	// Conterrà la somma dell'inverso delle distanze (serve per la normalizzazione)
    double inverSum = 0.0;
    for (const Individual &config : population) {
        inverSum += 1. / config.dist;
    }

    double r = rnd.Rannyu();	// numero casuale per la selezione
    double cumul_Sum = 0.; 		// somma cumulativa delle probabilità
    int sel_Index = 0; 			// indice selezionato
	
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

//===========================================================
// Funzione mutazione 1: Pair Permutation, seleziona casualmente due città consecutive e le scambia
void Mutazione1(Individual & sequence) {
	// numero città
	int N = sequence.pos.size()-1;
	// seleziono un indice casuale, dalla seconda alla terzultima città (non posso prendere la penultima altrimenti la scambierei conl'ultima, che però deve sempre essere la prima)
	int m = int(rnd.Rannyu(1, N-1));
	// scambio la coppia consecutiva
	swap(sequence.pos[m], sequence.pos[m+1]);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
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
void Mutazione2(Individual & sequence) {
	// numero città
	int N = sequence.pos.size()-1;
	// l'indice che indentifica il numero di elementi da spostare (m) va da 2 a N-1, ovvero posso shiftare minimo le prime 2 città oltre la prima e massimo le prime N-1 (infatti inclusa la prima fanno N città)
	int m = int( rnd.Rannyu(2, N) );
	// l'indice casuale di shift (n) va da 1 a m-1 (altrimeni non rimango nella sottosequenza dei primi m elementi o resterebbe invariata)
	int n = int( rnd.Rannyu(1, m) );
	// "rotate" prende come primo e terzo ingresso l'inizio e la fine (meno 1) del vector a cui applicare la rotazione (indice finale escluso, quindi lo incremento di 1)
	// il secondo ingresso rappresenta il valore che diventerà il primo posto dopo la rotazione (questo è incluso, ma voglio n_min=1 quindi lo incremento di 1)
    rotate(sequence.pos.begin() + 1, sequence.pos.begin() + m + 1 - n, sequence.pos.begin() + m + 1);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
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
void Mutazione3(Individual & sequence) {
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
void Mutazione4(Individual & sequence) {
	// numero città
	int N = sequence.pos.size()-1;
	// genero l'indice casuale da 2 a N-2
	int m = int(rnd.Rannyu(2, N-1));
	// inverto la sequenza dalla seconda città all'm-esima (l'ultimo indice non è preso)
	reverse(sequence.pos.begin()+1, sequence.pos.begin()+m+1);
	// calcolo la nuova distanza
	sequence.dist = L_1(sequence.pos);

	if (globalPrintMutation == true) {
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
// DEPRECATED: Modifica per riferimento anche l'indice che identifica il secondo genitore nella popolazione precedente. Se provo tutte le configurazioni della popolazione e non riesce a fare crossover lo modifica in -1
Individual Crossover(Individual & sequence, vector<Individual> pop, int & index) {
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
void NewGeneration(vector<Individual> & new_pop, vector<Individual> pop) {
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
			cout << "Individual " << i << " , selected sequence: " << sel_index + 1 << endl;
			cout << "-------" << endl;
		}

		// creo l'etichetta che identifica il genitore dalla quale ha mutato
	//	string parent_code = "M"+to_string(sel_index);
		
		// chiamo le mutazioni e il crossover con certe probabilità
		if (rnd.Rannyu() < p_mut) {
			Mutazione1(sel_config);
			// se la mutazione non ha ridotto la distanza torno al genitore con una certa prob di sopravvivenza
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} // se invece la distanza viene ridotta la salvo
			else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione2(sel_config);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione3(sel_config);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_mut) {
			Mutazione4(sel_config);
			if (sel_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
				sel_config = sel_config_c;
				if (globalPrint == true) {
					cout << "Distance unreduced, returning the parent" << endl;
				}
			} else { sel_config_c = sel_config; }
		};
		if (rnd.Rannyu() < p_cross) {
			// restituisce il secondo figlio Individual da crossover con il sel_config; se non riesce a fare crossover l'indice è negativo e la sequenza resta invariata. Cambia per riferimento l'indice della sequenza con cui l'ha fatto
			int partner_index = 0; // indice del secondo genitore
			Individual child_2_config = Crossover(sel_config, pop, partner_index);
			// solo se il crossover è avvenuto aggiungo il secondo figlio alla generazione
			
			if (partner_index >= 0) {
				// controllo le BC del partner				
				CheckBoundOk(child_2_config.pos);
				// salvo il secondo figlio nella nuova popolazione se è più corto del primo genitore (con una certa probabilità)
				if (child_2_config.dist > sel_config_c.dist && rnd.Rannyu() < p_surv) {
					child_2_config = sel_config_c;
					if (globalPrint == true) {
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
void PrintDistanceFile(vector<Individual> pop, bool sh) {
	// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Distances.dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Distances.dat";
	} else {
		cerr << "Error: shape must be 1 or 0; exit" << endl;
		exit(0);
	}
		
	const int wd = 15;
	int dimPop = pop.size();

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
			outFile << "D_" << i << setw(wd);
		}
		outFile << "D_" << dimPop << endl;
		outFile.close();
    }
	// chiudo lo stream di controllo
	isFile.close();

	// apro il file esistente in modalità append
	outFile.open(file_path, ios::app);
    if (outFile.is_open()) {
		// stampo la distanza di ogni sequenza della popolazione
		for (int i=0; i<=dimPop-2; i++) {
        	outFile << fixed << setprecision(5) << pop[i].dist  << setw(wd);
		}
        outFile << fixed << setprecision(5) << pop[dimPop-1].dist  << endl;
			
        outFile.close();
    } else {
        cerr << "Impossibile aprire il file al percorso " << file_path << endl;
    }
}

//==========================================
// Funzione che stampa le sequenze di ogni generazione (appendendo i valori) e per ciascuna l'indice della sequenza da cui è avvenuta la mutazione o la coppia di indici da cui è avvenuto il crossover
void PrintSequencesFile(vector<Individual> pop, bool sh) {
	// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/Sequences.dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/Sequences.dat";
	} else {
		cerr << "Error: shape must be 1 or 0; exit" << endl;
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
void PrintMinSequence(vector<Individual> seq, bool sh) {
// cambio percorso in base alla shape
	string file_path;
	if (sh == 0) {
		file_path = "./data/Square/minSequences.dat";
	} else if (sh == 1) {
		file_path = "./data/Circle/minSequences.dat";
	} else {
		cerr << "Error: shape must be 1 or 0; exit" << endl;
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