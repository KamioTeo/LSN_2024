#ifndef __functions__
#define __functions__

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
//Random numbers
#include "../RandomNumberGenerator/random.h"

using namespace std;

// Generazione del random number generator
Random rnd("../RandomNumberGenerator/Primes", "../RandomNumberGenerator/seed.in");

extern bool globalPrint; // Boolean che dice se stampare o no tutte le configurazioni
extern bool globalPrintMutation; // Boolean che dice se stampare o no tutte le configurazioni per ogni mutazione

int N_cities; 	// numero città da visitare
bool shape; 	// forma sulla quale sono disposte le città (vedi file input)
int dim_pop; 	// dimensione popolazione (configurazioni delle città)
int N_gen; 		// Numero di generazioni

// struct che identifica ogni città
struct Posizione {
	int index;
	float x;
	float y;
};

// struct che identifica una sequenza di città
struct Individual {
	vector<Posizione> pos;
	float dist;	// distanza di quella configurazione
//	string parent; // indice della configurazione genitore che l'ha generata
};

// DEPRECATED:
/* LEGENDA INDICI parents
la nuova configurazione è:
M{int} -> Mutazione della sequenza numero "int" della generazione precedente
E{int} -> Equivalente alla sequenza numero "int" della generazione precedente
c{int1}c{int2} -> crossover tra le sequenze "int1" e "int2" della generazione precedente
*/

void Input(void);
vector<Posizione> CityPlacer(int, bool); // crea una configurazione delle città
void PrintPositionsFile(vector<Posizione>, bool); // stampa il file con le coordinate delle varie città e il loro numero
vector<Individual> GeneratePopulation(int, vector<Posizione>); // crea la popolazione iniziale
float Distance(Posizione, Posizione); // distanza tra due città
float L_1(vector<Posizione>); // funzione costo
void CheckBoundOk(vector<Posizione>); //Check boundary conditions
void CheckPopBoundOk(vector<Individual> pop); //Check BC di una popolazione
int Selector(vector<Individual>); // seleziona una configurazione casuale (un indice nella popolazione), con probabilità proporzionale all'inverso della distanza

void Mutazione1(Individual &);
void Mutazione2(Individual &);
void Mutazione3(Individual &);
void Mutazione4(Individual &);

Individual Crossover(Individual &, vector<Individual>, int &);

void NewGeneration(vector<Individual> &, vector<Individual>);

void PrintDistanceFile(vector<Individual>, bool);
void PrintSequencesFile(vector<Individual>, bool);

Individual MinimumDistance(vector<Individual>); // estrae l'individual della popolazione per cui la distanza è minima

void PrintMinSequence(vector<Individual>, bool);
bool compareByDist(const Individual &, const Individual &);

#endif