#ifndef __functions__
#define __functions__

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

#include <sstream> // per la lettura di file

#include <mpi.h>
//Random numbers
#include "../RandomNumberGenerator/random.h"

using namespace std;

// Generazione del random number generator (sovrascritto nel main per averne uno diverso per ogni core)
Random rnd("../RandomNumberGenerator/primes32001.in", "../RandomNumberGenerator/seed.in");

extern bool globalPrint; // Boolean che dice se stampare o no tutte le configurazioni
extern bool globalPrintMutation; // Boolean che dice se stampare o no tutte le configurazioni per ogni mutazione

int N_cities; // numero città da visitare
int shape; // forma sulla quale sono disposte le città (vedi file input)
int dim_pop; // dimensione popolazione (configurazioni delle città)
int N_gen; // Numero di generazioni

int N_migr; // Numero di generazioni dopo la quale avviene lo swap tra core

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

// DEPRECATED
// ATTENZIONE: Attributo parents eliminato nelle versioni più recenti del codice
/* LEGENDA INDICI parents
la nuova configurazione è:
M{int} -> Mutazione della sequenza numero "int" della generazione precedente
E{int} -> Equivalente alla sequenza numero "int" della generazione precedente
c{int1}c{int2} -> crossover tra le sequenze "int1" e "int2" della generazione precedente
*/

void Input(int);
vector<Posizione> CityPlacer(int, int); // crea una configurazione delle città
void PrintPositionsFile(vector<Posizione>, int); // stampa il file con le coordinate delle varie città e il loro numero
vector<Individual> GeneratePopulation(int, vector<Posizione>); // crea la popolazione iniziale
float Distance(Posizione, Posizione, int); // distanza tra due città
float L_1(vector<Posizione>); // funzione costo
void CheckBoundOk(vector<Posizione>); //Check boundary conditions
void CheckPopBoundOk(vector<Individual> pop); //Check BC di una popolazione
int Selector(vector<Individual>); // seleziona una configurazione casuale (un indice nella popolazione), con probabilità proporzionale all'inverso della distanza

void Mutazione1(Individual &, int);
void Mutazione2(Individual &, int);
void Mutazione3(Individual &, int);
void Mutazione4(Individual &, int);

Individual Crossover(Individual &, vector<Individual>, int &, int);

void NewGeneration(vector<Individual> &, vector<Individual>, int);

void PrintDistancesFile(vector<vector<Individual>>, int, int);
void PrintSequencesFile(vector<Individual>, int, int);

Individual MinimumDistance(vector<Individual>); // estrae l'individual della popolazione per cui la distanza è minima

void PrintMinSequence(vector<Individual>, int, int);
bool compareByDist(const Individual &, const Individual &);

vector<pair<int,int>> GeneraCoppie(int);
vector<Posizione> CityLoader(int);
void CityCounter(int, int &);

#endif