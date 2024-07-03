/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
***************************************************************** 
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;

   cout << "cout 1" << endl;

   ifstream Primes("Primes");

   cout << "cout 2" << endl;

   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;

   cout << "cout 3" << endl;

   Primes.close();

   cout << "cout 4" << endl;

   ifstream input("seed.in");

   cout << "cout 5" << endl;

   string property;
   if (input.is_open()){

      cout << "input opened" << endl;

      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){

            cout << "RANDOMSEED found" << endl;

            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];

            cout << "seeds loaded" << endl;
            cout << seed[0] << seed[1] << seed[2] << seed[3] << endl;

            rnd.SetRandom(seed,p1,p2);
         }
      }

      cout << "eof" << endl;

      input.close();

      cout << "input closed" << endl;

   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   cout << "cout 6" << endl;

   ofstream out;
   out.open("number01.2.out");

   cout << "cout 7" << endl;

   out << "standard" << " Exponential" << " Lorentz" << endl;

   cout << "cout 8" << endl;

   for(int i=0; i<1E5; i++){
    out << rnd.Rannyu() << " " << rnd.Exponential(1) << " " << rnd.Lorentz(0, 1) << endl;
   }

   cout << "cout 9" << endl;

   out.close();

   cout << "cout 10" << endl;

   // for(int i=0; i<20; i++){
   //    cout << rnd.Rannyu() << endl;
   // }

   rnd.SaveSeed();

   cout << "cout 11" << endl;

   return 0;
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
