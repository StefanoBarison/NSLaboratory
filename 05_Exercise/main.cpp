#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

double error(vector<double> & av,vector<double> & av2, int n);
vector<double> Blocking(int nblocks,vector<double> & v);
void DataBlocking(vector <double> value, vector <double> &sum_prog, vector <double> &err_prog);

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


// Now we can star with the exercise
//We want to calculate the expectation value of r over the 1s orbital of hydrogen

int N= 100000; //Number of throws
int M= 1000; //Number of blocks
int L= N/M; // Throws per block


return 0;
}



// Useful functions

double error(vector<double> & av,vector<double> & av2, int n){  // Error function used to calculate
                                                                //statistical uncertainties
    if (n==0){
      return 0;
    }
    else{
      return sqrt((av2[n]-pow(av[n],2))/n);
    }
}

vector<double> Blocking(int nblocks,vector<double> & v){
 vector<double> b_v;
 int N= v.size();
 int L= N/nblocks;
 for(int i=0;i<nblocks;i++){
     double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += v[k];
      }
    b_v.push_back(sum/L);
  }
  return b_v;
}

void DataBlocking(vector <double> value, vector <double> &sum_prog, vector <double> &err_prog){
    int N=sum_prog.size();
    vector <double> su2_prog(N,0);
    
    for(int i=0;i<N;i++){
        for(int j=0;j<i+1;j++){
            sum_prog[i]+=value[j];
            su2_prog[i]+=pow(value[j],2);
        }
        sum_prog[i]=sum_prog[i]/(i+1);      //Cumulative average
        su2_prog[i]=su2_prog[i]/(i+1);      //Cumulative square average
        err_prog[i] = error(sum_prog,su2_prog,i);
    }
}