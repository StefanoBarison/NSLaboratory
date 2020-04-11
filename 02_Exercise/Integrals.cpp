#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

double error(vector<double> & av,vector<double> & av2, int n){  // Error function used to calculate
                                                                //statistical uncertainties
    if (n==0){
      return 0;
    }
    else{
      return sqrt((av2[n]-pow(av[n],2))/n);
    }
}
double Integrand(double x){
  double y= (M_PI/2)*cos(x*M_PI/2);
  return y;
}
double Imp_Integrand(double x){
  double y= (M_PI/4)*cos(x*M_PI/2)/(1-x);
  return y;
}

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


// Now we can star with the first exercise
// Exercise 02.1 - First point
// We are going to calculate an integral sampling from a uniform distribution in [0,1]

   int M= 1000000; //Number of throws
   int N= 100; //Number of blocks
   int L=M/N;

  vector<double> r;  //Vector of random numbers
  vector<int> x;     //Vector of throws
  vector<double> ave(N,0);
  vector<double> av2(N,0);
  vector<double> sum_prog(N,0);
  vector<double> su2_prog(N,0);
  vector<double> err_prog(N,0);

  ofstream data1;

//Fill the random number vector
  for (int i = 0; i < M; ++i)
  {
    r.push_back(rnd.Rannyu());
  }

//Fill the x vector with the number of throws
  for (int i = 0; i < N; i++) {
    x.push_back(i*L);
  }

//Calculate the average of the funciton for each bloch (the estimator of the integral)
  for (int i = 0; i < N; ++i){
    double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += Integrand(r[k]);
      ave[i]= sum/L;
      av2[i]= pow(ave[i],2);
    }
  }

//Calculate the cumulative mean and the error
  for(int i=0;i<N;i++){
      for(int j=0;j<i+1;j++){
          sum_prog[i] +=ave[j];
          su2_prog[i] += av2[j];
      }

    sum_prog[i] /= (i+1);
    su2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog,su2_prog,i);
  }

  data1.open("integral_uniform.txt");

  for(int i=0;i<N;i++){
      data1<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
  }
  data1.close();

//Now we can calculate the same integral with the importance sampling technique
// We can use the function 2(1-x) which is more similar to the integrand, to generate the random number
// The function to use to generate random number sampled from that probability is x=1-sqrt(1-y)

  //Let's clear all the vector and create a new data file
  r.clear();
  ave.clear();
  av2.clear();
  sum_prog.clear();
  su2_prog.clear();
  err_prog.clear();

  ave= vector<double>(N,0);
  av2= vector<double>(N,0);
  sum_prog=vector<double>(N,0);
  su2_prog=vector<double>(N,0);
  err_prog=vector<double>(N,0);

  ofstream data2;

  //Generate the random numbers
  for (int i = 0; i < M; ++i)
  {
    r.push_back(1-sqrt(1-rnd.Rannyu()));
  }

  //Calculate the average of the funciton for each bloch (the estimator of the integral)
  for (int i = 0; i < N; ++i){
    double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += Imp_Integrand(r[k]);
      ave[i]= sum/L;
      av2[i]= pow(ave[i],2);
    }
  }

//Calculate the cumulative mean and the error
  for(int i=0;i<N;i++){
      for(int j=0;j<i+1;j++){
          sum_prog[i] +=ave[j];
          su2_prog[i] += av2[j];
      }

    sum_prog[i] /= (i+1);
    su2_prog[i] /= (i+1);
    err_prog[i] = error(sum_prog,su2_prog,i);
  }

  data2.open("integral_importance.txt");

  for(int i=0;i<N;i++){
      data2<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
  }
  data2.close();

  return 0;
}