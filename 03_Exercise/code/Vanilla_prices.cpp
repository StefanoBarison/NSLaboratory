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

double S(double d_t, double mu ,double sigma,double S1, double rnd  ){

  return S1*exp((mu-pow(sigma,2)/2)*(d_t)+sigma*rnd*sqrt(d_t));

}

double Max(double a, double b){
  if(a>=b) return a;
  else return b;
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

//Now that the random number generator has been set, we can sample the option prices

//1) Sampling directly the price at time T

//Some parameters

   int M=1000000;
   int N=100;
   int L=M/N;  //Parameters for the blocking Monte Carlo

   double T=1.0; 
   double S_0= 100;
   double K=100;
   double rate=0.1;
   double sigma=0.25;   //Parameters for the process

//Now prepare the block estimation

   vector<double> r;
   vector<int> x;
   vector<double> ave(N,0);
   vector<double> av2(N,0);
   vector<double> sum_prog(N,0);
   vector<double> su2_prog(N,0);
   vector<double> err_prog(N,0);

   ofstream data1;

//Fill the vector with gaussian sampled random number

for(int i=0; i<M; i++){
  r.push_back(rnd.Gauss(0,T));
}

   //Fill the x vector with the number of throws
for (int i = 0; i < N; i++) {
  x.push_back(i*L);
}
// Calculate the average value for each block of throws
for(int i=0;i<N;i++){
    double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += exp(-1*rate*T)*Max(0,S(T,rate,sigma,S_0,r[k])-K);
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


//Print the data for direct call

 data1.open("Direct_call.txt");

 for(int i=0;i<N;i++){
  data1<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
 }

 data1.close();

 //Now let's make the direct put option

r.clear();
x.clear();
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

//Fill the vector with gaussian sampled random number

for(int i=0; i<M; i++){
  r.push_back(rnd.Gauss(0,T));
}

   //Fill the x vector with the number of throws
for (int i = 0; i < N; i++) {
  x.push_back(i*L);
}
// Calculate the average value for each block of throws
for(int i=0;i<N;i++){
    double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += exp(-1*rate*T)*Max(0,K-S(T,rate,sigma,S_0,r[k]));
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


//Print the data for direct pull

 data1.open("Direct_pull.txt");

 for(int i=0;i<N;i++){
  data1<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
 }

 data1.close();


// Now we can pass to the iterative resolution of the problem

// Iterative call

r.clear();
x.clear();
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

//Fill the x vector with the number of throws
for (int i = 0; i < N; i++) {
  x.push_back(i*L);
}

//Now the total time is the time of a step

double d_t=1/100.;

// Calculate the average value for each block of throws
for(int i=0;i<N;i++){
    double sum=0;
    for(int j=0;j<L;j++){

      for(int k=0;k<100;k++){
        r.push_back(rnd.Gauss(0,1));
      }

      double tS=S_0;
      for(int k=0;k<100;k++){
        tS= S(d_t,rate,sigma, tS,r[k]);
      }
      r.clear();

      sum += exp(-1*rate*T)*Max(0,tS-K);
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

//Print the data for iterative call

 data1.open("Iterative_call.txt");

 for(int i=0;i<N;i++){
  data1<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
 }

 data1.close();


 //Iterative pull

r.clear();
x.clear();
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


//Fill the x vector with the number of throws
for (int i = 0; i < N; i++) {
  x.push_back(i*L);
}

// Calculate the average value for each block of throws
for(int i=0;i<N;i++){
    double sum=0;
    for(int j=0;j<L;j++){

      for(int k=0;k<100;k++){
        r.push_back(rnd.Gauss(0,1));
      }

      double tS=S_0;
      for(int k=0;k<100;k++){
        tS=S(d_t,rate,sigma,tS,r[k]);
      }

      r.clear();
      sum += exp(-1*rate*T)*Max(0,K-tS);
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

//Print the data for iterative call

 data1.open("Iterative_pull.txt");

 for(int i=0;i<N;i++){
  data1<<x[i]<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
 }

 data1.close();
return 0;
}