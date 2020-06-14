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
double wave_1s(double x,double y, double z);
double wave_2p(double x,double y, double z);

int main(int argc, char *argv[]){

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

//The first thing to do is to "equilibrate" the system, namely to choose the correct delta for the uniform probability distribution
// in order to obtain an acceptance ratio of 0.5 or similar

int N_trials=1000;
double t_x=0;
double t_y=0;
double t_z=0;
double delta=0.1;
double ratio=0.5;
double acceptance=0.005;
double t_ratio=1;
int n_accept=0;


cout<<"First integral"<<endl;
cout<<"Equilibration"<<endl;

while(abs(ratio-t_ratio)>acceptance){ //This cycle let us equilibrate the system, in order to find the desired delta
  double r_x,r_y,r_z,r_p,r_r;
  n_accept=0;

  for(int i=0; i<N_trials;i++){
    //Standard equilibration 
    /*
    r_x=t_x-delta+2*delta*rnd.Rannyu();
    r_z=t_z-delta+2*delta*rnd.Rannyu();
    r_z=t_z-delta+2*delta*rnd.Rannyu();
    */

    //Gaussian equilibration
    r_x=t_x+rnd.Gauss(0,delta);
    r_y=t_y+rnd.Gauss(0,delta);
    r_z=t_z+rnd.Gauss(0,delta);

    r_p=wave_1s(r_x,r_y,r_z)/wave_1s(t_x,t_y,t_z);
    r_r=rnd.Rannyu();

    if(r_p>=r_r){
      t_x=r_x;
      t_y=r_y;
      t_z=r_z;

      n_accept++;
    }
  }

  t_ratio= (double)n_accept/(double)N_trials;

  //cout<<"ratio:"<<t_ratio<<endl;

  if(abs(ratio-t_ratio)>acceptance){
    if(t_ratio>ratio){
      delta+=0.01;
    }
    else if(t_ratio<ratio){
      delta-=0.01;
    }
  }

  //cout<<"delta= "<<delta<<endl;

}

//Now we have the right delta, without saving any big vector.

cout<< "System equilibrated with delta= "<<delta<<" and a ratio of acceptance of "<<t_ratio<<endl;

//We can proceed to the real experiment




int N= 1000000; //Number of throws
int M= 100; //Number of blocks
//int L= N/M; // Throws per block
double p; //Ratio between probabilities
double acc_r; //

vector<double> x(N,0);
vector<double> y(N,0);
vector<double> z(N,0);
vector<double> r(N,0);

//Initialize the vectors
cout<<"Initialize vectors"<<endl;

x[0]=0.0;
y[0]=0.0;
z[0]=0.0;
r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));

//Let's calculate the values for each step

cout<<"Simulate with metropolis algorithm"<<endl;

for(int i=0;i<N-1;i++){
  //cout<<"ok"<<i<<endl;
  /*
  t_x=x[i]-delta+ 2*delta*rnd.Rannyu();
  t_y=y[i]-delta+ 2*delta*rnd.Rannyu();
  t_z=z[i]-delta+ 2*delta*rnd.Rannyu();
  */
  t_x=x[i]+rnd.Gauss(0,delta);
  t_y=y[i]+rnd.Gauss(0,delta);
  t_z=z[i]+rnd.Gauss(0,delta);

  p=wave_1s(t_x,t_y,t_z)/wave_1s(x[i],y[i],z[i]);

  acc_r=rnd.Rannyu();

  if(p>acc_r){
    x[i+1]=t_x;
    y[i+1]=t_y;
    z[i+1]=t_z;
    r[i+1]= sqrt(pow(x[i+1],2)+pow(y[i+1],2)+pow(z[i+1],2));
  }
  else{
    x[i+1]=x[i];
    y[i+1]=y[i];
    z[i+1]=z[i];
    r[i+1]=r[i];
  }
}

//Now proceed to the data blocking in order to estimate the integral

cout<<"Blocking..."<<endl;

vector<double> ave=Blocking(M,r);
vector<double> sum_prog(M,0);
vector<double> err_prog(M,0);

DataBlocking(ave,sum_prog,err_prog);

// And print the data

cout<<"Printing data..."<<endl;
//ofstream points;
ofstream integral;

//points.open("1s_points.dat");
integral.open("1s_integral_gauss.dat");

/*
for(int i=0;i<N;i++){
  points<<x[i]<<","<<y[i]<<","<<z[i]<<endl;
}
*/
for(int i=0;i<M;i++){
  integral<<i+1<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
}

//points.close();
integral.close();


//Now we proceed to the evaluation of the second integral

// Let's clear all the importart variables

x.clear();
y.clear();
z.clear();
ave.clear();
sum_prog.clear();
err_prog.clear();
 
N_trials=1000;
t_x=0;
t_y=0;
t_z=0;
delta=0.1;
ratio=0.5;
acceptance=0.005;
t_ratio=1;
n_accept=0;

cout<<endl<<endl;
cout<<"Second integral"<<endl;
cout<<"Equilibration"<<endl;

while(abs(ratio-t_ratio)>acceptance){ //This cycle let us equilibrate the system, in order to find the desired delta
  double r_x,r_y,r_z,r_p,r_r;
  n_accept=0;

  for(int i=0; i<N_trials;i++){
    /*
    r_x=t_x-delta+2*delta*rnd.Rannyu();
    r_z=t_z-delta+2*delta*rnd.Rannyu();
    r_z=t_z-delta+2*delta*rnd.Rannyu();
    */

    r_x=t_x+rnd.Gauss(0,delta);
    r_y=t_y+rnd.Gauss(0,delta);
    r_z=t_z+rnd.Gauss(0,delta);

    r_p=wave_2p(r_x,r_y,r_z)/wave_2p(t_x,t_y,t_z);
    r_r=rnd.Rannyu();

    if(r_p>=r_r){
      t_x=r_x;
      t_y=r_y;
      t_z=r_z;

      n_accept++;
    }
  }

  t_ratio= (double)n_accept/(double)N_trials;

  //cout<<"ratio:"<<t_ratio<<endl;

  if(abs(ratio-t_ratio)>acceptance){
    if(t_ratio>ratio){
      delta+=0.01;
    }
    else if(t_ratio<ratio){
      delta-=0.01;
    }
  }

  //cout<<"delta= "<<delta<<endl;

}

//Now we have the right delta, without saving any big vector.

cout<< "System equilibrated with delta= "<<delta<<" and a ratio of acceptance of "<<t_ratio<<endl;

//We can proceed to the real experiment




N= 1000000; //Number of throws
M= 100; //Number of blocks

x=vector<double>(N,0);
y=vector<double>(N,0);
z=vector<double>(N,0);
r=vector<double>(N,0);

//Initialize the vectors
cout<<"Initialize vectors"<<endl;

x[0]=0.0;
y[0]=0.0;
z[0]=0.0;
r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));

//Let's calculate the values for each step

cout<<"Simulate with metropolis algorithm"<<endl;

for(int i=0;i<N-1;i++){
  //cout<<"ok"<<i<<endl;
  /*
  t_x=x[i]-delta+ 2*delta*rnd.Rannyu();
  t_y=y[i]-delta+ 2*delta*rnd.Rannyu();
  t_z=z[i]-delta+ 2*delta*rnd.Rannyu();
  */

  t_x=x[i]+rnd.Gauss(0,delta);
  t_y=y[i]+rnd.Gauss(0,delta);
  t_z=z[i]+rnd.Gauss(0,delta);

  p=wave_2p(t_x,t_y,t_z)/wave_2p(x[i],y[i],z[i]);

  acc_r=rnd.Rannyu();

  if(p>acc_r){
    x[i+1]=t_x;
    y[i+1]=t_y;
    z[i+1]=t_z;
    r[i+1]= sqrt(pow(x[i+1],2)+pow(y[i+1],2)+pow(z[i+1],2));
  }
  else{
    x[i+1]=x[i];
    y[i+1]=y[i];
    z[i+1]=z[i];
    r[i+1]=r[i];
  }
}

//Now proceed to the data blocking in order to estimate the integral

cout<<"Blocking..."<<endl;

ave=Blocking(M,r);
sum_prog=vector<double>(M,0);
err_prog=vector<double>(M,0);

DataBlocking(ave,sum_prog,err_prog);

// And print the data

cout<<"Printing data..."<<endl;

//points.open("2p_points.dat");
integral.open("2p_integral_gauss.dat");

/*
for(int i=0;i<N;i++){
  points<<x[i]<<","<<y[i]<<","<<z[i]<<endl;
}
*/
for(int i=0;i<M;i++){
  integral<<i+1<<","<<sum_prog[i]<<","<<err_prog[i]<<endl;
}

//points.close();
integral.close();

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

double wave_1s(double x,double y, double z){
    return exp(-2*sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/M_PI;
}

double wave_2p(double x,double y, double z){
    return pow(z,2)*exp(-sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/(32*M_PI);
}