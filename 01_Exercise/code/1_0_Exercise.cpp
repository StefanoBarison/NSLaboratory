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

// Now that we set the random generator we can star with the first exercise

//Exercise 01.1 - First point
//Initialize the useful quantities

int M=100000;
int N=100;
int L=M/N;

vector<double> r;
vector<int> x;
vector<double> ave(N,0);
vector<double> av2(N,0);
vector<double> sum_prog(N,0);
vector<double> su2_prog(N,0);
vector<double> err_prog(N,0);

ofstream data1;

//Fill the r vector with pseudo random generated numbers
for(int i=0;i<M;i++){
  r.push_back(rnd.Rannyu());
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
      sum += r[k];
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

//Now we will pass the data to a file
data1.open("r_mean.txt");

for(int i=0;i<N;i++){
  data1<<x[i]<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
}
data1.close();


//Exercise 01.1 - Second point

//Re-define some quantities
M=100000;
N=100;
L=M/N;

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

ofstream data2;

//Fill the r vector with pseudo random generated numbers
for(int i=0;i<M;i++){
  r.push_back(rnd.Rannyu());
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
      sum += pow(r[k]-0.5,2);
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

//Now we will pass the data to a file
data2.open("r_var.txt");

for(int i=0;i<N;i++){
  data2<<x[i]<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
}
data2.close();

//Exercise 01.1 - Third point

//Define some useful quantities

int n= 10000; //Number of throws for each block
M= 100;  // Number of blocks
r.clear(); //Clean up and re-initialize the random number vector
vector<double> succ(M,0);  //Vector of successes in each sub-interval
vector<double> chi(M,0); //Chi-squared
ofstream data3;

for(int i=0;i<M;i++){  //Repeat for every sub-interval

  double id=i; //Create a double from the interval number int

  for(int j=0; j<n;j++){ // Generate n number for each sub-interval
    r.push_back(rnd.Rannyu());
  }

  for(int k=0; k<r.size();k++){
    if(r[k] >= (id)/M && r[k]< (id+1.0)/M) {succ[i]=succ[i]+1;}
  }
  r.clear();
}

//Now we calculate the chi-squared

for(int i=0;i<M;i++){
  chi[i]=pow(succ[i]-(n/M),2)/(n/M);
}


//Let's print the data of the vector successes
data3.open("chi_data.txt");

for(int i=0; i<chi.size()-1;i++){
  data3<<chi[i]<<",";
}
data3<<chi[chi.size()-1];

data3.close();


//Exercise 01.2 Constructing different probability distributions

n=10000; //Number of throws for each N
double temp=0; //useful temporary variable
ofstream data4;
double Nv[]={1,2,10,100}; //Value of the different N

//Standard dice

// Generate the distribution of S_N with N=1,2,10,100 for a standard dice
vector<vector<double>> s(4);


for(int k=0; k<4;k++){
  for(int j=0;j<n;j++){
    temp=0;
    for(int i=0;i<Nv[k];i++){
      temp+=int(rnd.Rannyu(0,6))+1;
      }
    s[k].push_back(temp/Nv[k]);
  }
}

// Generate the distribution of S_N with N=1,2,10,100 for an exponential dice

vector<vector<double>> p(4);

for(int k=0; k<4;k++){
  for(int j=0;j<n;j++){
    temp=0;
    for(int i=0;i<Nv[k];i++){
      temp+=rnd.Exp(1.);
      }
    p[k].push_back(temp/Nv[k]);
  }
}

// Generate the distribution of S_N with N=1,2,10,100 for a Lorentzian dice

vector<vector<double>> l(4);

for(int k=0; k<4;k++){
  for(int j=0;j<n;j++){
    temp=0;
    for(int i=0;i<Nv[k];i++){
      temp+=rnd.Lorentz(0.,1.);
      }
    l[k].push_back(temp/Nv[k]);
  }
}

//Now, let's write all the generated data in a file

data4.open("dice_data.txt");

for(int i=0;i<4;i++){
  for(int j=0;j<s[i].size()-1;j++){
    data4<<s[i][j]<<",";
  }
  data4<<s[i][l[i].size()-1];
  data4<<endl;
}

for(int i=0;i<4;i++){
  for(int j=0;j<p[i].size()-1;j++){
    data4<<p[i][j]<<",";
  }
  data4<<p[i][l[i].size()-1];
  data4<<endl;
}

for(int i=0;i<4;i++){
  for(int j=0;j<l[i].size()-1;j++){
    data4<<l[i][j]<<",";
  }
  data4<<l[i][l[i].size()-1];
  data4<<endl;
}

data4.close();

return 0;
}
