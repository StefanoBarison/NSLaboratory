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

class Needle{
  private:
    double _x,_l,_theta;
  public:
    Needle(){
      _x=0;
      _l=0;
      _theta=0;
    }
    Needle(double x, double theta, double l){
      _x=x; //position of the mean point of the needle wrt the nearest line, so [0,d/2]
      _l=l; //lenght of the needle
      _theta=theta; //Angle wrt the x-axis
    }

    double Check(){ //This function returns 1 if the Needle crosses a line, 0 otherwise
      int count=0;
      if(_x<=abs(_l*cos(_theta)/2)){
          count++;
        }
      
      double dcount=count*1.;
      return dcount;
    }
    void Print(){
      cout<<"x mean= "<<_x<<", Angle "<<_theta<<endl;
    }
};

int main (int argc, char *argv[]){

//Parrallel random number generator

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

// Now we will initialize the useful variables for the Buffon experiment.
// In this case we have a 1x1 square divided by n-1 lines in n sub-interval.
// Each line is distant d=1/n from their nearest neighbors
int n=10; //Number of sub-interval
double d=1./n; //Distance of the bars
double l= 0.06; //Length of the needle
int M=1000000; //Number of throws
int N=100; //Number of outcome
int L=M/N; //Number of throws per outcome

//Now we generate a vector of N random x and a vector of N random _cos

vector<double> x;
vector<double> theta;
vector<Needle> t;

for(int i=0;i<M;i++){
  x.push_back(rnd.Rannyu(0,d/2));
}

for(int i=0;i<M;i++){
  double a=rnd.Rannyu();
  double b=rnd.Rannyu();

  while(pow(a,2)+pow(b,2)>=1){

    a=rnd.Rannyu();
    b=rnd.Rannyu();
  }

    //theta.push_back(rnd.Rannyu(0,M_PI/2));
    theta.push_back(acos( a / (sqrt(pow(a,2)+pow(b,2)))  ) );

}

//Now we create a vector of needles
for(int i=0;i<M;i++){
  t.push_back(Needle(x[i],theta[i],l));
}

//And then also a vector of Check
vector<double> check;

for (int i=0; i<M; i++) {
  check.push_back(t[i].Check());
}

//Now that we have a vector of checks, we can estimate pi for every block

vector<double> ave(N,0);
vector<double> av2(N,0);
vector<double> sum_prog(N,0);
vector<double> su2_prog(N,0);
vector<double> err_prog(N,0);

for(int i=0;i<N;i++){
    double sum=0;
    for(int j=0;j<L;j++){
      int k=j+i*L;
      sum += check[k];
      ave[i]= (2*L*l)/(sum*d);
      av2[i]= pow(ave[i],2);
    }
}

for(int i=0;i<N;i++){
    for(int j=0;j<i+1;j++){
        sum_prog[i] +=ave[j];
        su2_prog[i] += av2[j];
    }

  sum_prog[i] /= (i+1);
  su2_prog[i] /= (i+1);
  err_prog[i] = error(sum_prog,su2_prog,i);
}

//for(int i=0;i<N;i++){
  //t[i].Print();
  //cout<<sum_prog[i]<<endl;
//}

//Now let's create a file with all the data to plot in the jupyter notebook

ofstream data;
data.open("Buffon.txt");

//vector with the # of throws
vector<int> x_t;

for (int i = 0; i < N; i++) {
  x_t.push_back(i*L);
}

for(int i=0;i<N;i++){
  data<<x_t[i]<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
}

data.close();


return 0;
}
