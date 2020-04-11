#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
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

double Average(std::vector<double> const& v) {
    return 1.0 * std::accumulate(v.begin(), v.end(), 0LL) / v.size();
}

double Std_Error(std::vector<double> const& v){
  double av=Average(v);
  double sum=0;
  double var=0;
  int n=v.size();

  for(int i=0;i<n;i++){
    sum += pow((v[i]-av),2); 
  }
  var= sum/(n-1.);

  return sqrt(var);
}

class Unif_walker
{
private:
  double _x,_y,_z;

public:
  Unif_walker(){
    _x=0;
    _y=0;
    _z=0;
  }
  ~Unif_walker(){

  }

  Unif_walker(std::vector<double> pos){
    _x=pos[0];
    _y=pos[1];
    _z=pos[2];
  }

  vector<double> Position(){
    std::vector<double> v={_x,_y,_z};

    return v;
  }

  void Print(){
    cout<<_x<<","<<_y<<","<<_z<<endl;
  }

  double S_distance(){
    return pow(_x,2)+pow(_y,2)+pow(_z,2);
  }

  void Set_Position(vector<double> v){
    _x=v[0];
    _y=v[1];
    _z=v[2];
  }

  void Walk(double theta, double phi){
    _x= _x+ sin(theta)*cos(phi);
    _y= _y+ sin(theta)*sin(phi);
    _z= _z+ cos(theta);
  }
};

class Latt_Walker
{

private:
  double _a;
  int _x,_y,_z;

public:

  Latt_Walker(){
    _x=0;
    _y=0;
    _z=0;
    _a=1;
  }
  ~Latt_Walker(){

  }
  
  Latt_Walker(double a, std::vector<int> pos){
    _a=a;
    _x=pos[0];
    _y=pos[1];
    _z=pos[2];
  }

  void Walk_x(double sign){
    if(sign>= 0.5){
      _x++;
    }
    else{_x=_x-1;}
  }
  void Walk_y(double sign){
    if(sign>= 0.5){
      _y++;
    }
    else{_y=_y-1;}
  }

  void Walk_z(double sign){
    if(sign>= 0.5){
      _z++;
    }
    else{_z=_z-1;}
  }

  void Set_Position(vector<int> v){
    _x=v[0];
    _y=v[1];
    _z=v[2];
  }

  vector<int> Position(){
    std::vector<int> v={_x,_y,_z};

    return v;
  }

  void Print(){
    cout<<_x<<","<<_y<<","<<_z<<endl;
  }

  double S_distance(){
    return pow(_x,2)+pow(_y,2)+pow(_z,2);
  }

};


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


   //The first part of the exercise consist in creating a random walker on a lattice and then 
   //test the diffusive behavior

   int M=10000; // # of walkers;
   int l=100; // # of steps
   vector<int> origin={0,0,0} ;//The origin position
   vector<int> r; //Vector of the random direction chosen
   vector<double> s; //Vector of the random sign chosen
   double a=1; //Lattice length
   vector<vector<double>> rn2; //Vector of the different positions at different steps
   ofstream data1; //Out data of the final position


   for(int i=0;i<M;i++){
      for(int j=0;j<l;j++){
        r.push_back(int(rnd.Rannyu(0,3)+1));
      }

      for(int k=0;k<l;k++){
        s.push_back(rnd.Rannyu());
      }

      Latt_Walker w(a,origin);

      vector<double> temp;

      for(int j=0;j<l;j++){
        if(r[j]==1) w.Walk_x(s[j]);
        else if(r[j]==2) w.Walk_y(s[j]);
        else if(r[j]==3) w.Walk_z(s[j]);
        temp.push_back(w.S_distance());
      }


      rn2.push_back(temp);

      r.clear();
      s.clear();
      temp.clear();
    }
  
    
    //Now we calculate mean and std dev of every vector

    vector<double> mean;
    vector<double> std;

    for(int i=0;i<l;i++){
      vector<double> temp1;
      for(int j=0;j<M;j++){
        temp1.push_back(rn2[j][i]);
      }

      double t_av=Average(temp1);
      double t_std=Std_Error(temp1)/sqrt(M);

      mean.push_back(sqrt(t_av));
      std.push_back(t_std/(2*sqrt(t_av)));

      temp1.clear();
     }

     //Now print the data in an output file

    data1.open("Lattice_walk.txt");

    for(int j=0;j<l;j++){
      data1<<mean[j]<<","<<std[j]<<endl;
    }

    data1.close();


//Now we will concentrate on a random walker that moves uniformly in the 3D space

//Clear all the random number vector in case we want to re-use them
    r.clear();
    s.clear();
    for(int i=0;i<M;i++){
      rn2[i].clear();
    }
    rn2.clear();

//Introduce new variables

    vector<double> theta;
    vector<double> phi;
    vector<double> o={0,0,0};
    vector<vector<double>> rnu2;
    Unif_walker u(o);

    for(int i=0;i<M;i++){
      for(int j=0;j<l;j++){
        phi.push_back(rnd.Rannyu(0,2*M_PI));
      }
      for(int j=0;j<l;j++){
        theta.push_back(acos(1-2*rnd.Rannyu()));
      }

      vector<double> temp2;

      for(int j=0;j<l;j++){
        u.Walk(theta[j],phi[j]);
        temp2.push_back(u.S_distance());
      }

      rnu2.push_back(temp2);

      theta.clear();
      phi.clear();
      temp2.clear();
      u.Set_Position(o);
    }

    vector<double> mean2;
    vector<double> std2;

    for(int i=0;i<l;i++){
      vector<double> temp3;

      for(int j=0;j<M;j++){
        temp3.push_back(rnu2[j][i]);
      }

      double t_av=Average(temp3);
      double t_std=Std_Error(temp3)/sqrt(M);

      mean2.push_back(sqrt(t_av));
      std2.push_back(t_std/(2*sqrt(t_av)));

      temp3.clear();
     }

     //Now print the data in an output file

    data1.open("Uniform_walk.txt");

    for(int j=0;j<l;j++){
      data1<<mean2[j]<<","<<std2[j]<<endl;
    }

    data1.close();


   return 0;
}