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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  //Write final data
  //Metropolis sampling

  if(print==1 && metro ==1 && h==0.0){
    ofstream ene_f, heat_f, chi_f;
    ene_f.open("metro_ene.dat",ios::app);
    heat_f.open("metro_heat.dat",ios::app);
    chi_f.open("metro_chi.dat",ios::app);

    ene_f << temp<<","<<glob_av[iu]/(double)nblk<<","<< err_u<<endl;
    heat_f << temp<<","<<glob_av[ic]/(double)nblk<<","<< err_c<<endl;
    chi_f << temp<<","<<glob_av[ix]/(double)nblk<<","<< err_x<<endl;
  }
  else if(print==1 && metro ==1 && h==0.02){
    ofstream mag_f;
    mag_f.open("metro_mag.dat",ios::app);
    mag_f << temp<<","<<glob_av[im]/(double)nblk<<","<< err_m<<endl;
  }

  //Gibbs sampling
  else if(print==1 && metro ==0 && h==0.0){
    ofstream ene_f, heat_f, chi_f;
    ene_f.open("gibbs_ene.dat",ios::app);
    heat_f.open("gibbs_heat.dat",ios::app);
    chi_f.open("gibbs_chi.dat",ios::app);

    ene_f << temp<<","<<glob_av[iu]/(double)nblk<<","<< err_u<<endl;
    heat_f << temp<<","<<glob_av[ic]/(double)nblk<<","<< err_c<<endl;
    chi_f << temp<<","<<glob_av[ix]/(double)nblk<<","<< err_x<<endl;
  }
  else if(print==1 && metro ==0 && h==0.02){
    ofstream mag_f;
    mag_f.open("gibbs_mag.dat",ios::app);
    mag_f << temp<<","<<glob_av[im]/(double)nblk<<","<< err_m<<endl;
  }


  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> read;

  ReadInput >> print;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  if(read==0){  //Starting a new run with an initial random configuration
//initial configuration
    cout<<"Equilibration run"<<endl;
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  else if(read==1){ //Here we read froma  previous computed configuration
    cout<< "Measurement run"<<endl;
    ifstream ReadConfig;
    ReadConfig.open("config.final");
    for(int i=0; i<nspin; ++i){
      ReadConfig >> s[i];
    }
    ReadConfig.close();
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double alpha, p, energy_old, energy_new;// sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      energy_old=Boltzmann(s[o],o);
      energy_new=Boltzmann(-1*s[o],o);
      attempted++;

      alpha= exp(beta*(energy_old-energy_new));

      p=rnd.Rannyu();

      if(p<alpha){
        accepted++;
        s[o]=-s[o];
      }
    }
    else if(metro==0)//Gibbs sampling
    {
      attempted++;
      energy_up=Boltzmann(1,o);
      energy_down=Boltzmann(-1,o);
      alpha=1/(1+exp(beta*(energy_up-energy_down)));
      p=rnd.Rannyu();
      if(p <alpha) s[o]=1;
      else s[o]=-1;

      accepted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  //int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u,2); //We accumulate only the squared value, in order to calculate the real value only at the end
  walker[im] = m;
  walker[ix] = pow(m,2);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
   const int wb=3;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wb) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    stima_c = pow(beta,2)*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double) nspin;   //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wb) << iblk <<" "<<  setw(wd) << stima_c <<" "<< setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double) nspin;     //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wb) << iblk <<" "<<  setw(wd) << stima_m <<" "<< setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app);
    stima_x = beta*(blk_av[ix]/blk_norm- pow(blk_av[im]/blk_norm,2))/(double) nspin;    //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wb) << iblk <<" "<<  setw(wd) << stima_x <<" "<< setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
  if(iblk==1) return 0;

  else {
    return sqrt(abs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
  }
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
