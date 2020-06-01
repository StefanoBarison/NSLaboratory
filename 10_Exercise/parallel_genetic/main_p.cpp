#include <iostream>
#include <fstream>
#include <vector>

#include "tsp.h"
#include "Random.h"
#include "mpi.h"

using namespace std;

int main(int argc, char** argv){

	//Initialize MPi and get size and ranks
	MPI::Init(argc,argv);

	int size=MPI::COMM_WORLD.Get_size();
	int rank=MPI::COMM_WORLD.Get_rank();

	//Useful quantities for all cores

	int n_cities=32;
	int n_steps=2000;
	int n_ind=1000;
	int n_migr=100;
	string type ("Elitism");
	
	ifstream indata;
	indata.open("square.dat");

	//Now initialize random number generators with different primes for each core

	Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
  	for(int i=0;i<rank+1;i++){
  		Primes >> p1 >> p2 ;
  	}
  	cout<<p1<<" "<<p2<<endl;
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

  //Now we have initialized random number generators with different primes, we can start doing the GA

    ////Initialization

  	map m(n_cities);
  	m.File_initialize(indata);
  	m.Create_d_matrix();

  	population p(n_ind,m);
	p.Initialize(&rnd);
	p.Evaluate_all(m);

	/////Vector to save the results

	vector<double> best_result;
	vector<double> best_ave;

	////Vector to perform migration
	vector<int> travel;
	vector<int> emigrant;
	vector<int> immigrant;

	for(int j=0;j<size;j++){
		travel.push_back(j);
	}

	/////Perform the GA
	for(int i=0;i<n_steps;i++){
		
		//Here we have to perform the migration

		if(i>0 && i%n_migr==0){
			if(rank==0){
				cout<<"Permuting continents..."<<endl;
				random_shuffle(travel.begin(),travel.end());
			}
			
			MPI_Bcast(travel.data(),travel.size(),MPI_INTEGER,0,MPI::COMM_WORLD);
			
			//Preparing messages
			emigrant=p.Get_best(m)->Get_genes();
			immigrant.resize(n_cities);
			for(int k=0;k<n_cities;k++){
				immigrant[k]=0;
			}
			//Sending the messages
			for(int j=0;j<size;j++){
				if(j==size-1){
					if(rank==travel[j]){
						/*cout<<"Sending to "<<travel[0]<<endl;
						for(unsigned int l=0;l<emigrant.size();l++){
							cout<<emigrant[l]<<" ";
						}
						cout<<endl;*/
						MPI_Send(emigrant.data(),emigrant.size(),MPI_INT,travel[0],0,MPI_COMM_WORLD);

					}
				}
				else{
					if(rank==travel[j]){
						MPI_Send(emigrant.data(),emigrant.size(),MPI_INT,travel[j+1],0,MPI_COMM_WORLD);
					}
				}
			}
			// Receiving the messages

			for(int j=1;j<size;j++){
				if(rank==travel[j]){
					MPI_Recv(immigrant.data(),immigrant.size(),MPI_INT,travel[j-1],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			}
			if(rank==travel[0]){
				MPI_Recv(immigrant.data(),immigrant.size(),MPI_INT,travel[size-1],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				/*cout<<"Received from "<<travel[size-1]<<endl;
				for(unsigned int l=0;l<immigrant.size();l++){
							cout<<immigrant[l]<<" ";
						}
				cout<<endl;*/
			}

			//Now replacing the best of the population with the immigrant

			p.Get_best(m)->Set_genes(immigrant);
		}

		//cout<<"Step "<<i+1<<endl<<endl;
		best_result.push_back(p.Get_best(m)->Get_lenght());
		best_ave.push_back(p.Get_best_average(m));
		p.Evolutive_step(m,type,&rnd);

		if(i==n_steps-1){
			if(rank==0){
				p.Print_best(m);
			}
		}
	}


	////Now print the data on different files
	ofstream outdata1, outdata2;

	string best="best_square_"+to_string(rank)+".dat";
	string ave="average_square_"+to_string(rank)+".dat";

	outdata1.open(best);
	outdata2.open(ave);

	for(unsigned int i=0;i<best_result.size();i++){
		outdata1<<i<<","<<best_result[i]<<endl;
		outdata2<<i<<","<<best_ave[i]<<endl;
	}

	MPI::Finalize();

	return 0;
}