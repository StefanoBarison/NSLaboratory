
//parameters, observables
const int m_props=1000;
int n_props,igofr,iv,ik,it,iw,ie;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,stima_temp,stima_g[m_props],err_pot,err_press,err_temp,err_g;


//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk,seed,restart;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfOld(void);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
