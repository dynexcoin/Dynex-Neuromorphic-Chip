// Copyright (c) 2021-2022, The DYNEX Project
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, this list
//    of conditions and the following disclaimer in the documentation and/or other
//    materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its contributors may be
//    used to endorse or promote products derived from this software without specific
//    prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// REQUIREMENTS:
// brew install boost (macos)
//
// COMPILE:                
// g++ dynex.cc -o dynex -std=c++17 -Ofast -lpthread -fpermissive
// macos:                  g++ dynex.cc -o dynex -std=c++17 -Ofast -I /opt/homebrew/cellar/boost/1.78.0/include -L /opt/homebrew/cellar/boost/1.78.0/lib
//
// RUN:
// ./dynex -i cnf/transformed_barthel_n_100_r_8.000_p0_0.080_instance_001.cnf
// ./dynex -i cnf/transformed_barthel_n_200_r_8.000_p0_0.080_instance_001.cnf
// ./dynex -i cnf/transformed_barthel_n_500_r_8.000_p0_0.080_instance_001.cnf
// ./dynex -i cnf/transformed_barthel_n_1000_r_8.000_p0_0.080_instance_001.cnf
// ./dynex -i cnf/transformed_barthel_n_10000_r_8.000_p0_0.080_instance_001.cnf
// ./dynex -i cnf/transformed_barthel_n_100000_r_8.000_p0_0.080_instance_001.cnf


#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <signal.h>
#include <math.h>
#include <stdbool.h>
#include <locale.h>
#include <random>
#include <iostream>

/* BOOST */ /// 
#include <iostream>
#include <vector> //?
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/clamp.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace pt = boost::posix_time;

//--------------------------------------------------------------------------------------------------------------
#define         MAX_LITS_SYSTEM 25 //10

#define         PARAM_COUNT     9       //number of parameters
#define         LBH     0               //lower bound hard
#define         HBH     1               //higher bound hard
#define         LBS     2               //lower bound soft
#define         HBS     3               //higher bound soft
#define         ALPHA   0
#define         BETA    1
#define         GAMMA   2
#define         DELTA   3
#define         EPSILON 4
#define         ZETA    5
#define         ERR1    6
#define         ERR2    7
#define         INITDT  8

#define         ODE_CONSTANT        1
#define         ODE_CUSTOM_ADAPTIVE 2
#define         ODE_RUNGEKUTTA      3
#define         ODE_IMPLICIT        4
#define         ODE_ONE_STEP        5
#define         ODE_CONSTANT_NOBOOST 6

#define         ODE_LOGFILE         "log.csv"
#define         TUNE_LOGFILE        "tuninglog.csv"
#define         PARTABLE_FILE       "partable.txt"
#define         FLOWVECTOR_FILE     "flowvector.csv"

char            LOC_FILE[256];
char            SOLUTION_FILE[256];
char            TUNING_FILE[256];
char            ASSIGNMENT_FILE[256];
//--------------------------------------------------------------------------------------------------------------
int THREAD_COUNT = 8;
volatile int running_threads = 0;
//--------------------------------------------------------------------------------------------------------------
// PRECISION OF ODE INTEGRATION
//--------------------------------------------------------------------------------------------------------------
//typedef float TFloat;
typedef double TFloat;
//typedef long double TFloat;

#define USEFPRINTF
int digits = 15; // precision for ourputting numbers: 7 for float; 15 for double

//--------------------------------------------------------------------------------------------------------------
typedef TFloat value_type;
typedef std::vector< value_type > state_type; //we use a vector which has dynamic size
//--------------------------------------------------------------------------------------------------------------
/* SETTINGS & PARAMETERS */
bool   quiet           = false; //no screen output
bool   flowvector_log  = false; //output entire FLOWVECTOR_FILE (only thread 0)
bool   load_partable   = false;  //load partable.txt
bool   writelogfile    = false; //write logfile ODE_LOGFILE
bool   writelocfile    = false; //write all LOC to file
bool   writesolution   = false; //write solution to file during solving
bool   load_solution   = false; //true; //load solution.txt at startup
bool   load_assignment = false; //load .assignemnt.txt at startup (V,Xs,Xl,t_init)
bool   massive         = false; // multiple runs, each with different initial Xs, Xl
bool   constand_adaptive = true; // use adaptive timestep in ODE_CONSTANT_NOBOOST

// Equations - constants:
TFloat dmm_alpha      = TFloat(5.0);    // growth rate for long term memory Xl
TFloat                * dmm_alpha_cm;   // alpha for each clause - used for heuristics
TFloat dmm_beta       = TFloat(20.0);   // growth rate for short term memory Xs
TFloat dmm_gamma      = TFloat(TFloat(25)/TFloat(100)); // TFloat(0.25);   // restriction for Cm in short term memory Xs
TFloat dmm_delta      = TFloat(TFloat(5)/TFloat(100)); //TFloat(0.05);   // restriction for Cm in long term memory Xl
TFloat dmm_epsilon    = TFloat(TFloat(1)/TFloat(10)); //TFloat(0.1);    // remove spurious solution X,s,m = 0
TFloat dmm_zeta       = TFloat(TFloat(1)/TFloat(10)); //TFloat(0.1);    // reduction factor of rigidity G (learning rate) 10^-3; for ratio>=6: 10^-1 (0.1)
int    seed           = 0;      // random seed value (for initial assignment); defaults to 0,1,-1,rand,rand,rand...
int    xl_max         = 10000;  //10^4 M (x count clauses will be applied automatically - ODE should NEVER reach this value)

// ODE settings: 
int INTEGRATION_MODE = ODE_CONSTANT_NOBOOST; // ODE_CONSTANT; // ODE_RUNGEKUTTA;
bool forcebounds     = true; // update vector x: apply boundaries after every step 

// Constant integration params:
TFloat timeout              = TFloat(INT_MAX); //max simulated time; stops at reaching it
double walltime_timeout     = INT_MAX; //max walltime time; stops at reaching it

// Runge-Kutta Adaptive params:
TFloat rk_errorrate_1       = TFloat(TFloat(52)/TFloat(100)); //TFloat(0.52);  
TFloat rk_errorrate_2       = TFloat(TFloat(10)/TFloat(100)); //TFloat(0.10);
TFloat init_dt              = TFloat(TFloat(15)/TFloat(100)); // 0.15 TFloat(TFloat(78125)/TFloat(10000000)); //TFloat(0.0078125);; //2^-7
TFloat maxsteps             = TFloat(INT_MAX); 

// massive parallel runs:
int    massive_global;
int    massive_energy; 

// tuning options:
bool   tune                 = false;
double switchfraction       = 0.0001; 
int    tune_mode            = 0; // 0 = always from initial; 1 = from v_best, ...
int    tune_mode_params     = 2; // 0 = alpha..zeta, 1=ODE params, 2=all params
int    N_ITERATIONS         = INT_MAX; 

// Additional ODE heuristics:
bool         heuristics       = false;
bool         apply_heuristics = false; // apply heuristics at the current timestep?

//-------------------------------------------------------------------------------------------------------------------
/* VARIABLES */
bool            debug = false;
void            INThandler(int);       // ctrl-c handler

int * x_rem; // to calculate which vars flipped; not necessary for production

int             * cls;                 //stores the clauses (MAX_LITS_SYSTEM columns)
int             * unit_clause_vars;    //vars which are in unit clauses (one literal)
bool            unit_clauses = false;  //unit clauses existing?
int             * occurrence;          //occurences of each litint             * occurenceCounter;    //counter for allocation
int             * numOccurrenceT;      //number of occurence of a lit
int             * clauseSizes;         //number of literals of a clause
int             * occurenceCounter; 

int             maxNumOccurences = 0;  //max #occurence in clauses of a var
int             n;                     //number of variables
int             m;                     //number of clauses
int             solved;                //is formula solved? UNSAT = 0; SAT = 1;
int             global;                //current lowest loc over all threads
TFloat          global_energy;         //current lowest energy over all threads
TFloat          * v_best;              //x assignment of global local minima for each thread
char            input_filename[256];

/* thread specific vars */
int             * loc_thread;           //current local minima of thread (currently)
TFloat          * energy_thread;        //current energy of all clauses (sum of Cm)
TFloat          * energy_thread_min;    //minimum energy of thread
int             * global_thread;        //global local minima of thread
int             * global_all_runs_thread; //global local minimum of thread over all runs
TFloat          * time_thread;          //simtime of thread (better loc)
TFloat          * time_thread_actual;   //simtime of thread (currently)
double          * walltime_thread;      //walltime of thread (better loc)
TFloat          * initial_assignments;  //initial assignment for each thread
TFloat          * thread_params;        //parameters for each thread
double          * t_begin_thread;       //starting time of thread
double          * t_end_thread;         //current time of thread
int             global_best_thread;     //thread# which has currently best global;

double          * partable;             //if partable.txt is provided, here are the bounds
double          * defaults;             //if partable.txt is provided, here are the default values (max 128 threads)
bool            partable_loaded = false; 

int             * stepcounter;
TFloat          * t_init;               // starting time of integration

/* vars for parameters in CNF: */
char            c_key[256];
char            c_plain[256];
char            c_cipher[256];
int             *c_key_bin;
int             *c_plain_bin;
int             *c_cipher_bin;
bool            c_key_loaded = false;
bool            c_plain_loaded = false;
bool            c_cipher_loaded = false;
int             c_key_best = 0;
int             c_key_best_all_runs = 0;
int             key_best_thread = 0;
int             * c_key_thread;

struct node {
    int id;                 //thread-id
    int *model;             //current assignment
    int *temporal;          //temp assignment for oracle (not used in production)
    int *optimal;           //best assignment and solution afterwards
};

double t_begin;
double t_end; 

//-------------------------------------------------------------------------------------------------------------------
/* color definitions */
#define TEXT_DEFAULT  "\033[0m"
#define TEXT_YELLOW   "\033[1;33m"
#define TEXT_GREEN    "\033[1;32m"
#define TEXT_RED      "\033[1;31m"
#define TEXT_BLUE     "\033[1;34m"
#define TEXT_CYAN     "\033[1;36m"
#define TEXT_WHITE    "\033[1;37m"
#define TEXT_SILVER   "\033[1;315m" 

// test against solution:
int solmax = 0;

//-------------------------------------------------------------------------------------------------------------------
// p-Bit functions
//-------------------------------------------------------------------------------------------------------------------
// generate a random floating point number from min to max
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
//-------------------------------------------------------------------------------------------------------------------
//main p-bit function:
//input: weight (f.e. -4, -2, 0, 0.5, etc.)
//output: -1 / +1

std::mt19937 pbit_generator(seed);
std::uniform_real_distribution<double>  rand_pbit(-1.0, 1.0);

TFloat pbit(double p_in) {
    TFloat p_out = tanh(p_in)+rand_pbit(pbit_generator);
    TFloat res = 0;
    if (p_out>0) res = 1;
    if (p_out<0) res = -1;
    return res; //p_out
}
//-------------------------------------------------------------------------------------------------------------------
// -- / p-Bit
//-------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------------
/* threading */
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

//----------------------------------------------------------------------------------------------------------------
// handler for CTRL-C
void  INThandler(int sig)
{
    signal(sig, SIG_IGN);
    solved=-1; //this stops all the threads
    signal(SIGINT, INThandler);
}
// ---------------------------------------------------------------------------------------------------------------
// manual ode function (NO BOOST) :
// ---------------------------------------------------------------------------------------------------------------
std::vector<TFloat> dmm_dostep(int m_nodeid, std::vector<TFloat> _x, std::vector<TFloat> dxdt, TFloat h) {

    //update V:
    for (int i=0; i<n; i++) {
        _x[i] = _x[i] + h * dxdt[i];
        if (unit_clause_vars[i+1]!=0) _x[i] = unit_clause_vars[i+1]; //TODO: assign +1 or -1
        if (_x[i]<-1.0) _x[i]=-1.0;   
        if (_x[i]>1.0) _x[i]=1.0;
    }
    //update XS:
    for (int i=n; i<n+m; i++) {
        _x[i] = _x[i] + h * dxdt[i];
        if (_x[i]<0.0) _x[i]=0.0;
        if (_x[i]>1.0) _x[i]=1.0;
    }
    //update Xl:
    for (int i=n+m; i<n+m*2; i++) {
        _x[i] = _x[i] + h * dxdt[i];
        if (_x[i]<1.0) _x[i]=1.0;  
        if (_x[i]>xl_max) _x[i]=xl_max;
    }
    
    // increase stepcounter for this thread:
    stepcounter[m_nodeid]++;

    return _x;
}

// ---------------------------------------------------------------------------------------------------------------
// DXDT CALCULATION (WITHOUT BOOST):
// ---------------------------------------------------------------------------------------------------------------
std::vector<TFloat> dmm_generate_dxdt(int m_nodeid, std::vector<TFloat> x , TFloat t ) 
    {
        /* timer and stats */
        t_end_thread[m_nodeid] = pt::microsec_clock::local_time().time_of_day().total_milliseconds();//clock();
        double time_spent = (double)(t_end_thread[m_nodeid] - t_begin_thread[m_nodeid])/1000;// / CLOCKS_PER_SEC;
        time_thread_actual[m_nodeid] = t;

        // stop integration if solved or walltime timeout exceeded: 
        if (solved==1 || time_spent>walltime_timeout || (massive && time_spent < 0) || (tune && time_spent<0)) {
            pthread_mutex_lock(&mut);
                running_threads--;
            pthread_mutex_unlock(&mut);
            pthread_exit(&mut);
        }
        
        energy_thread[m_nodeid] = 0.0;

        std::vector<TFloat> dxdt(n+m*2); // dxdt vector
        std::vector< std::vector <TFloat >> dxdt_v(n); // vector for storing all voltages
        
        /* main ODE routine */
        if (solved!=1) {
            int loc;
            // Loop ---------------------------------------------------------------------------------------------------------------
            loc = m;
            // screen info variables:
            TFloat C_rem, R_rem, G_rem;
            int cnt_xl_pos = 0;
            int cnt_xs_pos = 0;
            // loop through each clause:
            for (int clause = 0; clause < m; clause++) {
                // number of literals in this clause:
                int ksat = clauseSizes[clause];
                // Xl & Xs:
                TFloat Xs = x[clause+n];   if (Xs<0.0) Xs = TFloat(0.0); if (Xs>1.0) Xs = TFloat(1.0); //Xs bounds
                TFloat Xl = x[clause+n+m]; if (Xl<1.0) Xl = TFloat(1.0); if (Xl>xl_max) Xl = TFloat(xl_max); //Xl bounds
                TFloat C  = TFloat(0.0);
                TFloat Ri, Rj, Rk, Gi, Gj, Gk;
                // 3-sat:
                if (ksat==3) {
                    int Qi = (cls[clause*MAX_LITS_SYSTEM+0]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qj = (cls[clause*MAX_LITS_SYSTEM+1]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qk = (cls[clause*MAX_LITS_SYSTEM+2]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int liti = abs(cls[clause*MAX_LITS_SYSTEM+0]);
                    int litj = abs(cls[clause*MAX_LITS_SYSTEM+1]);
                    int litk = abs(cls[clause*MAX_LITS_SYSTEM+2]);
                    TFloat Vi = x[liti-1]; if (Vi<-1.0) Vi = -1.0; if (Vi>1.0) Vi = 1.0; //V bounds
                    TFloat Vj = x[litj-1]; if (Vj<-1.0) Vj = -1.0; if (Vj>1.0) Vj = 1.0; //V bounds
                    TFloat Vk = x[litk-1]; if (Vk<-1.0) Vk = -1.0; if (Vk>1.0) Vk = 1.0; //V bounds
                    TFloat i = TFloat(1.0)-TFloat(Qi)*Vi;
                    TFloat j = TFloat(1.0)-TFloat(Qj)*Vj;
                    TFloat k = TFloat(1.0)-TFloat(Qk)*Vk;
                    C = TFloat(fmin(i, fmin(j, k)));
                    C = C / TFloat(2.0);
                    if (C<0.0) C=TFloat(0.0);
                    if (C>1.0) C=TFloat(1.0);
                    
                    // equation Gn,m(vn,vj,vk)= 1/2 qn,mmin[(1−qj,mvj),(1−qk,mvk)] (5.x):
                    Gi = TFloat(Qi) * fmin(j,k) / TFloat(2.0);
                    Gj = TFloat(Qj) * fmin(i,k) / TFloat(2.0);
                    Gk = TFloat(Qk) * fmin(i,j) / TFloat(2.0);
                    
                    // equation Rn,m (vn , vj , vk ) = 1/2(qn,m −vn), Cm(vn,vj,vk)= 1/2(1−qn,mvn), 0 otherwise (5.x):
                    if (C == TFloat(i/TFloat(2.0)) ) {Ri = (TFloat(Qi)-Vi)/2.0;} else {Ri = TFloat(0.0);} //Qi*i/2.0*-1;} //= 0.0
                    if (C == TFloat(j/TFloat(2.0)) ) {Rj = (TFloat(Qj)-Vj)/2.0;} else {Rj = TFloat(0.0);} //Qj*j/2.0*-1;} //= 0.0
                    if (C == TFloat(k/TFloat(2.0)) ) {Rk = (TFloat(Qk)-Vk)/2.0;} else {Rk = TFloat(0.0);} //Qk*k/2.0*-1;} //= 0.0
                    
                    // equation Vn = SUM xl,mxs,mGn,m + (1 + ζxl,m)(1 − xs,m)Rn,m (5.x):
                    TFloat _Vi = Xl * Xs * Gi + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * Ri ;
                    TFloat _Vj = Xl * Xs * Gj + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * Rj ;
                    TFloat _Vk = Xl * Xs * Gk + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * Rk ;

                    /*dxdt[liti-1] = dxdt[liti-1] + _Vi;
                    dxdt[litj-1] = dxdt[litj-1] + _Vj;
                    dxdt[litk-1] = dxdt[litk-1] + _Vk; */
                    //sum of vectors method:
                    if (_Vi!=0.0) dxdt_v[liti-1].push_back(_Vi);
                    if (_Vj!=0.0) dxdt_v[litj-1].push_back(_Vj);
                    if (_Vk!=0.0) dxdt_v[litk-1].push_back(_Vk);

                    // do not change unit_clauses:
                    if (unit_clauses) {
                        if (unit_clause_vars[liti]!=0) dxdt[liti-1] = 0.0;
                        if (unit_clause_vars[litj]!=0) dxdt[litj-1] = 0.0;
                        if (unit_clause_vars[litk]!=0) dxdt[litk-1] = 0.0;
                    }
                }
                // 2-sat:
                if (ksat==2) {
                    int Qi = (cls[clause*MAX_LITS_SYSTEM+0]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qj = (cls[clause*MAX_LITS_SYSTEM+1]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int liti = abs(cls[clause*MAX_LITS_SYSTEM+0]);
                    int litj = abs(cls[clause*MAX_LITS_SYSTEM+1]);
                    TFloat Vi = x[liti-1]; if (Vi<-1.0) Vi = -1.0; if (Vi>1.0) Vi = 1.0; 
                    TFloat Vj = x[litj-1]; if (Vj<-1.0) Vj = -1.0; if (Vj>1.0) Vj = 1.0;
                    TFloat i = TFloat(1.0)-TFloat(Qi)*Vi;
                    TFloat j = TFloat(1.0)-TFloat(Qj)*Vj;
                    C = TFloat(fmin(i, j));
                    C = C / TFloat(2.0) ;
                    if (C<0.0) C=TFloat(0.0);
                    if (C>1.0) C=TFloat(1.0);
                    //voltage:
                    Gi = TFloat(Qi) * j / TFloat(2.0);
                    Gj = TFloat(Qj) * i / TFloat(2.0);
                    
                    if (C == TFloat(i/TFloat(2.0)) ) {Ri = (TFloat(Qi)-Vi)/2.0;} else {Ri = TFloat(0.0);} //Qi*i/2.0*-1;} //= 0.0
                    if (C == TFloat(j/TFloat(2.0)) ) {Rj = (TFloat(Qj)-Vj)/2.0;} else {Rj = TFloat(0.0);} //Qj*j/2.0*-1;} //= 0.0

                    TFloat _Vi = Xl * Xs * Gi + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * Ri;
                    TFloat _Vj = Xl * Xs * Gj + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * Rj;

                    /*dxdt[liti-1] = dxdt[liti-1] + _Vi;
                    dxdt[litj-1] = dxdt[litj-1] + _Vj;*/
                    //sum of vectors method:
                    if (_Vi!=0.0) dxdt_v[liti-1].push_back(_Vi);
                    if (_Vj!=0.0) dxdt_v[litj-1].push_back(_Vj);
                    
                    // do not change unit_clauses:
                    if (unit_clauses) {
                        if (unit_clause_vars[liti]!=0) dxdt[liti-1] = 0.0;
                        if (unit_clause_vars[litj]!=0) dxdt[litj-1] = 0.0;
                    }
                }
                
                // k-sat:
                if (ksat!=1 && ksat!=2 && ksat!=3) {
                    int lit[MAX_LITS_SYSTEM], Q[MAX_LITS_SYSTEM];
                    TFloat V[MAX_LITS_SYSTEM], _i[MAX_LITS_SYSTEM], R[MAX_LITS_SYSTEM], G[MAX_LITS_SYSTEM];

                    TFloat c_min=INT_MAX;
                    for (int i=0; i<ksat; i++) {
                        Q[i] = (cls[clause*MAX_LITS_SYSTEM+i]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                        lit[i] = abs(cls[clause*MAX_LITS_SYSTEM+i]);
                        V[i] = x[lit[i]-1]; if (V[i]<-1.0) V[i]=-1.0; if (V[i]>1.0) V[i]=1.0; //boundary for v € [-1,1]:
                        _i[i] = TFloat(1.0)-TFloat(Q[i])*V[i];
                        // find min:
                        if (_i[i]<c_min) c_min = _i[i]; 
                    }
                    C = c_min / TFloat(2.0);
                    if (C<0.0) printf("*\n");//C=0.0; // never triggered?
                    if (C>1.0) printf("*\n");//C=1.0; // never triggered?
                    
                    for (int i=0; i<ksat; i++) {
                        //find min of others:
                        TFloat g_min = INT_MAX;
                        for (int ii=0; ii<ksat; ii++) {if (ii!=i && _i[ii]<g_min) g_min = _i[ii];}
                        G[i] = TFloat(Q[i]) * g_min / TFloat(2.0);
                        TFloat comp = _i[i]/TFloat(2.0);
                        if (comp<0.0) printf("*\n");//comp = 0.0; // never triggered?
                        if (comp>1.0) printf("*\n");//comp = 1.0; // never triggered?
                        if (C != comp) {R[i] = TFloat(0.0);} else {R[i] = (TFloat(Q[i])-V[i]) / TFloat(2.0);}
                        TFloat _V = Xl * Xs * G[i] + (TFloat(1.0) + dmm_zeta * Xl) * (TFloat(1.0) - Xs) * R[i];
                        
                        //dxdt[lit[i]-1] = dxdt[lit[i]-1] + _V;
                        //sum of vectors method:
                        if (_V!=0.0) dxdt_v[lit[i]-1].push_back(_V);
                    
                        // do not change unit_clauses:
                        if (unit_clauses) {
                            if (unit_clause_vars[lit[i]]!=0) dxdt[lit[i]-1] = 0.0;
                        }
                    }
                }

                //update energy:
                energy_thread[m_nodeid] += C;
                //update loc:
                if (C<0.5) loc--; //this clause is sat, reduce loc
                //if (C==0.0) loc--;
                
                // Calculate new Xs:
                dxdt[n+clause] = dmm_beta * (Xs + dmm_epsilon) * (C - dmm_gamma);
                
                // Calculate new Xl:
                dxdt[n+m+clause] = dmm_alpha_cm[m_nodeid*m+clause] * (C - dmm_delta);

                // update info variables:
                if (clause==0) {
                    C_rem = C; 
                    if (ksat<=3) {
                        G_rem = Gi; R_rem = Ri;
                    } else {
                        G_rem = 0.0;
                        R_rem = 0.0;
                    }
                    if (ksat==1) {
                        G_rem = 0.0; 
                        R_rem = 0.0;
                    }
                }
                if (x[n+m+clause]>1.0) cnt_xl_pos++;
                if (x[n+clause]>0.0) cnt_xs_pos++;
            } //---clause calculation loop

            // summation of voltages SUM dxdt_v[n] => dxdt[n]: ------------------------------------------------------
            for (int i=0; i<n; i++) {
                std::sort(dxdt_v[i].begin(), dxdt_v[i].end()); //summing with smallest first increases accuracy
                dxdt[i] = accumulate(dxdt_v[i].begin(), dxdt_v[i].end(), (TFloat) 0.0);
            }
            // -------------------------------------------------------------------------------------------------------

            // test against solution:
            int corr = 0;
            
            //update energy:
            if (tune) energy_thread[m_nodeid] = loc; //tuning on loc...
              
            // screen output: ----------------------------------------------------------------------------------------
            if (!c_key_loaded && !quiet && m_nodeid==global_best_thread) {
                int showglobal = global_thread[m_nodeid];
                if (massive) showglobal = massive_global;
                std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT
                << "\rc [" << m_nodeid << "] " << time_spent << "s "
                << "T=" << t
                << " GLOBAL=" << global
                << " (" << showglobal << ")" 
                << " (LOC=" << loc << ")" 
                << " (" << stepcounter[m_nodeid] << ")"
                //<< " α=" << m_alpha
                << " α=" << dmm_alpha_cm[m_nodeid*m+0]
                << " β=" << dmm_beta
                << " C=" << C_rem
                << " (" << cls[0] << "=" << x[abs(cls[0])-1]
                << ", " << cls[1] << "=" << x[abs(cls[1])-1]
                << ", " << cls[2] << "=" << x[abs(cls[2])-1] << ")"
                << " Xs=" << x[n]
                << " Xl=" << x[n+m] 
                //<< " R=" << R_rem 
                //<< " G=" << G_rem 
                << " #Xs>0:" << cnt_xs_pos 
                << " #Xl>1:" << cnt_xl_pos 
                << " Σe=" << energy_thread[m_nodeid] << " "
                // test against solution:
                << (n-corr) << "/" << n << " (best: " << (n-solmax) << ") ";
                fflush(stdout);
                
            }

            // update global_all_runs_thread & update v_best: --------------------------------------------------------
            if (loc <= global_all_runs_thread[m_nodeid]) {
                global_all_runs_thread[m_nodeid] = loc;
                t_init[m_nodeid] = t;
                pthread_mutex_lock(&mut);
                    //update v_best array:
                    for (int i=0; i<n+m*2; i++) v_best[m_nodeid*(n+m*2)+i] = x[i]; 
                pthread_mutex_unlock(&mut);
            }
            // update loc of thread: ---------------------------------------------------------------------------------
            loc_thread[m_nodeid] = loc;
            time_thread_actual[m_nodeid] = t;

            //new lower lock (global)? or lower energy (global)? -----------------------------------------------------
            if (loc<global || energy_thread[m_nodeid]<global_energy || energy_thread[m_nodeid]<energy_thread_min[m_nodeid] || loc<global_thread[m_nodeid]) {
                if (loc<global) {
                    pthread_mutex_lock(&mut);
                        global = loc;
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                    //if (!tune) t_begin_thread[m_nodeid] = t_end_thread[m_nodeid];
                }
                if (energy_thread[m_nodeid]<global_energy) {
                    pthread_mutex_lock(&mut);
                        global_energy = energy_thread[m_nodeid];
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                }
                if (energy_thread[m_nodeid]<energy_thread_min[m_nodeid]) {
                    pthread_mutex_lock(&mut);
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                }
                if (loc<global_thread[m_nodeid]) {
                    pthread_mutex_lock(&mut);
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                }

                if (!quiet) {
                    int showglobal = global_thread[m_nodeid];
                    if (massive) showglobal = massive_global;
                    std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT
                    << "\rc [" << m_nodeid << "] " << time_spent << "s "
                    << "T=" << t
                    << " GLOBAL=" << global
                    << " (" << showglobal << ")"
                    << " (LOC=" << loc << ")" 
                    << " (" << stepcounter[m_nodeid] << ")"
                    //<< " α=" << m_alpha
                    << " α=" << dmm_alpha_cm[m_nodeid*m+0]
                    << " β=" << dmm_beta
                    << " C=" << C_rem
                    << " (" << cls[0] << "=" << x[abs(cls[0])-1]
                    << ", " << cls[1] << "=" << x[abs(cls[1])-1]
                    << ", " << cls[2] << "=" << x[abs(cls[2])-1] << ")"
                    << " Xs=" << x[n]
                    << " Xl=" << x[n+m] 
                    //<< " R=" << R_rem 
                    //<< " G=" << G_rem 
                    << " #Xs>0:" << cnt_xs_pos 
                    << " #Xl>1:" << cnt_xl_pos 
                    << " Σe=" << energy_thread[m_nodeid] << " "
                    // test against solution:
                    << (n-corr) << "/" << n << " (best: " << (n-solmax) << ") " << std::endl;
                    fflush(stdout);
                } else {std::cout << std::endl;}
            }

            // update energy of thread: ------------------------------------------------------------------------------
            if (energy_thread[m_nodeid] < energy_thread_min[m_nodeid]) {
                energy_thread_min[m_nodeid] = energy_thread[m_nodeid];
            }

            //new thread global? --------------------------------------------------------------------------------------
            if (loc<global_thread[m_nodeid]) {
                // update global_thread, time_thread and walltime_thread:
                global_thread[m_nodeid] = loc;
                time_thread[m_nodeid] = t;
                walltime_thread[m_nodeid] = time_spent;
            }
            // loc file? (used for dashboard) --------------------------------------------------------------------------
            if (writelocfile && m_nodeid==global_best_thread && stepcounter[global_best_thread] % 100 == 0) {
                #ifdef USEFPRINTF
                    FILE *floc = fopen(LOC_FILE, "a");
                    fprintf(floc,"%d;%.5f;%.5f;%d;%.2f\n",global_best_thread,time_spent,t,global,global_energy);
                    fclose(floc);
                #endif
                
            }
            //solved? -------------------------------------------------------------------------------------------------
            //if (energy_thread[m_nodeid]<0.1) {
            if (loc==0) {
                quiet=1;
                pthread_mutex_lock(&mut);
                    printf(TEXT_YELLOW);
                    printf("\nc [%d] T=",m_nodeid); std::cout << t << " SOLUTION FOUND" << std::endl; 
                    for (int i=0; i<n; i++) {
                        if (x[i]>0) printf("%d ",(i+1));
                        if (x[i]<0) printf("%d ",(i+1)*-1);
                    }
                    fflush(stdout);
                    printf("\nc [%d] VERIFYING...\n",m_nodeid); if (!quiet) printf(TEXT_DEFAULT);
                    bool sat = true; bool clausesat;
                    for (int i=0; i<m; i++) {
                        for (int j=0; j<clauseSizes[i]; j++) {
                            clausesat = false;
                            int lit = abs(cls[i*MAX_LITS_SYSTEM+j]);
                            if ( (x[lit-1]>0 && cls[i*MAX_LITS_SYSTEM+j]>0) || (x[lit-1]<0 && cls[i*MAX_LITS_SYSTEM+j]<0) ) {
                                clausesat = true;
                                break;
                            }
                        }
                        if (!clausesat) {
                            sat = false;
                            break;
                        }
                    }
                    printf(TEXT_YELLOW);
                    if (sat)  {
                        printf(TEXT_YELLOW); printf("c [%d] SAT (VERIFIED)\n",m_nodeid); solved = 1;
                        //write solution to file:
                        #ifdef USEFPRINTF
                            FILE *fs = fopen(SOLUTION_FILE, "w");
                            for (int i=0; i<n; i++) {
                                if (x[i]>=0) fprintf(fs,"%d, ",i+1);
                                if (x[i]<0) fprintf(fs,"%d, ",(i+1)*-1);
                                //fprintf(fs,"%.5f, ",x[i]); // current G_field = solution
                            }
                            fclose(fs);
                        #endif
                        
                    }
                    if (!sat) {printf(TEXT_RED); printf("c [%d] UNSAT (VERIFIED)\n",m_nodeid);}
                    printf(TEXT_DEFAULT);
                pthread_mutex_unlock(&mut);

                // update locfile:
                if (writelocfile) {
                    #ifdef USEFPRINTF
                        FILE *floc = fopen(LOC_FILE, "a");
                        fprintf(floc,"%d;%.5f;%.5f;%d;%.2f\n",global_best_thread,time_spent,t,global,global_energy);
                        fclose(floc);
                    #endif
                    
                }
                quiet=0;
            } // ---output
        } 
        return dxdt;
}


//-------------------------------------------------------------------------------------------------------------------
// ode function as structure: (boost version)
struct odeS
{
    int    m_nodeid;
    TFloat m_alpha; // unused; we are using dmm_alpha_cm[], one alpha for each clause
    TFloat m_beta;
    TFloat m_gamma;
    TFloat m_delta;
    TFloat m_epsilon;
    TFloat m_zeta;

    odeS(int t_nodeid=0, TFloat t_alpha=5.0, TFloat t_beta=20.0, TFloat t_gamma=0.23, TFloat t_delta=0.04, TFloat t_epsilon=0.11, TFloat t_zeta=0.11)
    : m_nodeid(t_nodeid), m_alpha(t_alpha), m_beta(t_beta), m_gamma(t_gamma), m_delta(t_delta), m_epsilon(t_epsilon), m_zeta(t_zeta) {}

    void operator()(const state_type &x , state_type &dxdt , TFloat t ) const
    {
        /* timer and stats */
        t_end_thread[m_nodeid] = pt::microsec_clock::local_time().time_of_day().total_milliseconds();//clock();
        double time_spent = (double)(t_end_thread[m_nodeid] - t_begin_thread[m_nodeid])/1000;// / CLOCKS_PER_SEC;
        time_thread_actual[m_nodeid] = t;

        // stop integration if solved or walltime timeout exceeded: 
        if (solved==1 || time_spent>walltime_timeout || (tune && time_spent<0)) {
            pthread_mutex_lock(&mut);
                running_threads--;
            pthread_mutex_unlock(&mut);
            pthread_exit(&mut);
        }
        
        energy_thread[m_nodeid] = 0.0;

        /* main ODE routine */
        /* This bounding procedure effectively mimics the numerical technique that we use to bound the dynamics, where any outward pointing component of the flow field on the boundary is simply ignored.*/
        if (solved!=1) {
            int loc;
            // Loop ---------------------------------------------------------------------------------------------------------------
            loc = m;
            //set all V to 0.00:
            //fill( dxdt.begin() , dxdt.begin()+n , TFloat(0.0) ); //not required?
            // screen info variables:
            TFloat C_rem, R_rem, G_rem;
            int cnt_xl_pos = 0;
            int cnt_xs_pos = 0;
            // loop through each clause:
            for (int clause = 0; clause < m; clause++) {
            //for (int clause = m-1; clause >= 0; clause--) {
                // number of literals in this clause:
                int ksat = clauseSizes[clause];
                // Xl & Xs:
                TFloat Xs = x[clause+n];   if (Xs<0.0) Xs = TFloat(0.0); if (Xs>1.0) Xs = TFloat(1.0); //Xs bounds
                TFloat Xl = x[clause+n+m]; if (Xl<1.0) Xl = TFloat(1.0); if (Xl>xl_max) Xl = TFloat(xl_max); //Xl bounds
                TFloat C  = TFloat(0.0);
                TFloat Ri, Rj, Rk, Gi, Gj, Gk;
                // 3-sat:
                if (ksat==3) {
                    int Qi = (cls[clause*MAX_LITS_SYSTEM+0]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qj = (cls[clause*MAX_LITS_SYSTEM+1]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qk = (cls[clause*MAX_LITS_SYSTEM+2]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int liti = abs(cls[clause*MAX_LITS_SYSTEM+0]);
                    int litj = abs(cls[clause*MAX_LITS_SYSTEM+1]);
                    int litk = abs(cls[clause*MAX_LITS_SYSTEM+2]);
                    TFloat Vi = x[liti-1]; if (Vi<-1.0) Vi = -1.0; if (Vi>1.0) Vi = 1.0; //V bounds
                    TFloat Vj = x[litj-1]; if (Vj<-1.0) Vj = -1.0; if (Vj>1.0) Vj = 1.0; //V bounds
                    TFloat Vk = x[litk-1]; if (Vk<-1.0) Vk = -1.0; if (Vk>1.0) Vk = 1.0; //V bounds
                    TFloat i = TFloat(1.0)-TFloat(Qi)*Vi;
                    TFloat j = TFloat(1.0)-TFloat(Qj)*Vj;
                    TFloat k = TFloat(1.0)-TFloat(Qk)*Vk;
                    C = TFloat(fmin(i, fmin(j, k)));
                    C = C / TFloat(2.0);
                    if (C<0.0) C=TFloat(0.0);
                    if (C>1.0) C=TFloat(1.0);
                    
                    // equation Gn,m(vn,vj,vk)= 1/2 qn,mmin[(1−qj,mvj),(1−qk,mvk)] (5.x):
                    Gi = TFloat(Qi) * fmin(j,k) / TFloat(2.0);
                    Gj = TFloat(Qj) * fmin(i,k) / TFloat(2.0);
                    Gk = TFloat(Qk) * fmin(i,j) / TFloat(2.0);

                    // equation Rn,m (vn , vj , vk ) = 1/2(qn,m −vn), Cm(vn,vj,vk)= 12(1−qn,mvn), 0 otherwise (5.x):
                    if (C != i/TFloat(2.0) ) {Ri = TFloat(0.0);} else {Ri = (TFloat(Qi) - Vi) / TFloat(2.0);}
                    if (C != j/TFloat(2.0) ) {Rj = TFloat(0.0);} else {Rj = (TFloat(Qj) - Vj) / TFloat(2.0);}
                    if (C != k/TFloat(2.0) ) {Rk = TFloat(0.0);} else {Rk = (TFloat(Qk) - Vk) / TFloat(2.0);}
                    

                    // equation Vn = SUM xl,mxs,mGn,m + (1 + ζxl,m)(1 − xs,m)Rn,m (5.x):
                    TFloat _Vi = Xl * Xs * Gi + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * Ri ;
                    TFloat _Vj = Xl * Xs * Gj + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * Rj ;
                    TFloat _Vk = Xl * Xs * Gk + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * Rk ;

                    dxdt[liti-1] = dxdt[liti-1] + _Vi;
                    dxdt[litj-1] = dxdt[litj-1] + _Vj;
                    dxdt[litk-1] = dxdt[litk-1] + _Vk;  

                    // do not change unit_clauses:
                    if (unit_clauses) {
                        if (unit_clause_vars[liti]!=0) dxdt[liti-1] = 0.0;
                        if (unit_clause_vars[litj]!=0) dxdt[litj-1] = 0.0;
                        if (unit_clause_vars[litk]!=0) dxdt[litk-1] = 0.0;
                    }
                }
                // 2-sat:
                if (ksat==2) {
                    int Qi = (cls[clause*MAX_LITS_SYSTEM+0]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int Qj = (cls[clause*MAX_LITS_SYSTEM+1]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                    int liti = abs(cls[clause*MAX_LITS_SYSTEM+0]);
                    int litj = abs(cls[clause*MAX_LITS_SYSTEM+1]);
                    TFloat Vi = x[liti-1]; if (Vi<-1.0) Vi = -1.0; if (Vi>1.0) Vi = 1.0; 
                    TFloat Vj = x[litj-1]; if (Vj<-1.0) Vj = -1.0; if (Vj>1.0) Vj = 1.0;
                    TFloat i = TFloat(1.0)-TFloat(Qi)*Vi;
                    TFloat j = TFloat(1.0)-TFloat(Qj)*Vj;
                    C = TFloat(fmin(i, j));
                    C = C / TFloat(2.0) ;
                    if (C<0.0) C=TFloat(0.0);
                    if (C>1.0) C=TFloat(1.0);
                    //voltage:
                    Gi = TFloat(Qi) * j / TFloat(2.0);
                    Gj = TFloat(Qj) * i / TFloat(2.0);
                    if (C != i/TFloat(2.0) ) {Ri = TFloat(0.0);} else {Ri = (TFloat(Qi) - Vi) / TFloat(2.0);} 
                    if (C != j/TFloat(2.0) ) {Rj = TFloat(0.0);} else {Rj = (TFloat(Qj) - Vj) / TFloat(2.0);}
                    
                    TFloat _Vi = Xl * Xs * Gi + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * Ri;
                    TFloat _Vj = Xl * Xs * Gj + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * Rj;

                    dxdt[liti-1] = dxdt[liti-1] + _Vi;
                    dxdt[litj-1] = dxdt[litj-1] + _Vj;

                    // do not change unit_clauses:
                    if (unit_clauses) {
                        if (unit_clause_vars[liti]!=0) dxdt[liti-1] = 0.0;
                        if (unit_clause_vars[litj]!=0) dxdt[litj-1] = 0.0;
                    }
                }
                
                // k-sat:
                if (ksat!=1 && ksat!=2 && ksat!=3) {
                    int lit[MAX_LITS_SYSTEM], Q[MAX_LITS_SYSTEM];
                    TFloat V[MAX_LITS_SYSTEM], _i[MAX_LITS_SYSTEM], R[MAX_LITS_SYSTEM], G[MAX_LITS_SYSTEM];

                    TFloat c_min=INT_MAX;
                    for (int i=0; i<ksat; i++) {
                        Q[i] = (cls[clause*MAX_LITS_SYSTEM+i]>0)? 1:-1; // +1 if literal is >0, otherwise -1
                        lit[i] = abs(cls[clause*MAX_LITS_SYSTEM+i]);
                        V[i] = x[lit[i]-1]; if (V[i]<-1.0) V[i]=-1.0; if (V[i]>1.0) V[i]=1.0; //boundary for v € [-1,1]:
                        _i[i] = 1.0-Q[i]*V[i];
                        // find min:
                        if (_i[i]<c_min) c_min = _i[i]; 
                    }
                    C = c_min / TFloat(2.0);
                    if (C<0.0) printf("*\n");//C=0.0; // never triggered?
                    if (C>1.0) printf("*\n");//C=1.0; // never triggered?
                    
                    for (int i=0; i<ksat; i++) {
                        //find min of others:
                        TFloat g_min = INT_MAX;
                        for (int ii=0; ii<ksat; ii++) {if (ii!=i && _i[ii]<g_min) g_min = _i[ii];}
                        G[i] = Q[i] * g_min / TFloat(2.0);
                        TFloat comp = (1.0-Q[i]*V[i])/TFloat(2.0);
                        if (comp<0.0) printf("*\n");//comp = 0.0; // never triggered?
                        if (comp>1.0) printf("*\n");//comp = 1.0; // never triggered?
                        if (C != comp) {R[i] = 0.0;} else {R[i] = (Q[i] - V[i]) / TFloat(2.0);}
                        TFloat _V = Xl * Xs * G[i] + (TFloat(1.0) + m_zeta * Xl) * (TFloat(1.0) - Xs) * R[i];
                        dxdt[lit[i]-1] = dxdt[lit[i]-1] + _V;

                        // do not change unit_clauses:
                        if (unit_clauses) {
                            if (unit_clause_vars[lit[i]]!=0) dxdt[lit[i]-1] = 0.0;
                        }
                    }
                }
                //update energy:
                energy_thread[m_nodeid] += C;
                //update loc:
                if (C<0.5) loc--; //this clause is sat, reduce loc
                //if (C==0.0) loc--;
                
                // Calculate new Xs:
                dxdt[n+clause] = m_beta * (Xs + m_epsilon) * (C - m_gamma);
                
                // Calculate new Xl:
                dxdt[n+m+clause] = dmm_alpha_cm[m_nodeid*m+clause] * (C - m_delta);

                // update info variables:
                if (clause==0) {
                    C_rem = C; 
                    if (ksat<=3) {
                        G_rem = Gi; R_rem = Ri;
                    } else {
                        G_rem = 0.0;
                        R_rem = 0.0;
                    }
                    if (ksat==1) {
                        G_rem = 0.0; 
                        R_rem = 0.0;
                    }
                }
                if (x[n+m+clause]>1.0) cnt_xl_pos++;
                if (x[n+clause]>0.0) cnt_xs_pos++;
            }

            // test against solution:
            int corr = 0;
        
            //update energy:
            if (tune) energy_thread[m_nodeid] = loc;
              
            // screen output: ----------------------------------------------------------------------------------------
            if (!c_key_loaded && !quiet && m_nodeid==global_best_thread) {
                std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT
                << "\rc [" << m_nodeid << "] " << time_spent << "s "
                << "T=" << t
                << " GLOBAL=" << global
                << " (THREAD=" << global_thread[m_nodeid] << ")" 
                << " (LOC=" << loc << ")" 
                << " (" << stepcounter[m_nodeid] << ")"
                //<< " α=" << m_alpha
                << " α=" << dmm_alpha_cm[m_nodeid*m+0]
                << " β=" << m_beta
                << " C=" << C_rem
                << " (" << cls[0] << "=" << x[abs(cls[0])-1]
                << ", " << cls[1] << "=" << x[abs(cls[1])-1]
                << ", " << cls[2] << "=" << x[abs(cls[2])-1] << ")"
                << " Xs=" << x[n]
                << " Xl=" << x[n+m] 
                //<< " R=" << R_rem 
                //<< " G=" << G_rem 
                << " #Xs>0:" << cnt_xs_pos 
                << " #Xl>1:" << cnt_xl_pos 
                << " Σe=" << energy_thread[m_nodeid] << "  ";
                fflush(stdout);
                // test against solution:
                std::cout << (n-corr) << "/" << n << " (best: " << (n-solmax) << ") ";
                
            }

            //flowvector output: -------------------------------------------------------------------------------------
            /*
            if (m_nodeid==0 && flowvector_log) {
                FILE *fvec = fopen(FLOWVECTOR_FILE, "a");
                fprintf(fvec,"%.15f;",t);
                for (int i=0; i<n+m*2;i++) fprintf(fvec,"%.15f;",x[i]);
                fprintf(fvec,"\n");
                fclose(fvec);
            }
            */
            
            // update global_all_runs_thread & update v_best: --------------------------------------------------------
            if (loc <= global_all_runs_thread[m_nodeid]) {
                global_all_runs_thread[m_nodeid] = loc;
                t_init[m_nodeid] = t;
                pthread_mutex_lock(&mut);
                    //update v_best array:
                    for (int i=0; i<n+m*2; i++) v_best[m_nodeid*(n+m*2)+i] = x[i]; //G_field[i];
                pthread_mutex_unlock(&mut);
            }
            // update loc of thread: ---------------------------------------------------------------------------------
            loc_thread[m_nodeid] = loc;
            time_thread_actual[m_nodeid] = t;

            // update energy of thread: ------------------------------------------------------------------------------
            if (energy_thread[m_nodeid] < energy_thread_min[m_nodeid]) {
                energy_thread_min[m_nodeid] = energy_thread[m_nodeid];
            }
            
            //new lower lock (global)? or lower energy (global)? -----------------------------------------------------
            if (loc<global || energy_thread[m_nodeid]<global_energy) {
                if (loc<global) {
                    pthread_mutex_lock(&mut);
                        global = loc;
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                }
                if (energy_thread[m_nodeid]<global_energy) {
                    pthread_mutex_lock(&mut);
                        global_energy = energy_thread[m_nodeid];
                        global_best_thread = m_nodeid;
                    pthread_mutex_unlock(&mut);
                }

                if (!quiet) {
                    std::cout << std::setprecision(2) << std::fixed 
                << "\rc [" << m_nodeid << "] " << time_spent << "s " 
                << "T=" << t
                << " GLOBAL=" << global
                << " (THREAD=" << global_thread[m_nodeid] << ")"
                << " (LOC=" << loc << ")" 
                << " (" << stepcounter[m_nodeid] << ")"
                //<< " α=" << m_alpha
                << " α=" << dmm_alpha_cm[m_nodeid*m+0]
                << " β=" << m_beta
                << " C=" << C_rem
                << " (" << cls[0] << "=" << x[abs(cls[0])-1]
                << ", " << cls[1] << "=" << x[abs(cls[1])-1]
                << ", " << cls[2] << "=" << x[abs(cls[2])-1] << ")"
                << " Xs=" << x[n]
                << " Xl=" << x[n+m] 
                //<< " R=" << R_rem 
                //<< " G=" << G_rem 
                << " #Xs>0:" << cnt_xs_pos 
                << " #Xl>1:" << cnt_xl_pos 
                << " Σe=" << energy_thread[m_nodeid] << TEXT_DEFAULT << std::endl;
                fflush(stdout);
                // test against solution:
                std::cout << (n-corr) << "/" << n << " (best: " << (n-solmax) << ") ";
                }
                // reset timer and stepcounters:
                t_begin_thread[m_nodeid] = pt::microsec_clock::local_time().time_of_day().total_milliseconds();
            }
            //new thread global? --------------------------------------------------------------------------------------
            if (loc<global_thread[m_nodeid]) {
                // update global_thread, time_thread and walltime_thread:
                global_thread[m_nodeid] = loc;
                time_thread[m_nodeid] = t;
                walltime_thread[m_nodeid] = time_spent;
            }
            // loc file? (used for dashboard) --------------------------------------------------------------------------
            if (writelocfile && m_nodeid==global_best_thread && stepcounter[global_best_thread] % 100 == 0) {
                #ifdef USEFPRINTF
                    FILE *floc = fopen(LOC_FILE, "a");
                    fprintf(floc,"%d;%.5f;%.5f;%d;%.2f\n",global_best_thread,time_spent,t,global,global_energy);
                    fclose(floc);
                #endif
                
            }
            //solved? -------------------------------------------------------------------------------------------------
            //if (energy_thread[m_nodeid]<0.1) {
            if (loc==0) {
                quiet=1;
                pthread_mutex_lock(&mut);
                    printf(TEXT_YELLOW);
                    printf("\nc [%d] T=",m_nodeid); std::cout << t << " SOLUTION FOUND" << std::endl; 
                    for (int i=0; i<n; i++) {
                        if (x[i]>0) printf("%d ",(i+1));
                        if (x[i]<0) printf("%d ",(i+1)*-1);
                    }
                    fflush(stdout);
                    printf("\nc [%d] VERIFYING...\n",m_nodeid); if (!quiet) printf(TEXT_DEFAULT);
                    bool sat = true; bool clausesat;
                    for (int i=0; i<m; i++) {
                        for (int j=0; j<clauseSizes[i]; j++) {
                            clausesat = false;
                            int lit = abs(cls[i*MAX_LITS_SYSTEM+j]);
                            if ( (x[lit-1]>0 && cls[i*MAX_LITS_SYSTEM+j]>0) || (x[lit-1]<0 && cls[i*MAX_LITS_SYSTEM+j]<0) ) {
                                clausesat = true;
                                break;
                            }
                        }
                        if (!clausesat) {
                            sat = false;
                            break;
                        }
                    }
                    printf(TEXT_YELLOW);
                    if (sat)  {
                        printf(TEXT_YELLOW); printf("c [%d] SAT (VERIFIED)\n",m_nodeid); solved = 1;
                        //write solution to file:
                        #ifdef USEFPRINTF
                            FILE *fs = fopen(SOLUTION_FILE, "w");
                            for (int i=0; i<n; i++) {
                                if (x[i]>=0) fprintf(fs,"%d, ",i+1);
                                if (x[i]<0) fprintf(fs,"%d, ",(i+1)*-1);
                                //fprintf(fs,"%.5f, ",x[i]); // current G_field = solution
                            }
                            fclose(fs);
                        #endif
                        
                    }
                    if (!sat) {printf(TEXT_RED); printf("c [%d] UNSAT (VERIFIED)\n",m_nodeid);}
                    printf(TEXT_DEFAULT);
                pthread_mutex_unlock(&mut);

                // update locfile:
                if (writelocfile) {
                    #ifdef USEFPRINTF
                        FILE *floc = fopen(LOC_FILE, "a");
                        fprintf(floc,"%d;%.5f;%.5f;%d;%.2f\n",global_best_thread,time_spent,t,global,global_energy);
                        fclose(floc);
                    #endif
                    
                }
                quiet=0;
            } // ---output
        } 
    }
};
//-------------------------------------------------------------------------------------------------------------------
TFloat median(vector<TFloat> &v)
{
    size_t nq = v.size() / 2;
    nth_element(v.begin(), v.begin()+nq, v.end());
    return v[n];
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//-------------------------------------------------------------------------------------------------------------------
// ODE observer (called after each integration step):
struct write_ode
{
    int m_id; // thread ID
    write_ode( int &id ) : m_id( id ) { }
    void operator()( state_type &x , double t )
    {
        // increase stepcounter for this thread:
        stepcounter[m_id]++;

        // apply unit clauses:
        if (unit_clauses) {
            for (int i=0; i<n; i++) {
                if (unit_clause_vars[i+1]!=0) x[i] = unit_clause_vars[i+1]; //TODO: assign +1 or -1
            }
        }

        //boundary for v € [-1,1]:
        if (forcebounds) {
            for (int i=0;i<n;i++) {
                if (x[i]<-1.0) x[i]=-1.0;   
                if (x[i]>1.0) x[i]=1.0;
            }
            //boundary for xs € [0,1]:
            for (int i=n;i<n+m;i++) {
                if (x[i]<0.0) x[i]=0.0;
                if (x[i]>1.0) x[i]=1.0;
            }
            //boundary for xl € [1,10⁴M]:
            for (int i=n+m;i<n+m*2;i++) {
                if (x[i]<1.0) x[i]=1.0;  
                if (x[i]>xl_max) x[i]=xl_max;
            } 
        }

    }
};

std::mt19937 randsolution_generator(seed);
std::uniform_int_distribution<int>  rand_solution(0, 100);

//-------------------------------------------------------------------------------------------------------------------
// function to hash an assignment etc to find out if we where previously there...
std::map<uint64_t, bool> db;

uint64_t hashing(TFloat *x) {
    uint64_t hash{0};
    for (int i=0; i<n; i++) {
        int assign = 0; if (x[i]>0.0) assign = 1;
        hash ^= assign + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
}

//-------------------------------------------------------------------------------------------------------------------
// the main solver sequence
// 
// input:               node
// uses global:         initial_assignments[]
// uses global:         thread_params[]
// output: 1 if solution is found, -1 if no solution is found

int dmm(struct node *node) {
    if (!quiet) printf("c [%d] STARTING ODE...\n",node->id);
    bool sat;
    
    pthread_mutex_lock(&mut);
        solved = 0;
        global = m;    
        global_energy = m;
        c_key_best = 0;
        stepcounter[node->id] = 0;
    pthread_mutex_unlock(&mut);
    global_thread[node->id] = m;
    time_thread[node->id] = TFloat(0.0);
    walltime_thread[node->id] = 0.0;
    energy_thread[node->id] = m;
    energy_thread_min[node->id] = m;


    for (int i=0; i<m; i++) {
        dmm_alpha_cm[node->id*m+i] = thread_params[node->id*PARAM_COUNT+0];
    }

    //defining vector:
    int size_of_vector = n+m*2;
    state_type x(size_of_vector);
    //if (!quiet) printf("c [%d] STATE VECTOR CREATED.\n",node->id);

    // initial conditions:
    for (int i=0; i<n+m*2; i++) x[i] = initial_assignments[node->id*(n+m*2)+i]; 
    
    if (!quiet) {
        std::cout << TEXT_DEFAULT << "c [" << node->id << "] INITIAL ASSIGNMENTS SET: "
        << std::setprecision(2) << std::fixed << TEXT_CYAN
        << x[0] << " "
        << x[1] << " "
        << x[2] << " "
        << x[3] << " "
        << x[4] << " "
        << x[5] << " "
        << x[6] << " "
        << x[7] << " "
        << x[8] << " "
        << x[9] << " "
        << TEXT_DEFAULT << std::endl;
        fflush(stdout);
    }

    TFloat t = t_init[node->id]; 
    t_begin_thread[node->id] = pt::microsec_clock::local_time().time_of_day().total_milliseconds(); //clock();

    // ---------------------------------------------------------------------------------------------
    if (INTEGRATION_MODE==ODE_RUNGEKUTTA) {
        if (!quiet) printf("c [%d] ADAPTIVE RUNGE-KUTTA-CASH-KARP54 INTEGRATION...\n",node->id);

        /* abs_tol, rel_tol, (max_dt) , custom_stepper */
        auto stepper = make_controlled( thread_params[node->id*PARAM_COUNT+6] , thread_params[node->id*PARAM_COUNT+7] , runge_kutta_cash_karp54< state_type , value_type >() );

        if (!quiet) {
            printf("c [%d] PARAMETERS: ",node->id);
            std::cout << std::setprecision(digits) << std::fixed;
            std::cout << " α=" << thread_params[node->id*PARAM_COUNT];
            std::cout << " β=" << thread_params[node->id*PARAM_COUNT+1];
            std::cout << " γ=" << thread_params[node->id*PARAM_COUNT+2];
            std::cout << " ε=" << thread_params[node->id*PARAM_COUNT+3];
            std::cout << " δ=" << thread_params[node->id*PARAM_COUNT+4];
            std::cout << " ζ=" << thread_params[node->id*PARAM_COUNT+5];
            std::cout << " E1=" << thread_params[node->id*PARAM_COUNT+6];
            std::cout << " E2=" << thread_params[node->id*PARAM_COUNT+7];
            std::cout << " t=" << t_init[node->id];
            std::cout << " init_dt=" << thread_params[node->id*PARAM_COUNT+8] << std::endl;
        }

        auto rhs = odeS(node->id, thread_params[node->id*PARAM_COUNT+0], thread_params[node->id*PARAM_COUNT+1], thread_params[node->id*PARAM_COUNT+2], thread_params[node->id*PARAM_COUNT+3], thread_params[node->id*PARAM_COUNT+4], thread_params[node->id*PARAM_COUNT+5]);
        
        integrate_adaptive(stepper, rhs , x, t, maxsteps, thread_params[node->id*PARAM_COUNT+8], write_ode(node->id));
    }

    // ---------------------------------------------------------------------------------------------
    if (INTEGRATION_MODE==ODE_ONE_STEP) {
        if (!quiet) printf("c [%d] INTEGRATING ONE ODE STEP...\n",node->id);
        euler< state_type, value_type> stepper_euler; //define integration method
        
        if (!quiet) {
            printf("c [%d] PARAMETERS: ",node->id);
            std::cout << std::setprecision(digits) << std::fixed;
            std::cout << " α=" << thread_params[node->id*PARAM_COUNT];
            std::cout << " β=" << thread_params[node->id*PARAM_COUNT+1];
            std::cout << " γ=" << thread_params[node->id*PARAM_COUNT+2];
            std::cout << " ε=" << thread_params[node->id*PARAM_COUNT+3];
            std::cout << " δ=" << thread_params[node->id*PARAM_COUNT+4];
            std::cout << " ζ=" << thread_params[node->id*PARAM_COUNT+5];
            std::cout << " E1=" << thread_params[node->id*PARAM_COUNT+6];
            std::cout << " E2=" << thread_params[node->id*PARAM_COUNT+7];
            std::cout << " t=" << t_init[node->id];
            std::cout << " init_dt=" << thread_params[node->id*PARAM_COUNT+8] << std::endl;
        }
        //integrate 1 time step until we are done
        integrate_n_steps(stepper_euler, odeS(node->id, thread_params[node->id*PARAM_COUNT+0], thread_params[node->id*PARAM_COUNT+1], thread_params[node->id*PARAM_COUNT+2], thread_params[node->id*PARAM_COUNT+3], thread_params[node->id*PARAM_COUNT+4], thread_params[node->id*PARAM_COUNT+5]), x, t, thread_params[node->id*PARAM_COUNT+8], 1, write_ode(node->id));
        
    }

    // ---------------------------------------------------------------------------------------------
    if (INTEGRATION_MODE==ODE_IMPLICIT) {
        if (!quiet) printf("c [%d] IMPLICIT ODE INTEGRATION...\n",node->id);
        
        /*
        eps_abs         Absolute tolerance level.
        eps_rel         Relative tolerance level.
        factor_dxdt     Factor for the weight of the derivative.
        factor_x        Factor for the weight of the state. */
        bulirsch_stoer< state_type , value_type > stepper( thread_params[node->id*PARAM_COUNT+6] , thread_params[node->id*PARAM_COUNT+7] , 1.0 , 1.0 ); //1e-4, 0.0, 0.0, 1.0
        
        if (!quiet) {
            printf("c [%d] PARAMETERS: ",node->id);
            std::cout << std::setprecision(digits) << std::fixed;
            std::cout << " α=" << thread_params[node->id*PARAM_COUNT];
            std::cout << " β=" << thread_params[node->id*PARAM_COUNT+1];
            std::cout << " γ=" << thread_params[node->id*PARAM_COUNT+2];
            std::cout << " ε=" << thread_params[node->id*PARAM_COUNT+3];
            std::cout << " δ=" << thread_params[node->id*PARAM_COUNT+4];
            std::cout << " ζ=" << thread_params[node->id*PARAM_COUNT+5];
            std::cout << " E1=" << thread_params[node->id*PARAM_COUNT+6];
            std::cout << " E2=" << thread_params[node->id*PARAM_COUNT+7];
            std::cout << " t=" << t_init[node->id];
            std::cout << " init_dt=" << thread_params[node->id*PARAM_COUNT+8] << std::endl;
        }

        auto rhs = odeS(node->id, thread_params[node->id*PARAM_COUNT+0], thread_params[node->id*PARAM_COUNT+1], thread_params[node->id*PARAM_COUNT+2], thread_params[node->id*PARAM_COUNT+3], thread_params[node->id*PARAM_COUNT+4], thread_params[node->id*PARAM_COUNT+5]);
        
        integrate_adaptive(stepper, rhs, x, t, maxsteps, thread_params[node->id*PARAM_COUNT+8], write_ode(node->id) );
    }

    // ---------------------------------------------------------------------------------------------
    if (INTEGRATION_MODE==ODE_CUSTOM_ADAPTIVE) {
        if (!quiet) printf("c [%d] CUSTOM ADAPTIVE STEP INTEGRATION...\n",node->id);
        if (!quiet) {
            printf("c [%d] PARAMETERS: ",node->id);
            std::cout << std::setprecision(digits) << std::fixed;
            std::cout << " α=" << thread_params[node->id*PARAM_COUNT];
            std::cout << " β=" << thread_params[node->id*PARAM_COUNT+1];
            std::cout << " γ=" << thread_params[node->id*PARAM_COUNT+2];
            std::cout << " ε=" << thread_params[node->id*PARAM_COUNT+3];
            std::cout << " δ=" << thread_params[node->id*PARAM_COUNT+4];
            std::cout << " ζ=" << thread_params[node->id*PARAM_COUNT+5];
            std::cout << " t=" << t_init[node->id];
            std::cout << " stepsize=" << thread_params[node->id*PARAM_COUNT+8] << std::endl;
        }
        euler< state_type, value_type> stepper_euler; //define integration method
        auto rhs = odeS(node->id, thread_params[node->id*PARAM_COUNT+0], thread_params[node->id*PARAM_COUNT+1], thread_params[node->id*PARAM_COUNT+2], thread_params[node->id*PARAM_COUNT+3], thread_params[node->id*PARAM_COUNT+4], thread_params[node->id*PARAM_COUNT+5]);

        forcebounds = false; // we are enforcing them in this code right here...
        TFloat min_dt = thread_params[node->id*PARAM_COUNT+8];

        while (solved!=1 && t<timeout) {
        
            stepper_euler.do_step(rhs, x, t, thread_params[node->id*PARAM_COUNT+8]); // updates x
            t = t + thread_params[node->id*PARAM_COUNT+8];
            stepcounter[node->id]++;
            
            // apply boundary for V:
            for (int i=0;i<n;i++) {
                if (x[i]<-1.0) x[i]=-1.0;   
                if (x[i]>1.0) x[i]=1.0;
            }
            //boundary for xs € [0,1]:
            for (int i=n;i<n+m;i++) {
                if (x[i]<0.0) x[i]=0.0;
                if (x[i]>1.0) x[i]=1.0;
            }
            //boundary for xl € [1,10⁴M]:
            for (int i=n+m;i<n+m*2;i++) {
                if (x[i]<1.0) x[i]=1.0;  
                if (x[i]>xl_max) x[i]=xl_max;
            } 

        }
    }
    
    // ---------------------------------------------------------------------------------------------
    if (INTEGRATION_MODE==ODE_CONSTANT) {
        if (!quiet) printf("c [%d] CONSTANT STEP INTEGRATION...\n",node->id);
        euler< state_type, value_type> stepper_euler; //define integration method
        
        if (!quiet) {
            printf("c [%d] PARAMETERS: ",node->id);
            std::cout << " α=" << thread_params[node->id*PARAM_COUNT];
            std::cout << " β=" << thread_params[node->id*PARAM_COUNT+1];
            std::cout << " γ=" << thread_params[node->id*PARAM_COUNT+2];
            std::cout << " ε=" << thread_params[node->id*PARAM_COUNT+3];
            std::cout << " δ=" << thread_params[node->id*PARAM_COUNT+4];
            std::cout << " ζ=" << thread_params[node->id*PARAM_COUNT+5];
            std::cout << " t=" << t_init[node->id];
            std::cout << " stepsize=" << thread_params[node->id*PARAM_COUNT+8] << std::endl;
        }

        integrate_n_steps(stepper_euler, odeS(node->id, thread_params[node->id*PARAM_COUNT+0], thread_params[node->id*PARAM_COUNT+1], thread_params[node->id*PARAM_COUNT+2], thread_params[node->id*PARAM_COUNT+3], thread_params[node->id*PARAM_COUNT+4], thread_params[node->id*PARAM_COUNT+5]), x, t, thread_params[node->id*PARAM_COUNT+8], TFloat(maxsteps), write_ode(node->id));
        
    }
    // --------------------------------------------------------------------------------------------- 
    if (INTEGRATION_MODE==ODE_CONSTANT_NOBOOST) {
        if (!quiet) printf("c [%d] CONSTANT STEP INTEGRATION (NO BOOST)...\n",node->id);
        
        TFloat t = 0.0; // start time
        TFloat h = thread_params[node->id*PARAM_COUNT+8]; //stepsize

        while (solved!=1) {
            
            if (constand_adaptive) {
                /// ADAPTIVE TIMESTEP:
                auto dxdt = dmm_generate_dxdt(node->id, x, t);
                
                bool adaptive_accepted = false;
                TFloat h_min = 0.0078125;
                TFloat h_max = 10000; 
                h = 0.125; //1.0;
                auto x_tmp = x;
                int varchanges2;
                
                while (!adaptive_accepted) {
                    x_tmp = dmm_dostep(node->id, x, dxdt, h);

                    varchanges2 = 0; // # vars changing by x %
                    for (int i=0; i<n; i++) {
                        if (fabs(x_tmp[i]-x[i])>=1) varchanges2++; // maximum allowed voltage change
                    }
                    if (varchanges2>0) {h = h * 1/2; } else {adaptive_accepted = true;}
                    stepcounter[node->id]--;
                    if (h<=h_min) adaptive_accepted = true;
                }
                x = x_tmp;
                stepcounter[node->id]++;
            } else {
                x = dmm_dostep(node->id, x, dmm_generate_dxdt(node->id, x, t), h); //generate dxdt and do forward euler step    
            }
            
            t = t + h; // increase ODE time by stepsize

            // solved? if yes, we are done
            if (solved) {
                //move solution to node:
                for (int i=0; i<n; i++) {
                    if (x[i]>0) node->optimal[i] = 1;
                    if (x[i]<0) node->optimal[i] = 0;
                }
                return 1;
            }
            
            // test heuristics:
            if (heuristics) apply_heuristics = stepcounter[node->id] % 10000 == 0;
            if (apply_heuristics) {
                        apply_heuristics = false;
                        //alpha heuristic
                        
                        state_type Xl_vec(m);
                        int m_id = node->id;
                        for (int i=0; i<m; i++) {
                            TFloat _Xl = x[n+m+i];
                            Xl_vec[i] = _Xl; // Xl => Xl_vec
                        }
                        TFloat Xl_median = median(Xl_vec);
                        // adust alpha/cm:
                        int changecounter = 0, maxcounter = 0, upcounter = 0, downcounter = 0;
                        for (int i=0; i<m; i++) {
                            TFloat _Xl = x[n+m+i];
                            if (_Xl>Xl_median) {dmm_alpha_cm[m_id*m+i] *= TFloat(1.1); changecounter++; upcounter++;}
                            if (_Xl<=Xl_median) {dmm_alpha_cm[m_id*m+i] *= TFloat(0.9); changecounter++; downcounter++;}
                            if (dmm_alpha_cm[m_id*m+i]<thread_params[node->id*PARAM_COUNT+0]) dmm_alpha_cm[m_id*m+i] = thread_params[node->id*PARAM_COUNT+0];
                            if (x[n+m+i]>=xl_max) {
                                dmm_alpha_cm[m_id*m+i] = thread_params[m_id*PARAM_COUNT+0];
                                x[n+m+i] = dmm_alpha;
                                maxcounter++;
                            }
                        }
                        printf("\n%sc [%d] alpha heuristic changed %d up, %d down, %d max reached (median = %.2f)%s\n",TEXT_CYAN,node->id,upcounter,downcounter,maxcounter,Xl_median,TEXT_DEFAULT);

            }
        }
    }

    if (solved) {
        //move solution to node:
        for (int i=0; i<n; i++) {
            if (x[i]>0) node->optimal[i] = 1;
            if (x[i]<0) node->optimal[i] = 0;
        }
        return 1;
    }

    //free x:
    x = std::vector<TFloat>();
    
    return 0;
}
//-------------------------------------------------------------------------------------------------------------------
// this functions initiates a node, runs the solving sequence and prints the result
// usually also used for multi threading
// input: thread-id
//int apply(int id) {
void *apply(void *threadid) {
    uint64_t tid;
    tid = (uint64_t)threadid;
    int id = int(tid);
    int i,result;

    struct node node;
    node.id = id;
    node.temporal = (int *) calloc((size_t) n, sizeof(int));
    node.model = (int *) calloc((size_t) n, sizeof(int));
    node.optimal = (int *) calloc((size_t) n, sizeof(int));

    //run ODE:
    result = dmm(&node); // <== run ODE integration
    t_end_thread[id] = pt::microsec_clock::local_time().time_of_day().total_milliseconds(); //clock();
    
    if (!quiet) printf(TEXT_DEFAULT);
    double time_spent = (double)(t_end_thread[id] - t_begin_thread[id]) / 1000; //CLOCKS_PER_SEC;
    if (!quiet) {
        printf("c [%d] TIME SPENT: %.5fs\n",id,time_spent);
        if (global_thread[id]!=global) {
            printf("c [%d] BEST LOC = %d (GLOBAL %d)\n",id,global_thread[id],global);
        } else {
            printf("c [%d] BEST LOC = %d (GLOBAL %d) **\n",id,global_thread[id],global);
        }
    }

    // print solution:
    if (result == 1) {
        if (!quiet) printf("\ns [%d] SATISFIABLE",id);
        for (i = 0; i < n; i++) {
            if (i % 20 == 0) if (!quiet) printf("\nv ");
            if (!quiet) printf("%i ", node.optimal[i] ? +(i + 1) : -(i + 1));
        }
        if (!quiet) printf("0\n");
        fflush(stdout);

        
    }
    if (result == 0) {
        if (!quiet) printf("\ns [%d] UNKNOWN\n",id);
    }

    // free memory:
    free(node.temporal);
    free(node.model);
    free(node.optimal);

    pthread_mutex_lock(&mut);
        running_threads--;
    pthread_mutex_unlock(&mut);
    
    pthread_exit(NULL);
}

//-------------------------------------------------------------------------------------------------------------------
// runs all threads with the set parameters and wits until finished
// uses:
// parameters: thread_params[]
// initial values: initial_assignments[]
// returns:
// best loc of each thread => global_thread[]
void run_all_threads() {
    // run solving process multi-threaded:
    pthread_t threads[THREAD_COUNT];
    int rem_quiet = quiet;
    quiet = 1;
    for(int i = 0; i < THREAD_COUNT; i++ ) {
      int rc = pthread_create(&threads[i], NULL, apply, (void *)i);
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
      running_threads++;
    }
    // display loop while threads are running:
    while (running_threads > 0) {
        if (!rem_quiet) {
            printf("\rc [TUNE] %d THREADS: ",running_threads);
            for (int i=0; i<THREAD_COUNT; i++) {
                std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT;
                std::cout << time_thread_actual[i] << " [" << energy_thread_min[i] << "] ";
            }
            fflush(stdout);
        }
    }
    if (!rem_quiet) printf("\n");
    quiet = rem_quiet;
}

//-------------------------------------------------------------------------------------------------------------------
// tunes and find optimal settings
/* tune_mode:
   0 : parameter search (always from all-zero-V assignment to find ideal parameters)
   1 : solution search  (keeps best V-assignment at restart, pertubed with switchfraction) */
void tune_circuit() {

    /* Simulated Annealing */
    /* - we use 8 monte carlo chains, running for x iterations
       - each iterations start with an initial best_eval = loc
       - each monte carlo chain runs in one seperate thread
       - switchfraction = 0.0001 (percentage of assignments which are being pertubed
       - we will use hard and soft upper- and lower bounds
       - we change one or two params in each iteration
       - accepted params: 0;1 (after each iteration)
       - accepted assignments: 0;1 (after each iteration)
       )
    */

    /* tuning parameters */
    double  INITIAL_TEMP    = 1.00; //10

    if (!quiet) {
        printf("c [TUNE] MODE           : %d\n",tune_mode);
        printf("c [TUNE] N_ITERATIONS   : %d\n",N_ITERATIONS);
        printf("c [TUNE] INITIAL_TEMP   : %.2f\n",INITIAL_TEMP);
        printf("c [TUNE] SWITCHFRACTION : %.15f\n",switchfraction);
        if (tune_mode_params==0) printf("c [TUNE] PARAMS: TUNE PARAMS ALPHA-ZETA\n");
        if (tune_mode_params==1) printf("c [TUNE] PARAMS: TUNE PARAMS OF ODE\n");
        if (tune_mode_params==2) printf("c [TUNE] PARAMS: TUNE ALL PARAMS\n");
    }
     
    /* parameter ranges */
    double param_bounds[THREAD_COUNT][PARAM_COUNT][4];

    for (int i=0; i<THREAD_COUNT; i++) {
        if (1==1) {//!partable_loaded) { //ignore partable bounds for now
            //alpha
            param_bounds[i][ALPHA][LBH] = (0.000000000000001);    param_bounds[i][ALPHA][HBH] = (100.0);
            param_bounds[i][ALPHA][LBS] = (0.01);                 param_bounds[i][ALPHA][HBS] = (1.0);
            //beta
            param_bounds[i][BETA][LBH] = (0.000000000000001);     param_bounds[i][BETA][HBH] = (100.0);
            param_bounds[i][BETA][LBS] = (0.01);                  param_bounds[i][BETA][HBS] = (1.0);
            //gamma
            param_bounds[i][GAMMA][LBH] = (0.000000000000001);    param_bounds[i][GAMMA][HBH] = (0.5);
            param_bounds[i][GAMMA][LBS] = (0.01);                 param_bounds[i][GAMMA][HBS] = (0.25);
            //delta
            param_bounds[i][DELTA][LBH] = (0.000000000000001);    param_bounds[i][DELTA][HBH] = (1.0);
            param_bounds[i][DELTA][LBS] = (0.01);                 param_bounds[i][DELTA][HBS] = (0.25);
            //epsilon
            param_bounds[i][EPSILON][LBH] = (0.000000000000001);  param_bounds[i][EPSILON][HBH] = (0.5);
            param_bounds[i][EPSILON][LBS] = (0.01);               param_bounds[i][EPSILON][HBS] = (0.25);
            //zeta
            param_bounds[i][ZETA][LBH] = (0.000000000000001);     param_bounds[i][ZETA][HBH] = (1.0);
            param_bounds[i][ZETA][LBS] = (0.01);                  param_bounds[i][ZETA][HBS] = (0.25);
            //err1
            param_bounds[i][ERR1][LBH] = (0.000000000000001);     param_bounds[i][ERR1][HBH] = (1.00);
            param_bounds[i][ERR1][LBS] = (0.01);                  param_bounds[i][ERR1][HBS] = (1.00);
            //err2
            param_bounds[i][ERR2][LBH] = (0.000000000000001);     param_bounds[i][ERR2][HBH] = (1.00);
            param_bounds[i][ERR2][LBS] = (0.01);                  param_bounds[i][ERR2][HBS] = (1.00);
            //init_dt
            param_bounds[i][INITDT][LBH] = (0.000000000000001);   param_bounds[i][INITDT][HBH] = (1.00);
            param_bounds[i][INITDT][LBS] = (0.000000000000001);   param_bounds[i][INITDT][HBS] = (0.00001);
            if (!quiet) printf("c [TUNE] [%d] USING DEFAULT (LARGE) BOUNDS.\n",i);
        } else {
            //alpha
            param_bounds[i][ALPHA][LBH] = partable[0*4+2];   param_bounds[i][ALPHA][HBH] = partable[0*4+3];
            param_bounds[i][ALPHA][LBS] = partable[0*4+0];   param_bounds[i][ALPHA][HBS] = partable[0*4+1];
            //beta
            param_bounds[i][BETA][LBH] = partable[1*4+2];    param_bounds[i][BETA][HBH] = partable[1*4+3];
            param_bounds[i][BETA][LBS] = partable[1*4+0];    param_bounds[i][BETA][HBS] = partable[1*4+1];
            //gamma
            param_bounds[i][GAMMA][LBH] = partable[2*4+2];   param_bounds[i][GAMMA][HBH] = partable[2*4+3];
            param_bounds[i][GAMMA][LBS] = partable[2*4+0];   param_bounds[i][GAMMA][HBS] = partable[2*4+1];
            //delta
            param_bounds[i][DELTA][LBH] = partable[3*4+2];   param_bounds[i][DELTA][HBH] = partable[3*4+3];
            param_bounds[i][DELTA][LBS] = partable[3*4+0];   param_bounds[i][DELTA][HBS] = partable[3*4+1];
            //epsilon
            param_bounds[i][EPSILON][LBH] = partable[4*4+2]; param_bounds[i][EPSILON][HBH] = partable[4*4+3];
            param_bounds[i][EPSILON][LBS] = partable[4*4+0]; param_bounds[i][EPSILON][HBS] = partable[4*4+1];
            //zeta
            param_bounds[i][ZETA][LBH] = partable[5*4+2];    param_bounds[i][ZETA][HBH] = partable[5*4+3];
            param_bounds[i][ZETA][LBS] = partable[5*4+0];    param_bounds[i][ZETA][HBS] = partable[5*4+1];
            //err1
            param_bounds[i][ERR1][LBH] = partable[6*4+2];    param_bounds[i][ERR1][HBH] = partable[6*4+3];
            param_bounds[i][ERR1][LBS] = partable[6*4+0];    param_bounds[i][ERR1][HBS] = partable[6*4+1];
            //err2
            param_bounds[i][ERR2][LBH] = partable[7*4+2];    param_bounds[i][ERR2][HBH] = partable[7*4+3];
            param_bounds[i][ERR2][LBS] = partable[7*4+0];    param_bounds[i][ERR2][HBS] = partable[7*4+1];
            //init_dt
            param_bounds[i][INITDT][LBH] = partable[8*4+2];  param_bounds[i][INITDT][HBH] = partable[8*4+3];
            param_bounds[i][INITDT][LBS] = partable[8*4+0];  param_bounds[i][INITDT][HBS] = partable[8*4+1];
            if (!quiet) printf("c [TUNE] [%d] USING partable.txt BOUNDS.\n",i);
        }
    }

    /* variables */
    TFloat    best_eval[THREAD_COUNT];
    TFloat    curr_eval[THREAD_COUNT];  
    TFloat    candidate_eval[THREAD_COUNT];  
    TFloat    rem_params[THREAD_COUNT*PARAM_COUNT];
    
    int       num_change_vars[THREAD_COUNT];
    int       tune_param[THREAD_COUNT][PARAM_COUNT];
    int       rem_tune_param[THREAD_COUNT][PARAM_COUNT]; 
    for (int i=0; i<THREAD_COUNT;i++) for (int p=0; p<PARAM_COUNT;p++) rem_tune_param[i][p] = -1;
    
    TFloat  * rem_assignments; rem_assignments = (TFloat *) calloc((size_t) THREAD_COUNT*(n+m*2), sizeof(TFloat) );
    int       number_no_improvements[THREAD_COUNT]; 
    for  (int i=0; i<THREAD_COUNT;i++) number_no_improvements[i] = 0;
    
    TFloat    tune_global = m; 
    double    SMALL_STEP_FRACTION[THREAD_COUNT];
    for (int i=0; i<THREAD_COUNT; i++) SMALL_STEP_FRACTION[i]= 100;
    int       switch_count_vars;
    
    /* RNG */
    std::mt19937 generator(seed);
    
    /// depending on mode: 0=only params, 1=ODE integration params or 2=all:
    int p_low = 0;
    int p_high = PARAM_COUNT-1;
    if (tune_mode_params==0) {p_low=0; p_high=5;} 
    if (tune_mode_params==1) {p_low=6; p_high=PARAM_COUNT-1;} 
    if (tune_mode_params==2) {p_low=0; p_high=PARAM_COUNT-1;}
    if (tune_mode_params==1 && INTEGRATION_MODE==ODE_CONSTANT) {p_low=PARAM_COUNT-1; p_high=PARAM_COUNT-1;} 
    std::uniform_int_distribution<int>      rand_param(p_low,p_high); //which params to change 
    std::uniform_real_distribution<double>  rand_r(0.0, 1.0);
    std::uniform_int_distribution<int>      rand_switchfraction(0, n-1);
    std::uniform_int_distribution<int>      rand_v(-1, 1);
    std::uniform_real_distribution<double>  rand_Xs(0.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xl(1.0, 10.0); //xl_max); //10.00
    std::uniform_int_distribution<int>      rand_mode(0, 2);
    
    /* move assigments => rem_assignments[] */
    for (int i=0; i<THREAD_COUNT; i++)
        for (int j=0; j<n+m*2; j++)
            rem_assignments[i*(n+m*2)+j] = initial_assignments[i*(n+m*2)+j];
    if (!quiet) printf("c [TUNE] STARTING WITH initial_assignments ASSIGNMENTS.\n"); 
    
    /* find initial best_eval for each thread by performing one step:*/
    int rem_INT_MODE = INTEGRATION_MODE;
    INTEGRATION_MODE=ODE_ONE_STEP;
    if (!quiet) printf("\rc [TUNE] INITIAL OBJECTIVE...\n");
    run_all_threads();
    //output best_eval/loc for each thread:
    if (!quiet) printf("c [TUNE] E: "); 
    for (int i=0; i<THREAD_COUNT; i++) {
        if (!quiet) std::cout << "[" << energy_thread_min[i] << "]"; // printf("[%.2f] ",energy_thread_min[i]); 
        best_eval[i] = energy_thread_min[i];
        curr_eval[i] = best_eval[i];
    }
    if (!quiet) printf("               \n");
    
    /* MAIN TUNING ITERAZION LOOP ------------------------------------------------------------------------------------------------ */

    bool change_params = false;

    for (int iter=0; iter<N_ITERATIONS; iter++) {

        if (quiet) std::cout << "c [TUNE] ITERATION " << iter << ". GLOBAL=" << tune_global << std::endl; // printf(" c [TUNE] INTERATION %d. GLOBAL=%.2f\n",iter,tune_global);

        int min_number_no_improvements = N_ITERATIONS; 

        INTEGRATION_MODE=rem_INT_MODE;
        
        /* temp store data and apply switchfraction: */
        for (int i=0; i<THREAD_COUNT; i++) {
            if (!quiet) printf("c [TUNE] [%d] TUNING ITERATION %d - NO CANDIDATE ACCEPTED SINCE: %d IT.\n",i,iter,number_no_improvements[i]);
            // calculate min of all number_no_improvements:
            if (number_no_improvements[i]<min_number_no_improvements) min_number_no_improvements = number_no_improvements[i];

            // current parameters => rem_params[]:
            for (int j=0; j<PARAM_COUNT; j++) rem_params[i*PARAM_COUNT+j] = thread_params[i*PARAM_COUNT+j];
            
            // mode 0: always start from original assignment ---------------------------------------------------------------------
            if (tune_mode==0) {
                for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
                if (!quiet) {
                    printf("c [TUNE] [%d] USING STARTING ASSIGNMENT (V, Xs, Xl) ",i);
                    std::cout << initial_assignments[i*(n+m*2)+0] << " "
                    << initial_assignments[i*(n+m*2)+1] << " "
                    << initial_assignments[i*(n+m*2)+2] << " "
                    << initial_assignments[i*(n+m*2)+3] << " "
                    << initial_assignments[i*(n+m*2)+4] << std::endl;
                }
                t_init[i] = 0.0;
            }
            
            // mode 1: re-start from best Xl ----------------------------------------------------------------------------------
            if (tune_mode==1 && iter>0) {
                
                for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
                for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];//0.0;// reset Xs
                for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = v_best[i*(n+m*2)+j];//1.0;// reset Xl
                if (!quiet) {
                    printf("c [TUNE] [%d] USING BEST Xl ASSIGNMENT ",i);
                    std::cout << initial_assignments[i*(n+m*2)+0] << " "
                    << initial_assignments[i*(n+m*2)+1] << " "
                    << initial_assignments[i*(n+m*2)+2] << " "
                    << initial_assignments[i*(n+m*2)+3] << " "
                    << initial_assignments[i*(n+m*2)+4] << std::endl;
                }
            }
            if (tune_mode==1 && iter==0) {
                for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
            }

            // mode 2: re-start from best params over all threads ------------------------------------------------
            if (tune_mode==2 && iter>0) {
                
                if (global_best_thread!=-1) {
                    // best global params (not just thread):
                    for (int j=0; j<PARAM_COUNT; j++) thread_params[i*PARAM_COUNT+j] = thread_params[global_best_thread*PARAM_COUNT+j];
                    // best globlal assignment (not just thread):
                    //for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = v_best[global_best_thread*(n+m*2)+j];
                    for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = v_best[global_best_thread*(n+m*2)+j];
                } 
                if (!quiet) {
                    printf("c [TUNE] [%d] USING OVERALL BEST PARAMETERS, ORIGINAL ASSIGNMENT ",i);
                    std::cout << initial_assignments[i*(n+m*2)+0] << " "
                    << initial_assignments[i*(n+m*2)+1] << " "
                    << initial_assignments[i*(n+m*2)+2] << " "
                    << initial_assignments[i*(n+m*2)+3] << " "
                    << initial_assignments[i*(n+m*2)+4] << std::endl;
                }
                t_init[i] = t_init[global_best_thread];
            }
            if (tune_mode==2 && iter==0) {
                for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
                t_init[i] = 0.0;
            }

            // mode 3: re-start with inital assignment and best params over all threads -----------------------------------------
            if (tune_mode==3 && iter>0) {
                
                if (global_best_thread!=-1) {
                    // best bloal params (not just thread):
                    for (int j=0; j<PARAM_COUNT; j++) thread_params[i*PARAM_COUNT+j] = thread_params[global_best_thread*PARAM_COUNT+j];
                    for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
                } 
                for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] = 0.0;               // reset Xs
                for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = 1.0;           // reset Xl
                if (!quiet) {
                    printf("c [TUNE] [%d] USING OVERALL BEST ASSIGNMENT ",i);
                    std::cout << initial_assignments[i*(n+m*2)+0] << " "
                    << initial_assignments[i*(n+m*2)+1] << " "
                    << initial_assignments[i*(n+m*2)+2] << " "
                    << initial_assignments[i*(n+m*2)+3] << " "
                    << initial_assignments[i*(n+m*2)+4] << std::endl;
                }
                t_init[i] = 0.0;
            }
            if (tune_mode==3 && iter==0) {
                for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
                t_init[i] = 0.0;
            }

            // mode 4: restart with lowest loc if it got better, otherwise just continue ---------------------------------------
            if (tune_mode==4 && iter>0) {
                if (global_thread[i]>global_all_runs_thread[i]) {
                    for (int j=0; j<n+m*2; j++) {
                        initial_assignments[i*(n+m*2)+j] = v_best[i*(n+m*2)+j];
                        rem_assignments[i*(n+m*2)+j] = v_best[i*(n+m*2)+j];
                    }
                    if (!quiet) {
                        printf("c [TUNE] [%d] CURRENT LOC IS WORSE, RESET TO v_best. ",i);
                        std::cout << initial_assignments[i*(n+m*2)+0] << " "
                        << initial_assignments[i*(n+m*2)+1] << " "
                        << initial_assignments[i*(n+m*2)+2] << " "
                        << initial_assignments[i*(n+m*2)+3] << " "
                        << initial_assignments[i*(n+m*2)+4] << std::endl;
                    }
                } else {
                    if (!quiet) {
                        printf("c [TUNE] [%d] CURRENT LOC IS EQUAL OR BETTER, CONTINUE WITH NO CHANGES ",i);
                        std::cout << initial_assignments[i*(n+m*2)+0] << " "
                        << initial_assignments[i*(n+m*2)+1] << " "
                        << initial_assignments[i*(n+m*2)+2] << " "
                        << initial_assignments[i*(n+m*2)+3] << " "
                        << initial_assignments[i*(n+m*2)+4] << std::endl;
                    }
                }
            }
            if (tune_mode==4 && iter==0) {
                for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rem_assignments[i*(n+m*2)+j];
            }

            // apply switchfraction: -------------------------------------------------------------------------------------------
            if (iter>0) {
                switch_count_vars = abs(n*switchfraction);//abs((n+m*2)*switchfraction);
                if (switchfraction>0 && switch_count_vars==0) switch_count_vars = 1;
                for (int j=0; j<switch_count_vars; j++) {
                    int change_var = rand_switchfraction(generator);
                    initial_assignments[i*(n+m*2)+change_var] = rand_v(generator);
                    //printf("%d->%.2f ",change_var,initial_assignments[i*(n+m*2)+change_var]);
                }
            }
            
        }
        if (iter>0 && !quiet) printf("c [TUNE] SWITCHFRACTION CHANGED %d ASSIGNMENTS.\n",switch_count_vars);

        // new: increase walltime_timeout if no one improved since XXX steps:
        if (min_number_no_improvements>=25) {
            //walltime_timeout += 1;
            // NEW:
            walltime_timeout *= 2;

            maxsteps = maxsteps * 2;
            if (!quiet) {
                printf(TEXT_RED);
                printf("c [TUNE] NOTHING IMPROVED SINCE 25 ITERATIONS. INCREASED walltime_timeout TO %.2f maxsteps TO ",walltime_timeout); std::cout << maxsteps << std::endl;

                // NEW:
                if (tune_mode==4) {
                    // we set all threads to params and assignments from the current best one:
                    // find best thread:
                    int global_best_thread_loc = m;
                    for (int i=0; i<THREAD_COUNT; i++) {
                        if (curr_eval[i]<global_best_thread_loc) {
                            global_best_thread_loc=(int)curr_eval[i]; 
                            global_best_thread=i;
                        }
                    }
                    for (int i=0; i<THREAD_COUNT; i++) {
                        if (i!=global_best_thread) {
                            // move best thread's params => thread:
                            for (int j=0; j<PARAM_COUNT; j++) {
                                thread_params[i*PARAM_COUNT+j] = rem_params[global_best_thread*PARAM_COUNT+j];
                                rem_params[i*PARAM_COUNT+j] = rem_params[global_best_thread*PARAM_COUNT+j];
                            }
                            // move best thread's assignment => thread:
                            for (int j=0; j<n+m*2;j++) {
                                initial_assignments[i*(n+m*2)+j] = rem_assignments[global_best_thread*(n+m*2)+j];
                                rem_assignments[i*(n+m*2)+j] = rem_assignments[global_best_thread*(n+m*2)+j];
                            }
                            // copy best thread's init_t => thread:
                            t_init[i] = t_init[global_best_thread];
                        }
                    }
                    printf("c [TUNE] ALSO ASSIGNED BEST PARAMS AND ASSIGNMENTS FROM [%d] TO ALL THREADS.\n",global_best_thread); 
                    }
                // --NEW

                printf(TEXT_DEFAULT);
            }
            for (int i=0; i<THREAD_COUNT; i++) number_no_improvements[i] = 0;
            for (int i=0; i<THREAD_COUNT; i++) SMALL_STEP_FRACTION[i] = 100.0;
            change_params = false; // re-run each thread with best params and more time
        }

        if (change_params) {
            if (!quiet) printf("c [TUNE] walltime_timeout=%.2f change_params=YES\n",walltime_timeout);
        } else {
            if (!quiet) printf("c [TUNE] walltime_timeout=%.2f change_params=NO\n",walltime_timeout);
        }
        
        // loop through each thread and pick a new param:
        int mode[THREAD_COUNT];
        for (int i=0; i<THREAD_COUNT; i++) {
            num_change_vars[i] = 0;
            for (int p=0; p<PARAM_COUNT; p++) rem_tune_param[i][p] = -1;
        }

        if (change_params) {
            
            for (int i=0; i<THREAD_COUNT; i++) {

                /* adjust SMALL_STEP_FRACTION 1e1, 1e2, 1e3 - 1e14*/
                if (number_no_improvements[i]>0 && number_no_improvements[i] % 10 == 0) {
                    SMALL_STEP_FRACTION[i] = SMALL_STEP_FRACTION[i] * 10;
                    if (SMALL_STEP_FRACTION[i]>=1000000000000) SMALL_STEP_FRACTION[i] = 1000000000000; //100;
                    if (!quiet) {
                        printf(TEXT_GREEN);
                        printf("C [TUNE] [%d] SMALL_STEP_FRACTION ADJUSTED TO %.0f AFTER %d UNSUCCESSFUL STEPS\n",i, SMALL_STEP_FRACTION[i],number_no_improvements[i]);
                        printf(TEXT_DEFAULT);
                    }
                }
                
                // choose mode:
                TFloat candidate;
                TFloat frac;
                mode[i] = rand_mode(generator);
                //if (tune_mode==4) mode[i] = 0; // in tune-mode 4 we only do small step moves

                // choose how many vars to change (log-normal distributed, mostly 1)
                std::lognormal_distribution<double> rand_num_vars(0.0,0.5);
                TFloat _num_change_vars = rand_num_vars(generator);
                num_change_vars[i] = (int)_num_change_vars; 
                if (num_change_vars[i]<1) num_change_vars[i] = 1;
                if (num_change_vars[i]>PARAM_COUNT) num_change_vars[i] = PARAM_COUNT;
                // loop through number of params to change:
                for (int p=0; p<num_change_vars[i]; p++) {
                    //pick a param (if previosuly a param change was successful, try it again)
                    if (rem_tune_param[i][p]==-1) {
                        tune_param[i][p] = rand_param(generator);
                    } else {
                        tune_param[i][p] = rem_tune_param[i][p];
                    }
                    rem_tune_param[i][p] = tune_param[i][p]; //will be set to -1 if candidate rejected
                    
                    // small random move: uses lognormal distribution and changes the last 4 digits (mostly) -------------------------
                    if (mode[i]==0) {
                        //stepsize of param:
                        double soft_bound_range = param_bounds[i][tune_param[i][p]][HBS] - param_bounds[i][tune_param[i][p]][LBS];
                        double param_stepsize = soft_bound_range / 100;
                        //adjust param:
                        std::uniform_real_distribution<double> rand_p(0,soft_bound_range);
                        double _rand = rand_p(generator);
                        _rand = _rand * 100 / SMALL_STEP_FRACTION[i];
                        candidate = thread_params[i*PARAM_COUNT+tune_param[i][p]] + _rand * param_stepsize;

                        std::uniform_real_distribution<double> rand_move(-1.5,1.5);
                        double move = rand_move(generator);
                        candidate = thread_params[i*PARAM_COUNT+tune_param[i][p]] + thread_params[i*PARAM_COUNT+tune_param[i][p]] * move/SMALL_STEP_FRACTION[i];

                        // bounds violated:
                        if (candidate<=0.0) candidate = 1e-8;

                        if (!quiet) {
                            printf("%sc [TUNE] [%d] SMALL MOVE: p_%d ",TEXT_CYAN,i,tune_param[i][p]);
                            std::cout << std::setprecision(18) << std::fixed
                             << rem_params[i*PARAM_COUNT+tune_param[i][p]] << " > " << candidate;
                            std::cout << std::setprecision(5) << std::fixed
                             << TEXT_DEFAULT << " SMALL_STEP_FRACTION=" << SMALL_STEP_FRACTION[i] << std::endl;
                        }
                    }
                    
                    // random point in soft bounds: ---------------------------------------------------------------------------------
                    if (mode[i]==1) {
                        //random generators:
                        std::uniform_real_distribution<double>  rand_alpha(param_bounds[i][ALPHA][LBS], param_bounds[i][ALPHA][HBS]);
                        std::uniform_real_distribution<double>  rand_beta(param_bounds[i][BETA][LBS], param_bounds[i][BETA][HBS]);
                        std::uniform_real_distribution<double>  rand_gamma(param_bounds[i][GAMMA][LBS], param_bounds[i][GAMMA][HBS]);
                        std::uniform_real_distribution<double>  rand_delta(param_bounds[i][DELTA][LBS], param_bounds[i][DELTA][HBS]);
                        std::uniform_real_distribution<double>  rand_epsilon(param_bounds[i][EPSILON][LBS], param_bounds[i][EPSILON][HBS]);
                        std::uniform_real_distribution<double>  rand_zeta(param_bounds[i][ZETA][LBS], param_bounds[i][ZETA][HBS]);
                        std::uniform_real_distribution<double>  rand_err1(param_bounds[i][ERR1][LBS], param_bounds[i][ERR1][HBS]);
                        std::uniform_real_distribution<double>  rand_err2(param_bounds[i][ERR2][LBS], param_bounds[i][ERR2][HBS]);
                        std::uniform_real_distribution<double>  rand_initdt(param_bounds[i][INITDT][LBS], param_bounds[i][INITDT][HBS]);
                        if (tune_param[i][p]==0) candidate = rand_alpha(generator);
                        if (tune_param[i][p]==1) candidate = rand_beta(generator);
                        if (tune_param[i][p]==2) candidate = rand_gamma(generator);
                        if (tune_param[i][p]==3) candidate = rand_delta(generator);
                        if (tune_param[i][p]==4) candidate = rand_epsilon(generator);
                        if (tune_param[i][p]==5) candidate = rand_zeta(generator);
                        if (tune_param[i][p]==6) candidate = rand_err1(generator);
                        if (tune_param[i][p]==7) candidate = rand_err2(generator);
                        if (tune_param[i][p]==8) candidate = rand_initdt(generator);
                        
                        if (tune_param[i][p]==4 && candidate >= thread_params[i*PARAM_COUNT+2]) candidate = thread_params[i*PARAM_COUNT+2]-0.1;
                        if (!quiet) {
                            printf("%sc [TUNE] [%d] RANDOM POINT IN SOFT BOUNDS: p_%d ",TEXT_CYAN,i,tune_param[i][p]);
                            std::cout << rem_params[i*PARAM_COUNT+tune_param[i][p]] << " > " << candidate;
                            std::cout << TEXT_DEFAULT << std::endl;
                        }
                    }
                    
                    // random point in hard bounds: --------------------------------------------------------------------------------
                    if (mode[i]==2) {
                        //random generators:
                        std::uniform_real_distribution<double>  rand_alpha_h(param_bounds[i][ALPHA][LBH], param_bounds[i][ALPHA][HBH]);
                        std::uniform_real_distribution<double>  rand_beta_h(param_bounds[i][BETA][LBH], param_bounds[i][BETA][HBH]);
                        std::uniform_real_distribution<double>  rand_gamma_h(param_bounds[i][GAMMA][LBH], param_bounds[i][GAMMA][HBH]);
                        std::uniform_real_distribution<double>  rand_delta_h(param_bounds[i][DELTA][LBH], param_bounds[i][DELTA][HBH]);
                        std::uniform_real_distribution<double>  rand_epsilon_h(param_bounds[i][EPSILON][LBH], param_bounds[i][EPSILON][HBH]);
                        std::uniform_real_distribution<double>  rand_zeta_h(param_bounds[i][ZETA][LBH], param_bounds[i][ZETA][HBH]);
                        std::uniform_real_distribution<double>  rand_err1_h(param_bounds[i][ERR1][LBH], param_bounds[i][ERR1][HBH]);
                        std::uniform_real_distribution<double>  rand_err2_h(param_bounds[i][ERR2][LBH], param_bounds[i][ERR2][HBH]);
                        std::uniform_real_distribution<double>  rand_initdt_h(param_bounds[i][INITDT][LBH], param_bounds[i][INITDT][HBH]);
                        if (tune_param[i][p]==0) candidate = rand_alpha_h(generator);
                        if (tune_param[i][p]==1) candidate = rand_beta_h(generator);
                        if (tune_param[i][p]==2) candidate = rand_gamma_h(generator);
                        if (tune_param[i][p]==3) candidate = rand_delta_h(generator);
                        if (tune_param[i][p]==4) candidate = rand_epsilon_h(generator);
                        if (tune_param[i][p]==5) candidate = rand_zeta_h(generator);
                        if (tune_param[i][p]==6) candidate = rand_err1_h(generator);
                        if (tune_param[i][p]==7) candidate = rand_err2_h(generator);
                        if (tune_param[i][p]==8) candidate = rand_initdt_h(generator);
                        
                        if (tune_param[i][p]==4 && candidate >= thread_params[i*PARAM_COUNT+2]) candidate = thread_params[i*PARAM_COUNT+2]-0.1;
                        if (!quiet) {
                            printf("%sc [TUNE] [%d] RANDOM POINT IN HARD BOUNDS: p_%d ",TEXT_CYAN,i,tune_param[i][p]);
                            std::cout << rem_params[i*PARAM_COUNT+tune_param[i][p]] << " > " << candidate;
                            std::cout << TEXT_DEFAULT << std::endl;
                        }
                    }

                    // set param (only after first round):
                    if (mode[i]!=-1) thread_params[i*PARAM_COUNT+tune_param[i][p]] = candidate;
                }    
                
                // output params: 
                if (!quiet) {
                    std::cout << "c [TUNE] [" << i << "]";
                    std::cout << std::setprecision(digits) << std::fixed;
                    std::cout << " -c=" << thread_params[i*PARAM_COUNT+0];
                    std::cout << " -b=" << thread_params[i*PARAM_COUNT+1];
                    std::cout << " -n=" << thread_params[i*PARAM_COUNT+2];
                    std::cout << " -h=" << thread_params[i*PARAM_COUNT+3];
                    std::cout << " -j=" << thread_params[i*PARAM_COUNT+4];
                    std::cout << " -k=" << thread_params[i*PARAM_COUNT+5];
                    std::cout << " -x=" << thread_params[i*PARAM_COUNT+6];
                    std::cout << " -y=" << thread_params[i*PARAM_COUNT+7];
                    std::cout << " -d=" << thread_params[i*PARAM_COUNT+8] << std::endl;
                }
                
            }
        } else {
            if (!quiet) printf("%sc [TUNE] ALL PARAMETERS UNCHANGED.%s\n",TEXT_CYAN,TEXT_DEFAULT);
            for (int i=0; i<THREAD_COUNT; i++) {
                // output params:
                if (!quiet) {
                    std::cout << "c [TUNE] [" << i << "]";
                    std::cout << std::setprecision(digits) << std::fixed;
                    std::cout << " -c=" << thread_params[i*PARAM_COUNT+0];
                    std::cout << " -b=" << thread_params[i*PARAM_COUNT+1];
                    std::cout << " -n=" << thread_params[i*PARAM_COUNT+2];
                    std::cout << " -h=" << thread_params[i*PARAM_COUNT+3];
                    std::cout << " -j=" << thread_params[i*PARAM_COUNT+4];
                    std::cout << " -k=" << thread_params[i*PARAM_COUNT+5];
                    std::cout << " -x=" << thread_params[i*PARAM_COUNT+6];
                    std::cout << " -y=" << thread_params[i*PARAM_COUNT+7];
                    std::cout << " -d=" << thread_params[i*PARAM_COUNT+8] << std::endl;
                }
            }
        }

        //if (iter % 2 == 0) {change_params = true;} else {change_params = false;}
        change_params = true;

        // run simulation: ---------------------------------------------------------------------------------------------------
        t_begin = clock();
        run_all_threads();
        t_end = clock();
        double time_spent = (double)(t_end - t_begin) / CLOCKS_PER_SEC;

        //output global energy for each thread:
        if (!quiet) {
            printf("c [TUNE] E (GLOBAL): "); 
            for (int i=0; i<THREAD_COUNT; i++) {
                std::cout << "[" << curr_eval[i] << "] "; // printf("[%.2f] ",curr_eval[i]); 
            }
            printf("               \n");
        }

        //output best_eval/loc for each thread: -------------------------------------------------------------------------------
        int tune_global_best = -1;
        if (!quiet) printf("c [TUNE] E         : "); 
        for (int i=0; i<THREAD_COUNT; i++) {
            candidate_eval[i] = energy_thread_min[i];
            if (candidate_eval[i]<tune_global) {
                tune_global = candidate_eval[i];
                // we found a better overvall global - output parameters:
                if (!quiet) printf("%s NEW GLOBAL %s",TEXT_YELLOW,TEXT_DEFAULT);
                if (!quiet) std::cout << "[" << TEXT_SILVER << energy_thread_min[i] << TEXT_DEFAULT << "] "; //  printf("[%s%.2f%s] ",TEXT_SILVER,energy_thread_min[i],TEXT_DEFAULT); 
                #ifdef USEFPRINTF
                    FILE *ftune = fopen(TUNING_FILE, "a");
                    fprintf(ftune,"%d;%d;%.15f;%.15f;%.15f;%.15f;%.15f;%.15f;%.15f;%.15f\n",iter,i,thread_params[i*PARAM_COUNT+0],thread_params[i*PARAM_COUNT+1],thread_params[i*PARAM_COUNT+2],thread_params[i*PARAM_COUNT+3],thread_params[i*PARAM_COUNT+4],thread_params[i*PARAM_COUNT+5],thread_params[i*PARAM_COUNT+8],tune_global);
                    fclose(ftune);
                    //write current best assignment to file, including Xs, Xl + init time (so we can continue from there):
                    FILE *fs = fopen(ASSIGNMENT_FILE, "w");
                    for (int j=0; j<n+m*2; j++) fprintf(fs,"%.15f, ",v_best[i*(n+m*2)+j]); 
                    fprintf(fs,"%.15f ",t_init[i]);
                    fclose(fs);
                #endif
                
                    
            } else {
                if (!quiet) std::cout << "[" << energy_thread_min[i] << "] "; // printf("[%.2f] ",energy_thread_min[i]); 
            }
        }
        
        // output currently best params:
        if (!quiet) std::cout << " GLOBAL: " << TEXT_SILVER << tune_global << TEXT_DEFAULT << std::endl; //printf(" GLOBAL: %s%.2f%s\n",TEXT_SILVER, tune_global, TEXT_DEFAULT);
        
        // check for new best solution: ---------------------------------------------------------------------------------------
        int cnt_accepted = 0;
        for (int i=0; i<THREAD_COUNT; i++) {

            if (candidate_eval[i] < best_eval[i]) {
                //store new best point:
                best_eval[i] = energy_thread_min[i];
                //sat? we exit
                if (global_thread[i]==0) goto tuning_found_solution;
                
            } 
            // difference between candidate and current point:
            TFloat add = 1.0; 
            //int add = 0; if (tune_mode!=0) add = 1; //if we restart with v_best, we cant accept similar loc
            TFloat diff = candidate_eval[i] - curr_eval[i] + add;
            // temperature for current epoch:
            double t = INITIAL_TEMP / double(iter+1);
            // calculate metropolis acceptance criterion:
            double metropolis = exp(-(double)diff / t);
            // check if we should keep that point:
            if (diff < 0 || rand_r(generator)<metropolis) {
                if (!quiet) {
                    std::cout << TEXT_YELLOW << "c [TUNE] [" << i << "] " << num_change_vars[i] << " CANDIDATES ACCEPTED ";
                    std::cout << "(E=" << candidate_eval[i] << ", LOC=" << global_thread[i] << ", temp=" << t << TEXT_DEFAULT << std::endl;
                }
                curr_eval[i] = candidate_eval[i];
                cnt_accepted++;
                number_no_improvements[i] = 0;
            } else {
                for (int p=0; p<num_change_vars[i]; p++) {
                    thread_params[i*PARAM_COUNT+tune_param[i][p]] = rem_params[i*PARAM_COUNT+tune_param[i][p]];
                    if (!quiet) {
                        std::cout << "c [TUNE] ["<<i<<"] CANDIDATE REJECTED (E="<<candidate_eval[i]<<", LOC="<<global_thread[i]<<") p_"<<tune_param[i][p]<<" set back to ";
                        //printf("c [TUNE] [%d] CANDIDATE REJECTED (E=%.2f, LOC=%d) - p_%d set back to ", i,candidate_eval[i], global_thread[i], tune_param[i][p]);
                        std::cout << rem_params[i*PARAM_COUNT+tune_param[i][p]] << " temperature=" << t << std::endl;
                    }
                    rem_tune_param[i][p] = -1; //if not -1, we try it next time again...
                }
                number_no_improvements[i]++;
            }
        }

    }
    goto tuning_finished;
        printf("c [TUNE] FINISHED.\n");

    tuning_found_solution:
        printf(TEXT_YELLOW);
        printf("c [TUNE] SOLUTION FOUND.\n");

    tuning_finished:
        printf("c [TUNE] FINISHED.\n");

    return;    
}
//-------------------------------------------------------------------------------------------------------------------
//parse command line options // free: s
int scan_opt(int argc, char **argv, const char *opt) {
    char c;
    while ((c = getopt (argc, argv, opt)) != -1)
        switch (c) {
            case 'a': if (optarg[0]=='=') optarg[0] = ' '; timeout   =atof(optarg); break;
            case 'b': if (optarg[0]=='=') optarg[0] = ' '; dmm_beta=atof(optarg); break;
            case 'c': if (optarg[0]=='=') optarg[0] = ' '; dmm_alpha=atof(optarg); break;
            case 'd': if (optarg[0]=='=') optarg[0] = ' '; init_dt=atof(optarg); break;
            case 'e': if (optarg[0]=='=') optarg[0] = ' '; load_partable=atoi(optarg); break;
            case 'f': if (optarg[0]=='=') optarg[0] = ' '; forcebounds=atoi(optarg); break;
            case 'g': if (optarg[0]=='=') optarg[0] = ' '; tune_mode_params   =atoi(optarg); break;
            case 'h': if (optarg[0]=='=') optarg[0] = ' '; dmm_delta=atof(optarg); break;
            case 'i': if (optarg[0]=='=') optarg[0] = ' '; strcpy(input_filename, optarg); break;
            case 'j': if (optarg[0]=='=') optarg[0] = ' '; dmm_epsilon=atof(optarg); break;
            case 'k': if (optarg[0]=='=') optarg[0] = ' '; dmm_zeta=atof(optarg); break;
            case 'l': if (optarg[0]=='=') optarg[0] = ' '; seed=atoi(optarg); break;
            case 'm': if (optarg[0]=='=') optarg[0] = ' '; xl_max=atoi(optarg); break;
            case 'n': if (optarg[0]=='=') optarg[0] = ' '; dmm_gamma=atof(optarg); break;
            case 'o': if (optarg[0]=='=') optarg[0] = ' '; INTEGRATION_MODE=atoi(optarg); break;
            case 'p': if (optarg[0]=='=') optarg[0] = ' '; heuristics=atoi(optarg); break;
            case 'q': if (optarg[0]=='=') optarg[0] = ' '; quiet=atoi(optarg); break;
            case 'r': if (optarg[0]=='=') optarg[0] = ' '; tune_mode   =atoi(optarg); break;
            case 't': if (optarg[0]=='=') optarg[0] = ' '; tune=atoi(optarg); break;
            case 'u': if (optarg[0]=='=') optarg[0] = ' '; walltime_timeout   =atof(optarg); break;
            case 'v': if (optarg[0]=='=') optarg[0] = ' '; switchfraction   =atof(optarg); break;
            case 'w': if (optarg[0]=='=') optarg[0] = ' '; THREAD_COUNT=atoi(optarg); break;
            case 'x': if (optarg[0]=='=') optarg[0] = ' '; rk_errorrate_1=atof(optarg); break;
            case 'y': if (optarg[0]=='=') optarg[0] = ' '; rk_errorrate_2=atof(optarg); break;
            case 'z': if (optarg[0]=='=') optarg[0] = ' '; maxsteps=atof(optarg); break;
            case 's': if (optarg[0]=='=') optarg[0] = ' '; N_ITERATIONS=atof(optarg); break;
            
            default: return(-1);
        }
    return(0);
}

//-------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {

    int i, j;
    char buffer[512];
    global_best_thread = 0;

    char *syntax =
    "c          GENERAL OPTIONS:\n"
    "c          -i file      : INPUT FILE (CNF)\n"
    "c          -o [1;2;3;4] : ODE INTEGRATION: 1=CONSTANT, 3=RUNGE KUTTA, 4=IMPLICIT (DEFAULT:1)\n"
    "c          -f [0;1]     : APPLY BOUNDS AFTER EACH ODE STEP (DEFAULT: 1)\n"
    "c          -q [0;1]     : QUIET MODE; 0=OFF 1=ON (DEFAULT:0)\n"
    "c          -u [TFloat]  : WALL TIME TIMEOUT AFTER x s NO BETTER LOC FOUND (DEFAULT: INT_MAX)\n"
    "c          -w [int]     : NUMBER OF PARALLEL THREADS (DEFAULT: 8)\n"
    "c          -e [0;1]     : LOAD partable.txt (DEFAULT: 0)\n"
    "c\n"
    "c          -c [TFloat]  : ALPHA (GROWTH RATE FOR LONG TERM MEMORY Xl)\n"
    "c          -b [TFloat]  : BETA (GROWTH RATE FOR SHORT TERM MEMORY Xs)\n"
    "c          -n [TFloat]  : GAMMA (RESTRICTION FOR CM IN Xs)\n"
    "c          -h [TFloat]  : DELTA (GROWTH RATE FOR CM IN Xl)\n"
    "c          -j [TFloat]  : EPSILON (REMOVE SPURIOUS SOLUTION Xm=0)\n"
    "c          -k [TFloat]  : ZETA (REDUCTION FACTOR OF RIGIDITY G)\n"
    "c\n"
    "c          -l [int]     : RANDOM SEED; (DEFAULT [0]: 0.0;1.0;-1.0;rand [-1]: ALL 0.0; [x]: rand)\n"
    "c          -m [int]     : MAX VALUE FOR Xl (DEFAULT 10000)\n"
    "c\n"
    "c          -p [0;1]     : ALPHA HEURISTICS; 0=OFF 1=ON (DEFAULT:1)\n"
    "c\n"
    "c          TUNING OPTIONS:\n"
    "c          -t [0;1]     : TUNE CIRCUIT; 0=OFF 1=ON (DEFAULT:0)\n"
    "c          -g [0;1;2]   : PARAMS TO TUNE: 0=α..δ 1=ODE 2=ALL PARAMETERS (DEFAULT: 0)\n"
    "c          -r [0;1]     : TUNING MODE: 0=AWAYS FROM INIT; 1=CONTINOUS FROM BEST; 2=ALWAYS FROM INIT, BEST PARAMS; etc...\n"
    "c          -v [TFloat]  : SWITCH FRACTION (DEFAULT 0.0001)\n"
    "c          -s [Int]     : MAX NUMBER OF TUNING ITERATIONS (DEFAULT INT_MAX)\n"
    "c\n"
    "c          ODE - CONSTANT OPTIONS:\n"
    "c          -s [TFloat]  : STEP SIZE\n"
    "c          -a [TFloat]  : TIME OUT (INTEGRATION TIME)\n"
    "c\n"
    "c          ODE - RUNGE KUTTA OPTIONS:\n"
    "c          -x [TFloat]  : ERROR RATE 1 (ABSOLUTE ERROR TOLLERANCE)\n"
    "c          -y [TFloat]  : ERROR RATE 2 (RELATIVE ERROR TOLLERANCE)\n"
    "c          -d [TFloat]  : INITIAL dt VALUE\n"
    "c          -z [TFloat]  : MAXIMUM INTEGRATION STEPS (DEFAULT: INT_MAX)\n"
    ;
    goto on_continue;

on_break:
    printf("c Syntax: %s <... Args ...>\n", argv[0]);
    printf("c Args:\n");
    printf("%s", syntax);
    printf("\n");
    return EXIT_FAILURE;
    
on_continue:
    if(scan_opt(argc, argv, "r:v:w:u:q:p:b:c:n:h:j:k:l:m:e:f:g:s:a:z:x:y:d:o:t:i:")) goto on_break;

    //init standard radomizer:
    srand(seed);

    if (!quiet) printf("c --------------------------------------------\n");
    if (!quiet) printf("c TURINGX NEUROMORPHIC SAT-SOLVER (C) 2022    \n");
    if (!quiet) printf("c --------------------------------------------\n");
    if (!quiet) printf("c INSTANCE  : %s\n", input_filename);

    /// load CNF header:
    int ret;
    FILE *file = fopen(input_filename, "r");
    if (strcmp(buffer, "c") == 0) {
        while (strcmp(buffer, "\n") != 0) {
            ret = fscanf(file, "%s", buffer);
        }
    }
    while (strcmp(buffer, "p") != 0) {
        ret = fscanf(file, "%s", buffer);
        // parse for special keywords:
        if (strstr(buffer,"key=") != NULL) {
            strncpy(c_key, buffer + 4, 256-4);
            int c_len = strlen(c_key);
            printf("%sc CNF CONTAINS key=%s (%d BIT)%s\n",TEXT_CYAN,c_key,c_len,TEXT_DEFAULT);
            c_key_bin = (int *) calloc((size_t) c_len, sizeof(int));
            for (int i=0; i<c_len; i++) {
                if (c_key[i]=='0') c_key_bin[i] = 0;
                if (c_key[i]=='1') c_key_bin[i] = 1;
            }
            c_key_loaded = true;
        }
        if (strstr(buffer,"plain=") != NULL) {
            strncpy(c_plain, buffer + 6, 256-6);
            int c_len = strlen(c_plain);
            printf("%sc CNF CONTAINS plain=%s (%d BIT)%s\n",TEXT_CYAN,c_plain,c_len,TEXT_DEFAULT);
            c_plain_bin = (int *) calloc((size_t) c_len, sizeof(int));
            for (int i=0; i<c_len; i++) {
                if (c_plain[i]=='0') c_plain_bin[i] = 0;
                if (c_plain[i]=='1') c_plain_bin[i] = 1;
            }
            c_plain_loaded = true;
        }
        if (strstr(buffer,"cipher=") != NULL) {
            strncpy(c_cipher, buffer + 7, 256-7);
            int c_len = strlen(c_cipher);
            printf("%sc CNF CONTAINS cipher=%s (%d BIT)%s\n",TEXT_CYAN,c_cipher,c_len,TEXT_DEFAULT);
            c_cipher_bin = (int *) calloc((size_t) c_len, sizeof(int));
            for (int i=0; i<c_len; i++) {
                if (c_cipher[i]=='0') c_cipher_bin[i] = 0;
                if (c_cipher[i]=='1') c_cipher_bin[i] = 1;
            }
            c_cipher_loaded = true;
        }
        //--
    }
    ret = fscanf(file, " cnf %i %i", &n, &m);

    if (!quiet) printf("c VARIABLES : %'d\n", n);
    if (!quiet) printf("c CLAUSES   : %'d\n", m);
    if (!quiet) printf("c RATIO     : %lf\n", (double) m / n);

    xl_max = xl_max * m;
    if (xl_max<=0) xl_max = INT_MAX;

    /// reserve  memory - needs to be done before anything else:
    cls = (int *) calloc((size_t) m*MAX_LITS_SYSTEM, sizeof(int));
    for (int i=0; i<m*MAX_LITS_SYSTEM; i++) cls[i] = 0; 
    numOccurrenceT = (int *) calloc((size_t) n+1, sizeof(int));
    clauseSizes = (int *) calloc((size_t) m, sizeof(int));
        
    /// read CNF: /////////////////////////////////////////
    int lit;
    for (i = 0; i < m; i++) {
        j = 0; 
        do {
            ret = fscanf(file, "%s", buffer);
            if (strcmp(buffer, "c") == 0) {
                ret = fscanf(file, "%512[^\n]\n", buffer);
                continue;
            }
            lit = atoi(buffer);
            cls[i*MAX_LITS_SYSTEM+j] = lit;
            // increase number of Occurence of the variable, max number of occurences
            if (lit!=0) {
                numOccurrenceT[abs(lit)]++;
                if (numOccurrenceT[abs(lit)]>maxNumOccurences) {maxNumOccurences=numOccurrenceT[abs(lit)];}
                clauseSizes[i] = j+1;
            }
            j++;
        } while (strcmp(buffer, "0") != 0);
        j--;
        if (j > MAX_LITS_SYSTEM) {
            printf("c ERROR: CLAUSE %d HAS MORE THAN %d LITERALS.\n",i,MAX_LITS_SYSTEM);
            return EXIT_FAILURE;
        }
    }
    
    if (!quiet) printf("c MAX VARIABLE OCCURENCE: %'d\n", maxNumOccurences);

    if (!quiet) printf("c FIRST 10 CLAUSES:\n");
    for (i = 0; i < 11; i++) {
        if (!quiet) printf("c CLAUSE %i: ",i);
        for (j = 0; j < clauseSizes[i]; j++) {if (!quiet) printf(" %d",cls[i*MAX_LITS_SYSTEM+j]);}
        if (!quiet) printf(" (%d)",clauseSizes[i]);
        if (!quiet) printf("\n");
    }

    if (!quiet) printf("c LAST 10 CLAUSES:\n");
    for (i = m-1; i > (m-10); i--) {
        if (!quiet) printf("c CLAUSE %i: ",i);
        for (j = 0; j < clauseSizes[i]; j++) {if (!quiet) printf(" %d",cls[i*MAX_LITS_SYSTEM+j]);}
        if (!quiet) printf(" (%d)",clauseSizes[i]);
        if (!quiet) printf("\n");
    }

    //build occurence array: [var][cls...] /////////////////////////////////////////
    occurrence = (int *) calloc((size_t) (n+1)*maxNumOccurences, sizeof(int));
    occurenceCounter = (int *) calloc((size_t) n+1, sizeof(int));
    
    for (i=0; i<m; i++) {
        for (j = 0; j < clauseSizes[i]; j++) {
            lit = abs(cls[i*MAX_LITS_SYSTEM+j]);
            occurrence[lit*maxNumOccurences+occurenceCounter[lit]] = i;
            occurenceCounter[lit]++;
        }
    }

    //output:
    printf("c OCCURENCES: ");
    for (i=1; i<=20; i++) printf("%d->%d ",i,occurenceCounter[i]);
    printf("\n");

    /// initialize arrays:
    v_best = (TFloat *) calloc((size_t) THREAD_COUNT*(n+m*2), sizeof(TFloat));
    loc_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    unit_clause_vars = (int *) calloc((size_t) n+1, sizeof(int));
    energy_thread = (TFloat *) calloc((size_t) THREAD_COUNT, sizeof(TFloat));
    energy_thread_min = (TFloat *) calloc((size_t) THREAD_COUNT, sizeof(TFloat));
    global_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    global_all_runs_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    time_thread = (TFloat *) calloc((size_t) THREAD_COUNT, sizeof(TFloat));
    time_thread_actual = (TFloat *) calloc((size_t) THREAD_COUNT, sizeof(TFloat));
    walltime_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    t_init = (TFloat *) calloc((size_t) THREAD_COUNT, sizeof(TFloat));
    t_begin_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    t_end_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    stepcounter = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    initial_assignments = (TFloat *) calloc((size_t) THREAD_COUNT*(n+m*2), sizeof(TFloat));
    thread_params = (TFloat *) calloc((size_t) THREAD_COUNT*PARAM_COUNT, sizeof(TFloat));
    partable = (double *) calloc((size_t) PARAM_COUNT*4, sizeof(double));
    defaults = (double *) calloc((size_t) PARAM_COUNT*128, sizeof(double));
    dmm_alpha_cm = (TFloat *) calloc((size_t) m*THREAD_COUNT, sizeof(TFloat));
    if (c_key_loaded) {
        int keylen = strlen(c_key);
        c_key_thread = (int *) calloc((size_t) THREAD_COUNT*keylen, sizeof(int));
    }

    x_rem = (int *) calloc((size_t) n, sizeof(int));

    if (!quiet) printf("c ARRAYS INITIALISED\n");

    /// detect (and assign) unit clauses (clauses with one literal): ----------------------------
    int cnt_unit_clauses = 0;
    for (int i=0; i<m; i++) {
        if (clauseSizes[i]==1) {
            cnt_unit_clauses++;
            unit_clauses = true;
            int lit = cls[i*MAX_LITS_SYSTEM];
            if (lit>0) unit_clause_vars[abs(lit)] = 1;
            if (lit<0) unit_clause_vars[abs(lit)] = -1;
        }
    }
    if (!quiet) printf("c %d UNIT CLAUSES DETECTED.\n",cnt_unit_clauses);

    /// set _all_runs vars: ---------------------------------------------------------------------
    for (int i=0; i<THREAD_COUNT; i++) global_all_runs_thread[i] = m;

    /// set t_init: ---------------------------------------------------------------------
    for (int i=0; i<THREAD_COUNT; i++) t_init[i] = 0.0;

    /// prepare SOLUTION_FILE:
    strcpy(SOLUTION_FILE, input_filename);
    char APPENDS[128]; strcpy(APPENDS,".solution.txt");
    strcat(SOLUTION_FILE,APPENDS);

    /// load solution: --------------------------------------------------------------------------
    bool solution_loaded = false;
    if (load_solution==true) {
        FILE *fs;
        fs = fopen(SOLUTION_FILE, "r");
        if (fs == NULL) {
            if (!quiet) printf("c %s NOT PROVIDED. USING RANDOM ASSIGNMENTS.\n",SOLUTION_FILE);
        } else {
            // set everything to 0:
            for (int i=0; i<n; i++) v_best[i] = 0.0;
            // now load solution:
            float ivar;
            for (int i=0; i<n;i++) {
            //for (int i=0; i<13000;i++) {
                fscanf(fs,"%f, ",&ivar);
                if (ivar>0) v_best[i] = 1.0;
                if (ivar<0) v_best[i] =-1.0;
                
            }
            fclose(fs);

            solution_loaded = true;
            if (!quiet) {
                printf("c SOLUTION %s LOADED.\n",SOLUTION_FILE);
                printf("c ASSIGNMENT: "); for (int i=0; i<50; i++) std::cout<<v_best[i]<<" "; std::cout<<std::endl; // printf("%.15f ",v_best[i]); printf("\n");
            }
        }
    } else {
        if (!quiet) printf("c SOLUTION %s NOT LOADED (load_solution=false)\n",SOLUTION_FILE);
        for (int i = 0; i<n+m*2; i++) v_best[i] = TFloat(0.00);
    }

    /// prepare ASSIGNMENT_FILE:
    strcpy(ASSIGNMENT_FILE, input_filename);
    char APPENDA[128]; strcpy(APPENDA,".assignment.txt");
    strcat(ASSIGNMENT_FILE,APPENDA);

    /// load assignment: ------------------------------------------------------------------------
    bool assignment_loaded = false;
    if (load_assignment) {
        FILE *fa;
        fa = fopen(ASSIGNMENT_FILE, "r");
        if (fa!= NULL) {
            float ivar;
            for (int i=0; i<n+m*2;i++) {
                fscanf(fa,"%f, ",&ivar);
                v_best[i] = ivar;
            }
            fscanf(fa,"%f",&ivar);
            for (int i=0; i<THREAD_COUNT; i++) t_init[i] = ivar;
            fclose(fa);
            assignment_loaded = true;
            if (!quiet) {
                printf("c ASSIGNMENT %s LOADED, WILL START AT t=",ASSIGNMENT_FILE);
                std::cout << t_init[0] << std::endl;
                printf("c "); for (int i=0; i<10; i++) std::cout<<v_best[i]<<" "; std::cout<<std::endl; // printf("%.15f ",v_best[i]); printf("\n");
            }
        }   
    }

    /// load partable: --------------------------------------------------------------------------
    if (load_partable) {
        FILE *fp;
        fp = fopen(PARTABLE_FILE, "r");
        if (fp == NULL) {
            if (!quiet) printf("c partable.txt NOT PROVIDED. USING DEFAULT SETTINGS.\n");
        } else {
            // first line: number of default values:
            int num_defaults;
            fscanf(fp,"%d",&num_defaults);
            // num_defaults x params:
            double ipar;
            for (int j=0; j<num_defaults; j++) {
                for (int i=0; i<PARAM_COUNT; i++) {
                    fscanf(fp,"%lf;",&ipar);
                    defaults[j*PARAM_COUNT+i] = ipar;
                    // the #0 default will be set here, too:
                    if (j==0) {
                        switch (i) {
                            case 0: dmm_alpha = ipar; break;
                            case 1: dmm_beta = ipar; break;
                            case 2: dmm_gamma = ipar; break;
                            case 3: dmm_delta = ipar; break;
                            case 4: dmm_epsilon = ipar; break;
                            case 5: dmm_zeta = ipar; break;
                            case 6: rk_errorrate_1 = ipar; break;
                            case 7: rk_errorrate_2 = ipar; break;
                            case 8: init_dt = ipar; break;
                        }
                    }
                }
            }
            //now bounds for the vars:
            for (int i=0; i<PARAM_COUNT; i++) {
                double softlb; fscanf(fp,"%lf;",&softlb); partable[i*4+0] = softlb;
                double softub; fscanf(fp,"%lf;",&softub); partable[i*4+1] = softub;
                double hardlb; fscanf(fp,"%lf;",&hardlb); partable[i*4+2] = hardlb;
                double hardub; fscanf(fp,"%lf;",&hardub); partable[i*4+3] = hardub;
                //printf("c P%d: soft: %.10f-%.10f hard: %.10f-%.10f\n",i,softlb,softub,hardlb,hardub);
            }
            fclose(fp);
            if (!quiet) printf("c PARTABLE LOADED %d PARAMS AND %d DEFAULTS.\n",PARAM_COUNT,num_defaults);
            for (int i=0; i<num_defaults; i++) {
                printf("c DEFAULTS SET %d: ",i);
                for (int j=0; j<PARAM_COUNT; j++) printf("%.32f ",defaults[i*PARAM_COUNT+j]);
                printf("\n");
            }
            partable_loaded = true;
        }
    } else {
        if (!quiet) printf("c partable.txt NOT LOADED.\n");
    }

    /// OUPUT SETTINGS: ---------------------------------------------------------------------------
    if (!quiet) {
        printf(TEXT_CYAN);
        printf("c SETTINGS:\n");
        printf("c #THREADS        : %d\n",THREAD_COUNT);
        if (INTEGRATION_MODE==ODE_RUNGEKUTTA) {
            printf("c INTGRATION MODE : ADAPTIVE RUNGE KUTTA\n");
            std::cout << std::setprecision(digits) << std::fixed;
            printf("c ERROR RATE 1    : "); std::cout << rk_errorrate_1 << std::endl;
            printf("c ERROR RATE 2    : "); std::cout << rk_errorrate_2 << std::endl;
            printf("c INITAL dt       : "); std::cout << init_dt << std::endl;
            printf("c MAX STEPS       : "); std::cout << maxsteps << std::endl;
        }
        if (INTEGRATION_MODE==ODE_CONSTANT) {
            std::cout << std::setprecision(digits) << std::fixed;
            printf("c INTGRATION MODE : CONSTANT\n");
            //printf("c STEPSIZE        : "); std::cout << stepsize << std::endl;
            printf("c TIMEOUT         : "); std::cout << timeout << std::endl;
            printf("c INITAL dt       : "); std::cout << init_dt << std::endl;
            printf("c MAX STEPS       : "); std::cout << maxsteps << std::endl;
        }
        printf("c HEURISTICS      : %d\n",heuristics);
        printf("c TUNE CIRCUIT    : %d\n",tune);
        std::cout << std::setprecision(digits) << std::fixed;
        printf("c ALPHA           : "); std::cout << dmm_alpha << std::endl; // %.17f\n",dmm_alpha);
        printf("c BETA            : "); std::cout << dmm_beta << std::endl;
        printf("c GAMMA           : "); std::cout << dmm_gamma << std::endl;
        printf("c DELTA           : "); std::cout << dmm_delta << std::endl;
        printf("c EPSILON         : "); std::cout << dmm_epsilon << std::endl;
        printf("c ZETA            : "); std::cout << dmm_zeta << std::endl;
        printf("c SEED            : %d\n",seed);
        printf("c XL_MAX          : %.d\n",xl_max);

        printf(TEXT_DEFAULT);
    }
    

    /// init thread specific ODE parameters //////////////////////
    
    /* RNG */
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double>  rand_v_double(-1.0, 1.0); // voltage continous between -1 and +1 
    std::uniform_int_distribution<int>      rand_v(-1, 1); // voltage -1, 0 or +1
    std::uniform_real_distribution<double>  rand_Xs(0.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xl(1.0, 100.0); //xl_max); //10.00

    //TFloat init_dt_matrix[8] = {init_dt, 0.9, 0.8, 0.7, 0.6, 0.15, 0.10, 0.05};
    
    for(int i = 0; i < THREAD_COUNT; i++ ) {
        
        //params:
        thread_params[i*PARAM_COUNT+0]    = dmm_alpha;
        thread_params[i*PARAM_COUNT+1]    = dmm_beta;
        thread_params[i*PARAM_COUNT+2]    = dmm_gamma;
        thread_params[i*PARAM_COUNT+3]    = dmm_delta;
        thread_params[i*PARAM_COUNT+4]    = dmm_epsilon;
        thread_params[i*PARAM_COUNT+5]    = dmm_zeta;
        thread_params[i*PARAM_COUNT+6]    = rk_errorrate_1;
        thread_params[i*PARAM_COUNT+7]    = rk_errorrate_2;
        thread_params[i*PARAM_COUNT+8]    = init_dt; //init_dt_matrix[i]; 

        // set declining alpha params for each threads: ---- TEST ------
        if (massive==0 && i>0) thread_params[i*PARAM_COUNT+0] = thread_params[i*PARAM_COUNT+0] / i;
        // ---------
        
        //initial assignments:
        if (i<=2 && seed==0) {
            //special case: seed 0: all V zeros:
            if (i==0) for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = TFloat(0.0);
            if (i==1) for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = TFloat(1.0);
            if (i==2) for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = TFloat(-1.0);
        } else {
            int set_seed = seed + (i*4096);
            for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = TFloat(rand_v(generator));
            for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] =  0.0; //TFloat(0.0); //rand_Xs(generator); //; 
            for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = 1.0; //TFloat(1.0); //rand_Xl(generator); // //TFloat(1.0); //
        }
        // seed -1 ? all to zero
        if (seed==-1) {
            for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = TFloat(0.0);
        }
        // seed -2 ? all to -1
        if (seed==-2) {
            for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = TFloat(-1.0);
        }
        // solution loaded: Assign V, set Xs,Xv default:
        if (solution_loaded) {
            for (int j=0; j<n+m*2; j++) {
                initial_assignments[i*(n+m*2)+j] = v_best[j];
            }
            for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] = rand_Xs(generator);
            for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = rand_Xl(generator);
        }
        // assignments loaded: Assign V,Xs and Xl:
        if (tune==1 && assignment_loaded) {
            for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = v_best[j];
            for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] = v_best[j];
            for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = v_best[j];
        }
        // override unit clauses for initial assignments:
        if (unit_clauses) {
            for (int j=0; j<n; j++) {
                if (unit_clause_vars[j+1]!=0) initial_assignments[i*(n+m*2)+j] = unit_clause_vars[j+1];
            }
        }
    }

    if (tune) writelocfile=false;

    /// prepare LOC file?
    if (writelocfile) {
        strcpy(LOC_FILE, input_filename);
        char APPEND[128]; strcpy(APPEND, ".loc.txt");
        strcat(LOC_FILE,APPEND);
        #ifdef USEFPRINTF
            FILE *floc = fopen(LOC_FILE, "w");
            fprintf(floc,"INSTANCE;WALLTIME;ODE_TIME;LOC;ENERGY\n");
            fprintf(floc,"%d;%.5f;%.5f;%d;%d\n",0,0.0,0.0,m,m);
            fclose(floc);
        #endif
        
    }

    /// prepare TUNING_FILE:
    if (tune) {
        strcpy(TUNING_FILE, input_filename);
        char APPENDT[128]; strcpy(APPENDT,".tuning.txt");
        strcat(TUNING_FILE,APPENDT);
        #ifdef USEFPRINTF
            FILE *ftune = fopen(TUNING_FILE, "w");
            fprintf(ftune,"ITERATION;INSTANCE;PARAM_01;PARAM_02;PARAM_03;PARAM_04;PARAM_05;PARAM_06;PARAM_07;ENERGY\n");
            fclose(ftune);
        #endif
        
    }

    /* ------------------------------------- SOLVER -------------------------------------------------------------- */
    /// tune circuit? ////////////////////////////////////////////////////////////////////////////////////////////
    if (tune) {
        if (!quiet) printf("c [TUNE] TUNING CIRCUIT...\n");
        tune_circuit();
    } else {
        /// start solver (multi-thread) //////////////////////////////////////////////////////////////////////////
        if (!massive) {
                // normal mode, 1 run: ---------------------------------------------------------------------------
                pthread_t threads[THREAD_COUNT];
                for(uint64_t i = 0; i < THREAD_COUNT; i++ ) {
                    pthread_create(&threads[i], NULL, apply, (void *)i);
                  /*int rc = pthread_create(&threads[i], NULL, apply, (void *)i);
                  if (rc) {
                     cout << "Error:unable to create thread," << rc << endl;
                     exit(-1);
                  }*/
                  running_threads++;
                }
                // display loop while threads are running:
                while (running_threads > 0) {
                    if (quiet) {
                        
                        printf("\rc %d THREADS: ",running_threads);
                        for (int i=0; i<THREAD_COUNT; i++) {
                            std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT;
                            std::cout << time_thread_actual[i] << " [" << global_thread[i] << "] ";
                        }
                        fflush(stdout);
                        
                    } else {
                        // special case: key=xxx in CNF:
                        if (c_key_loaded) {
                            int keylen = strlen(c_key);
                            // run through all threads:
                            key_best_thread = -1;
                            c_key_best = 0;
                            int correct;
                            for (int t=0; t<THREAD_COUNT; t++) {
                                correct = 0;
                                for (int i=0; i<keylen; i++) {
                                    if (c_key_thread[t*keylen+i]==c_key_bin[i]) correct++;
                                }
                                if (correct>c_key_best) {
                                    c_key_best = correct;
                                    key_best_thread = t;
                                }
                            }
                            // screen output:
                            std::cout << "\rc ";
                            for (int i=0; i<keylen; i++) {
                                int _assign = c_key_thread[key_best_thread*keylen+i];
                                if (_assign==c_key_bin[i]) {
                                    std::cout << _assign; 
                                } else {
                                    std::cout << TEXT_RED << _assign << TEXT_DEFAULT;
                                }
                                std::cout << " ";
                            }
                            std::cout << " => " << c_key_best << "/" << strlen(c_key) <<" ([" << key_best_thread << "] TOP:"<< c_key_best_all_runs <<") ";
                            //std::cout << std::endl;
                            if (c_key_best>c_key_best_all_runs) {
                                c_key_best_all_runs = c_key_best;
                                std::cout << std::endl;
                            }
                            fflush(stdout);
                            
                        }
                        //---special case
                    } //---!quiet
                } //---running_threads
        } else {
            // run in massive mode: ---------------------------------------------------------------------------------
            // needs -u 5 to be set, otherwise runs until MAX_INT second no improvement
            printf("c MASSIVE MODE\n");
            solved = false;
            massive_global = m;
            massive_energy = INT_MAX;
            int massive_iterations = seed;
            int massive_best_seed = 0;

            while (!solved) {
                //set seed:
                srand(massive_iterations*256);
                //std::mt19937 generator_massive(massive_iterations*256);
                //set random Xs, Xl:
                for(int i = 0; i < THREAD_COUNT; i++ ) {
                    
                    //set xx vars random:
                    for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = 0.0;
                    for (int j=0; j<1000; j++) {
                        int rr = rand() % 2;
                        int rrvar = rand() % n;
                        if (rr==0) initial_assignments[i*(n+m*2)+rrvar] = -1;
                        if (rr==1) initial_assignments[i*(n+m*2)+rrvar] = +1;
                    }

                }
                printf("%sc [MASSIVE] ITERATION %d - SEED %d - NEW RANDOM ASSIGNMENTS SET.%s\n",TEXT_YELLOW,massive_iterations,massive_iterations*256, TEXT_DEFAULT);


                //walltime_timeout = 5.0; //interrupt after xx s
                pthread_t threads[THREAD_COUNT];
                for(uint64_t i = 0; i < THREAD_COUNT; i++ ) {
                    pthread_create(&threads[i], NULL, apply, (void *)i);
                    running_threads++;
                }
                while (running_threads > 0) {
                    if (quiet) {
                        printf("\rc %d THREADS: ",running_threads);
                        for (int i=0; i<THREAD_COUNT; i++) {
                            std::cout << std::setprecision(2) << std::fixed << TEXT_DEFAULT;
                            std::cout << time_thread_actual[i] << " [" << global_thread[i] << "] ";
                        }
                        std::cout << "GLOBAL:" << global << " [" << massive_global << "] best sol: " << (n-solmax) << "   ";
                        fflush(stdout);
                    }
                }
                if (global<massive_global) {
                    massive_global = global;
                    massive_best_seed = massive_iterations*256;
                    printf("\n%sc [MASSIVE] ITERATION %d - NEW GLOBAL OPTIMUM %d WITH SEED %d%s\n",TEXT_YELLOW,massive_iterations,massive_global,massive_best_seed,TEXT_DEFAULT);
                } else {
                    printf("\nc [MASSIVE] ITERATION %d - CURRENT BEST GLOBAL=%d WITH SEED %d\n",massive_iterations, massive_global, massive_best_seed);
                }
                massive_iterations++;

            }
        } // end run in massive mode: ----------------------------------------------------------
        

        

        printf("%d",global); // output lowest global for google advisor param optimizer
        pthread_exit(NULL);
    }
   
    // free memory:
    free(cls);
    free(v_best);
    free(occurrence);          //occurences of each lit
    free(occurenceCounter);    //counter for allocation
    free(numOccurrenceT);      //number of occurence of a lit
    free(clauseSizes); 
    free(thread_params);
    free(initial_assignments);
    free(loc_thread);
    free(global_thread);
    free(time_thread);
    free(time_thread_actual);
    free(unit_clause_vars);

    
    printf("c EXIT\n");
    return EXIT_SUCCESS;
}

