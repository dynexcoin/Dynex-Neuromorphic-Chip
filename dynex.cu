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

// COMPILE w. CUDA:      nvcc dynex.cu -o dynex_gpu -std=c++17 -O4 

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
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
#include "boost/date_time/posix_time/posix_time.hpp"

void     INThandler(int);
bool     apply_restart = false;

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
#define         ODE_LOGFILE         "log.csv"
#define         TUNE_LOGFILE        "tuninglog.csv"
#define         PARTABLE_FILE       "partable.txt"
#define         SOLUTION_FILE       "solution.txt"
#define         FLOWVECTOR_FILE     "flowvector.csv"

int THREAD_COUNT = 8;

//-------------------------------------------------------------------------------------------------------------------
// precision (attention: change from double needs to update AtomicAddd, too!)
typedef double TFloat; //precision of cuda vars
//typedef float TFloat; //precision of cuda vars

typedef TFloat value_type;
typedef std::vector< value_type > state_type; //we use a vector which has dynamic size

//--------------------------------------------------------------------------------------------------------------
/* SETTINGS & PARAMETERS */
bool   quiet          = false; //no screen output
bool   flowvector_log = false; //output entire FLOWVECTOR_FILE (only thread 0)
bool   load_partable  = true;  //load partable.txt
bool   writelogfile   = false; //write logfile ODE_LOGFILE
bool   coupledsystem  = false; //create a coupled circuit  !!! TODO: NEEDS TO BE FIXED AND CHECKED
bool   loadsolution   = false; 

// Equations - constants:
__device__ TFloat vmin = -1.0;  // voltage lower bound - solved by ./dmm_gpu -i red_8bit_10rds_cut3.cnf -s 0.025 -c 1.25 -k 0.0125 (sometimes)
__device__ TFloat vmax =  1.0;  // voltage upper bound

double dmm_alpha      = 5.0;    // -c  growth rate for long term memory Xl
double dmm_beta       = 20.0;   // -b  growth rate for short term memory Xs
double dmm_gamma      = 0.25;   // -n  restriction for Cm in short term memory Xs
double dmm_delta      = 0.05;   // -h  restriction for Cm in long term memory Xl
double dmm_epsilon    = 0.1;    // -j  remove spurious solution X,s,m = 0
double dmm_zeta       = 0.1;    // -k  reduction factor of rigidity G (learning rate) 10^-3; for ratio>=6: 10^-1 (0.1)
int    seed           = 1;      // -l  random seed value (for initial assignment)
int    xl_max         = 10000;  // -m  10^4 M (x count clauses will be applied automatically - ODE should NEVER reach this value)

// ODE settings:
int INTEGRATION_MODE = ODE_CONSTANT;

// Constant integration params:
double stepsize             = 0.15; //15; // 0.15; //0.015; //0.15; // 0.0078125;        //2^-7
double timeout              = INT_MAX; //max simulated time; stops at reaching it
double walltime_timeout     = INT_MAX; 
double walltime_abs_timeout = INT_MAX;

// Runge-Kutta Adaptive params:
double rk_errorrate_1       = 0.52; // 1.0e-5; //both 0.1 for CBS_k3_n100_m403_b10_3.cnf 
double rk_errorrate_2       = 0.10;
double init_dt              = 0.0078125;; //2^-7
double maxsteps             = INT_MAX; 

// tuneing options:
bool   tune                 = false;
double switchfraction       = 0.0001; 
int    tune_mode            = 0; // 0 = always from -1 assignments; 1 = continous 
int    tune_mode_params     = 2; // 0 = alpha..zeta, 1=ODE params, 2=all params
int    tune_global;

// custom adaptive (experiemental) params:
double adaptive_min         = 0.0078125; //2^-7   0.0078125
double adaptive_max         = 1000; //1000;        //10^3
double step_error           = 5.5;//5.5; //0.001;

// Additional ODE heuristics:
bool         heuristics     = true;
double       alpha_increase = 1.1; //1.1; //increase factor for dmm_alpha_m
double       alpha_decrease = 0.9; //0.9; //decrease factor for dmm_alpha_m
int          alpha_correction = 10000; //10000; //reallocation of dmm_alpha_m every alpha_correction steps
int          alpha_resetvalue = 5.00; //1.00; //after maximum reached, this will be the new alpha

//-------------------------------------------------------------------------------------------------------------------
/* VARIABLES */
bool            debug = false;
void            INThandler(int);       // ctrl-c handler
#define         max_lits_system        10 //10   //3;   //max k-SAT

int             * cls;                 //stores the clauses (max_lits_system columns)
int             * occurrence;          //occurences of each litint             * occurenceCounter;    //counter for allocation
int             * numOccurrenceT;      //number of occurence of a lit
int             * clauseSizes;         //number of literals of a clause
int             * occurenceCounter; 

int             maxNumOccurences = 0;  //max #occurence in clauses of a var
int             n;                     //number of variables
int             m;                     //number of clauses
int             solved;                //is formula solved? UNSAT = 0; SAT = 1;
int             global;                //current lowest loc over all threads
int             global_all_runs;       //best global over all runs (while tuning or restarting f.e.)
TFloat          * v_best;              //G_field assignment of global local minima for each thread
char            input_filename[256];

/* thread specific vars */
int             * loc_thread;           //current local minima of thread (currently)
int             * global_thread;        //global local minima of thread
int             * global_all_runs_thread; //global local minimum of thread over all runs
double          * time_thread;          //simtime of thread (better loc)
double          * time_thread_actual;   //simtime of thread (currently)
double          * walltime_thread;      //walltime of thread (better loc)
TFloat          * initial_assignments;  //initial assignment for each thread
double          * thread_params;        //parameters for each thread
double          * t_begin_thread;       //starting time of thread
double          * t_end_thread;         //current time of thread
int             global_best_thread;     //thread# which has currently best global;

double          * partable;             //if partable.txt is provided, here are the bounds
double          * defaults;             //if partable.txt is provided, here are the default values (max 128 threads)
bool            partable_loaded = false; 

int             stepcounter;

struct node {
    int id;                 //thread-id
    int *model;             //current assignment
    int *temporal;          //temp assignment for oracle (not used in production)
    int *optimal;           //best assignment and solution afterwards
};

double t_begin;
double t_end; 
double t_abs_begin;
double t_abs_end;

TFloat t_rem;
//-------------------------------------------------------------------------------------------------------------------
/* CUDA vars */
TFloat * d_x;// state variables (V,Xs,Xl) 
TFloat *h_v; 
TFloat * d_x_tmp; // temporary state variables (V,Xs,Xl) for adaptive step
TFloat * d_dxdt; // derivative
int * d_cls;
int * d_clauseSizes;
int * d_loc; 
int h_loc[1];
int * d_varchanges;
int h_varchanges[1];
float * d_energy; 
float h_energy[1];
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


//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------
//                                                    CUDA                                                         //
//-------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------

//
#define gpuErrorCheck(ans, abort) { gpuAssert((ans), __FILE__, __LINE__, abort); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if(code != cudaSuccess) {
        fprintf(stderr,"assert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if(abort) {
        exit(code);
        }
    }
}
//
// GPU Helper function: atomicAdd for double:
__device__ double atomicAddd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, 
                        __double_as_longlong(val + 
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

//-------------------------------------------------------------------------------------------------------------------
// -o 1 forward Euler
// GPU kernel: x=>dxdt=>x 
__global__
void gpu_euler(TFloat * d_x_tmp, TFloat * d_x, TFloat * d_dxdt, int size, TFloat t, double h, int n, int m, int xl_max, int * d_varchanges) {
    
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < size; 
         i += blockDim.x * gridDim.x) 
      {

    	// euler step:
        d_x_tmp[i] = d_x[i] + h * d_dxdt[i]; 
        
        // bounded variables:
        if (i<n) {
            if (d_x_tmp[i]<-1.0) d_x_tmp[i] = -1.0;
            if (d_x_tmp[i]> 1.0) d_x_tmp[i] =  1.0;
        }
        if (i>=n && i<n+m) {
            if (d_x_tmp[i]<0.0) d_x_tmp[i] = 0.0;
            if (d_x_tmp[i]>1.0) d_x_tmp[i] = 1.0;   
        }
        if (i>=n+m*2) {
            if (d_x_tmp[i]<1.0) d_x_tmp[i] = 1.0;
            if (d_x_tmp[i]>xl_max) d_x_tmp[i] = xl_max;   
        }
        
	// change of var > x? increase counter d_varchanges
        if (fabs(d_x[i]-d_x_tmp[i])>=1.0) {
      		atomicAdd(&d_varchanges[0], 1); 
      	}
    }
    
}
//-------------------------------------------------------------------------------------------------------------------
// GPU kernel: reset d_loc (set to m), d_energy and d_varchanges
__global__
void gpu_reset_vars(int * d_loc, float * d_energy, int * d_varchanges, int m) {
    d_loc[0] = m;
    d_energy[0] = 0.0;
    d_varchanges[0] = 0;
}
//-------------------------------------------------------------------------------------------------------------------
// GPU kernel: reset dxdf, set all values to 0.00
__global__
void gpu_reset_dxdt(TFloat * d_dxdt, int size) {
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < size; 
         i += blockDim.x * gridDim.x) 
      {
        d_dxdt[i] = 0.0;
      }
}

//-------------------------------------------------------------------------------------------------------------------
// GPU kernel: reset x, set all values to 0.00
__global__
void gpu_reset_x(TFloat * d_x, int size) {
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < size; 
         i += blockDim.x * gridDim.x) 
      {
        d_x[i] = 0.0;
      }
}

//-------------------------------------------------------------------------------------------------------------------
// GPU kernel: do a step => generate dxdt & update d_loc
__global__
void gpu_step(TFloat * d_x, TFloat * d_dxdt, int * d_cls, int * d_loc, float * d_energy, int m, int n, double xl_max, double m_alpha, double m_beta, double m_gamma, double m_delta, double m_epsilon, double m_zeta, int * d_clauseSizes) {
    for (int clause = blockIdx.x * blockDim.x + threadIdx.x; 
         clause < m; 
         clause += blockDim.x * gridDim.x) 
      {
        
        TFloat Qi = (d_cls[clause*max_lits_system+0]>0)? 1.0:-1.0; // +1 if literal is >0, otherwise -1
        TFloat Qj = (d_cls[clause*max_lits_system+1]>0)? 1.0:-1.0; // +1 if literal is >0, otherwise -1
        TFloat Qk = (d_cls[clause*max_lits_system+2]>0)? 1.0:-1.0; // +1 if literal is >0, otherwise -1
        TFloat C;
        TFloat Xs = d_x[clause+n]; if (Xs<0.0) Xs=0.0; if (Xs>1.0) Xs=1.0; //boundary for xs € [0,1]:
        TFloat Xl = d_x[clause+n+m]; if (Xl<1.0) Xl=1.0; if (Xl>xl_max) Xl=xl_max; //boundary for xl € [1,10⁴M]:

        //k-sat implementation:
        int k = d_clauseSizes[clause];
        
        if (k==3) {
            int liti = abs(d_cls[clause*max_lits_system+0]);
            int litj = abs(d_cls[clause*max_lits_system+1]);
            int litk = abs(d_cls[clause*max_lits_system+2]);
            TFloat Vi = d_x[liti-1]; if (Vi<vmin) Vi=vmin; if (Vi>vmax) Vi=vmax; //boundary for v € [-1,1]:
            TFloat Vj = d_x[litj-1]; if (Vj<vmin) Vj=vmin; if (Vj>vmax) Vj=vmax; //boundary for v € [-1,1]:
            TFloat Vk = d_x[litk-1]; if (Vk<vmin) Vk=vmin; if (Vk>vmax) Vk=vmax; //boundary for v € [-1,1]:
            TFloat i = 1.0-Qi*Vi;
            TFloat j = 1.0-Qj*Vj;
            TFloat k = 1.0-Qk*Vk;
            C = fmin(i, fmin(j, k));
            C = C / 2.0 ;
            if (C<0.0) C=0.0;
            if (C>1.0) C=1.0;
            //voltages:
            TFloat Gi = Qi * fmin(j,k) / 2.0;
            TFloat Gj = Qj * fmin(i,k) / 2.0;
            TFloat Gk = Qk * fmin(i,j) / 2.0;
            TFloat Ri, Rj, Rk;
            if (C != i/2.0 ) {Ri = 0.0;} else {Ri = (Qi - Vi) / 2.0;}
            if (C != j/2.0 ) {Rj = 0.0;} else {Rj = (Qj - Vj) / 2.0;}
            if (C != k/2.0 ) {Rk = 0.0;} else {Rk = (Qk - Vk) / 2.0;}
            atomicAddd(&d_dxdt[liti-1], (Xl * Xs * Gi + (1.0 + m_zeta * Xl) * (1.0 - Xs) * Ri) );
            atomicAddd(&d_dxdt[litj-1], (Xl * Xs * Gj + (1.0 + m_zeta * Xl) * (1.0 - Xs) * Rj) );
            atomicAddd(&d_dxdt[litk-1], (Xl * Xs * Gk + (1.0 + m_zeta * Xl) * (1.0 - Xs) * Rk) );

        }
        if (k==2) {
            int liti = abs(d_cls[clause*max_lits_system+0]);
            int litj = abs(d_cls[clause*max_lits_system+1]);
            TFloat Vi = d_x[liti-1]; if (Vi<vmin) Vi=vmin; if (Vi>vmax) Vi=vmax; //boundary for v € [-1,1]:
            TFloat Vj = d_x[litj-1]; if (Vj<vmin) Vj=vmin; if (Vj>vmax) Vj=vmax; //boundary for v € [-1,1]:
            TFloat i = 1.0-Qi*Vi;
            TFloat j = 1.0-Qj*Vj;
            C = fmin(i, j);
            C = C / 2.0;
            if (C<0.0) C=0.0;
            if (C>1.0) C=1.0;
            //voltages:
            TFloat Gi = Qi * j / 2.0;
            TFloat Gj = Qj * i / 2.0;
            TFloat Ri, Rj;
            if (C != i/ 2.0 ) {Ri = 0.0;} else {Ri = (Qi - Vi) / 2.0;}
            if (C != j/ 2.0 ) {Rj = 0.0;} else {Rj = (Qj - Vj) / 2.0;}
            atomicAddd(&d_dxdt[liti-1], (Xl * Xs * Gi + (1.0 + m_zeta * Xl) * (1.0 - Xs) * Ri) );
            atomicAddd(&d_dxdt[litj-1], (Xl * Xs * Gj + (1.0 + m_zeta * Xl) * (1.0 - Xs) * Rj) );
        }
        if (k!=3 && k!=2) {
            int lit[max_lits_system];
            TFloat Q[max_lits_system], V[max_lits_system], _i[max_lits_system], R[max_lits_system], G[max_lits_system];
            TFloat c_min=INT_MAX;
            for (int i=0; i<k; i++) {
                Q[i] = (d_cls[clause*max_lits_system+i]>0)? 1.0:-1.0; // +1 if literal is >0, otherwise -1
                lit[i] = abs(d_cls[clause*max_lits_system+i]);
                V[i] = d_x[lit[i]-1]; if (V[i]<vmin) V[i]=vmin; if (V[i]>vmax) V[i]=vmax; //boundary for v € [-1,1]:
                _i[i] = 1.0-Q[i]*V[i];
                // find min:
                if (_i[i]<c_min) c_min = _i[i]; 
            }
            C = c_min / 2.0;
            //voltages:            
            for (int i=0; i<k; i++) {
                //find min of others:
                TFloat g_min = INT_MAX;
                for (int x=0; x<k; x++) {if (x!=i && _i[x]<g_min) g_min = _i[x];}
                G[i] = Q[i] * g_min / 2.0;
                TFloat comp = (1.0-Q[i]*V[i])/2.0;
                if (C != comp) {R[i] = 0.0;} else {R[i] = (Q[i] - V[i]) / 2.0;}
                atomicAddd(&d_dxdt[lit[i]-1], (Xl * Xs * G[i] + (1.0 + m_zeta * Xl) * (1.0 - Xs) * R[i]) );
            }    
        }

        //update #satsified? 
        if (C<0.5) atomicAdd(&d_loc[0], -1); //this clause is sat
        //update energy:
        atomicAdd(&d_energy[0], C);
        // Calculate Xs:
        d_dxdt[n+clause] = m_beta * (Xs + m_epsilon) * (C - m_gamma);
        
        // Calculate Xl:
        d_dxdt[n+m+clause] = m_alpha * (C - m_delta);

    }
    
}

//-------------------------------------------------------------------------------------------------------------------
// this functions solves the CNF with GPU
int solveGPU() {

    printf("c [GPU] STARTING INTEGRATION AT t=%.5f\n",t_rem);
    int threadsPerBlock    = 256;
    int blocksPerGrid      =  ( (n+m*2) + threadsPerBlock - 1 ) / threadsPerBlock;// + 1;
    int threadsPerBlock_m  = 256;
    int blocksPerGrid_m    = ( m + threadsPerBlock_m - 1 ) / threadsPerBlock_m;// + 1;
    int threadsPerBlock_n  = 256;
    int blocksPerGrid_n    = ( n + threadsPerBlock_n - 1 ) / threadsPerBlock_n;// + 1;

    // initiate d_x;
    gpuErrorCheck(cudaMemcpy(d_x, initial_assignments, (n+m*2)*sizeof(TFloat), cudaMemcpyHostToDevice),true); //only if seed is not 0; seed 0 sets everything to zero
    printf("c [GPU] INITIAL ASSIGNMENTS: %.5f, %.5f, %.5f, %.5f, %.5f\n",initial_assignments[0], initial_assignments[1], initial_assignments[2], initial_assignments[3], initial_assignments[4]);
     
    int integration_steps;
    int no_improvement_since;
    TFloat t;
    double t_begin, t_end;
    double global_energy = m;

    printf("c [GPU] STARTING INTEGRATION...\n");    
    printf("c [GPU] PARAMETERS: α=%.15f β=%.15f γ=%.15f ε=%.15f δ=%.15f ζ=%.15f\n",dmm_alpha, dmm_beta, dmm_gamma, dmm_delta, dmm_epsilon, dmm_zeta);
    printf("c [GPU] INTEGRATION MODE: %d\n",INTEGRATION_MODE);
    printf("c [GPU] STARTING STEPSIZE=%.15f\n",stepsize);

    solved = 0;
    global = m;
    t = t_rem;

    integration_steps = 0;
    no_improvement_since = 0;

    t_begin = clock(); 
    t_abs_begin = t_begin;

    int corr_best = 0;
    
    gpu_reset_vars<<<1,1>>>(d_loc, d_energy, d_varchanges, m); // init loc = m:
    
    /// main integration routine: ///
    while (solved == 0) {
        
        // do one step (updates d_dxdt, d_loc[0]):
        gpu_reset_dxdt<<<blocksPerGrid,threadsPerBlock>>>(d_dxdt, n+m*2); // kernel reset dxdf = all zero:
        gpu_reset_vars<<<1,1>>>(d_loc, d_energy, d_varchanges, m); // init loc = m:
        gpu_step<<<blocksPerGrid_m,threadsPerBlock_m>>>(d_x, d_dxdt, d_cls, d_loc, d_energy, m, n, xl_max, dmm_alpha, dmm_beta, dmm_gamma, dmm_delta, dmm_epsilon, dmm_zeta, d_clauseSizes);
        
        // only periodically update screen etc:
        if (integration_steps % 1 == 0)  {
		
		// update time:
		t_end = clock(); 
		double time_spent = (double)(t_end - t_begin)/CLOCKS_PER_SEC;//1000;// / 
		t_abs_end = clock();
		double time_abs_spent = (double)(t_abs_end - t_abs_begin)/CLOCKS_PER_SEC;//1000;// / 

		//get loc:
		cudaMemcpy(h_loc, d_loc, sizeof(int), cudaMemcpyDeviceToHost);

		//get energy:
		cudaMemcpy(h_energy, d_energy, sizeof(float), cudaMemcpyDeviceToHost);
		
		//screen output:
		if (tune) {
		    printf("\rc [GPU] %.2fs \tt=%.5f \tglobal=%6d (%6d) [T:%6d] steps=%9d (+%.5f) α=%.5f β=%.5f γ=%.5f ε=%.5f δ=%.5f ζ=%.5f E=%.5f ", time_spent,t,global,h_loc[0],tune_global,integration_steps,stepsize,dmm_alpha, dmm_beta, dmm_gamma, dmm_delta, dmm_epsilon, dmm_zeta, h_energy[0]);    
		} else {
		   printf("\rc [GPU] %.2fs \tt=%.5f \tglobal=%6d (%6d) steps=%9d (+%.5f) α=%.5f β=%.5f γ=%.5f ε=%.5f δ=%.5f ζ=%.5f E=%.5f (%.5f) ", time_spent,t,global,h_loc[0],integration_steps,stepsize,dmm_alpha, dmm_beta, dmm_gamma, dmm_delta, dmm_epsilon, dmm_zeta, h_energy[0], global_energy);    
		}
		
		// solved?
		if (h_loc[0]==0) {
		    solved = 1;
		    break;
		}
		// better loc?
		if (h_loc[0]<global) {
		    global = h_loc[0];
		    if (!tune) printf("\n");
		    if (global<2000 & (!tune || (tune && global<tune_global))) {
		        //write solution to file:
		        cudaMemcpy(h_v, d_x, n*sizeof(TFloat), cudaMemcpyDeviceToHost);
		        FILE *fs = fopen(SOLUTION_FILE, "w");
		        for (int i=0; i<n; i++) fprintf(fs,"%.32f, ",h_v[i]); // current G_field = solution
		        fclose(fs);       
		    }
		    if (tune && tune_mode==1) {
		        cudaMemcpy(initial_assignments, d_x, n*sizeof(TFloat), cudaMemcpyDeviceToHost);
		    }                
		    // restart counters reset after better loc found:
		    no_improvement_since = 0;
		    if (!tune) t_begin = clock(); 
		} else {
			// better energy?
			if (h_energy[0] < global_energy) {
				global_energy = h_energy[0];
				if (!tune) printf("\n");
			}
	    	}

		// exit on: walltime_timeout, maxsteps (no improvement since)
		if (time_spent>walltime_timeout || no_improvement_since>maxsteps || (tune && time_abs_spent>walltime_abs_timeout)){
		    printf(TEXT_SILVER);
		    if (no_improvement_since>maxsteps) printf("\nc [GPU] BREAK ON NO IMPROVEMENTS SINCE LIMIT %d\n",no_improvement_since);
		    if (time_spent>walltime_timeout) printf("\nc [GPU] BREAK ON WALLTIME LIMIT %.0fs\n",walltime_timeout);
		    if (time_abs_spent>walltime_abs_timeout) printf("\nc [GPU] BREAK ON WALLTIME ABS LIMIT %.0fs\n",walltime_abs_timeout);
		    printf(TEXT_DEFAULT);
		    std::cout << "c integration steps=" << integration_steps << std::endl;
		    return global;
		}
	}

        /// apply rhs: ------------------------------------------------------------------------------------------------------------------------------
        t = t + stepsize;

        /// CUSTOM FORWARD ADAPTIVE EULER: ////////////////////////////////////////
        bool adaptive_accepted = false;
        TFloat h_min = 0.0078125;
        TFloat h_max = 1.0;
        stepsize = 0.125;
        while (!adaptive_accepted) {
	        gpu_reset_vars<<<1,1>>>(d_loc, d_energy, d_varchanges, m); // reset d_varchanges -> WE SHOULD HAVE ITS OWN ROUTINE HERE!
        	gpu_euler<<<blocksPerGrid,threadsPerBlock>>>(d_x_tmp, d_x, d_dxdt, n+m*2, t, stepsize, n, m , xl_max, d_varchanges); // MOVE ENTIRELY TO CUDA
        	//get #changed vars:
		cudaMemcpy(h_varchanges, d_varchanges, sizeof(int), cudaMemcpyDeviceToHost);
		if (h_varchanges[0]==0) adaptive_accepted=true;
		stepsize = stepsize * 1/2;
		if (stepsize<=h_min) {
			stepsize = h_min;
			adaptive_accepted=true;
		}
        }
        /// --- CUSTOM FORWARD //////////////////////////////////////////////////////

        // do the final step:
        gpu_euler<<<blocksPerGrid,threadsPerBlock>>>(d_x, d_x, d_dxdt, n+m*2, t, stepsize, n, m , xl_max, d_varchanges); //we can optimize -> reuse from adaptive
        integration_steps++;
        no_improvement_since++;
    	}
    /// --- 
    std::cout << "c integration steps=" << integration_steps << std::endl;
    
    // output solution:
    cudaMemcpy(h_v, d_x, n*sizeof(TFloat), cudaMemcpyDeviceToHost);

    printf(TEXT_YELLOW); printf("v [GPU] ");
    for (int i=0; i<n; i++) {
        if (h_v[i]>0) printf("%d ",i+1);
        if (h_v[i]<0) printf("%d ",(i+1)*-1);
        if (h_v[i]==0) printf("%d ",(i+1));
    }
    printf("\n"); printf(TEXT_DEFAULT);
    
    // verify solution:
    printf("\nc [GPU] VERIFYING...\n"); printf(TEXT_DEFAULT);
    bool sat = true; bool clausesat;
    for (int i=0; i<m; i++) {
        for (int j=0; j<clauseSizes[i]; j++) {
            clausesat = false;
            int lit = abs(cls[i*max_lits_system+j]);
            if ( (h_v[lit-1]>0 && cls[i*max_lits_system+j]>0) || (h_v[lit-1]<0 && cls[i*max_lits_system+j]<0) || (h_v[lit-1]==0 && cls[i*max_lits_system+j]>0) ) {
                clausesat = true;
                break;
            }
        }
        if (!clausesat) {
            sat = false;
            //output wrong assignment:
            printf("CLAUSE %d [",i);
            for (int j=0; j<clauseSizes[i]; j++) printf("%d ",cls[i*max_lits_system+j]);
            printf("] -> ");
        	for (int j=0; j<clauseSizes[i]; j++) printf("%.2f ",h_v[abs(cls[i*max_lits_system+j])]);
        	printf(" IS NOT SAT.\n");
            
            break;
        }
    }
    if (sat)  {
        printf(TEXT_YELLOW); printf("c [GPU] SAT (VERIFIED)\n"); solved = 1;
        //write solution to file:
        FILE *fs = fopen(SOLUTION_FILE, "w");
        for (int i=0; i<n; i++) fprintf(fs,"%.15f, ",h_v[i]); // current G_field = solution
        fclose(fs);
    }
    if (!sat) {printf(TEXT_RED); printf("c [GPU] UNSAT (VERIFIED)\n");}
    printf(TEXT_DEFAULT);

    return 0;
}
//-------------------------------------------------------------------------------------------------------------------
// run on GPU and tune on the go:
void tuneGPU() {
    tune_global = m;
    int tune_loc = m;
    int PC = 7;
    TFloat rem_tune_param[PC];
    int    rem_tune_param_idx[PC];
    int tune_iteration = 0;
    int tune_no_improvement_since = 0;
    bool tune_change_params = false;

    /* RNG std::uniform_real_distribution<double> */
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int>      rand_param(0,PC-1); //which params to change TEST: DONT CHANGE STEPSIZE - BREAKS THE FLOW FIELD?
    std::lognormal_distribution<double>     rand_num_vars(0.0,0.5);
    std::uniform_real_distribution<double>      rand_alpha(1e-16, 1e2);
    std::uniform_real_distribution<double>      rand_beta(1e-16, 1e2);
    std::uniform_real_distribution<double>      rand_gamma(1e-16, 0.5);
    std::uniform_real_distribution<double>      rand_delta(1e-16, 1.0);
    std::uniform_real_distribution<double>      rand_epsilon(1e-16, 0.5);
    std::uniform_real_distribution<double>      rand_zeta(1e-16, 1.0);
    std::uniform_real_distribution<double>      rand_stepsize(1e-8, 0.15);
    
    std::uniform_int_distribution<int>  rand_v(-1.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xs(0.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xl(1.0, 10.0); //xl_max); //10.00

    std::uniform_int_distribution<int>      rand_choosevar(0,n); //+m*2); //variable choosen from switchfraction
    std::uniform_int_distribution<int>      rand_choosstrategy(0,1); //tuning strategy
    std::uniform_real_distribution<double>  rand_smallstep(0.985111111, 1.015111111);

    while (solved==0) {
        printf("c ------------------------------------------------------------------------------------\n");
        printf("c ITERATION %d (NO IMPROVEMENT SINCE %d)\n",tune_iteration, tune_no_improvement_since);

        // move params => rem_tune_param:
        rem_tune_param[0] = dmm_alpha;
        rem_tune_param[1] = dmm_beta;
        rem_tune_param[2] = dmm_gamma;
        rem_tune_param[3] = dmm_delta;
        rem_tune_param[4] = dmm_epsilon;
        rem_tune_param[5] = dmm_zeta;
        rem_tune_param[6] = stepsize;
        //empty rem_tune_param_idx:
        for (int i=0; i<PC; i++) rem_tune_param_idx[i] = 0;

        if (tune_change_params) {
            // how many parameters do we change?
            double _num_change_vars = rand_num_vars(generator);
            int num_change_vars = (int)_num_change_vars; 
            if (num_change_vars<1) num_change_vars = 1;
            if (num_change_vars>PC) num_change_vars = PC;

            for (int p=0; p<num_change_vars; p++) {
                int tune_param = rand_param(generator);
                // mark the tuned param:
                rem_tune_param_idx[tune_param] = 1;
                // choose tuning strategy:
                int strategy = rand_choosstrategy(generator);
                if (strategy>0) {
                    // small move:
                    double smallstep = rand_smallstep(generator);
                    switch (tune_param) {
                        case 0: dmm_alpha =     dmm_alpha*smallstep; break;
                        case 1: dmm_beta =      dmm_beta*smallstep; break;
                        case 2: dmm_gamma =     dmm_gamma*smallstep; break;
                        case 3: dmm_delta =     dmm_delta*smallstep; break;
                        case 4: dmm_epsilon =   dmm_epsilon*smallstep; if (dmm_epsilon>=dmm_gamma) dmm_epsilon = dmm_gamma - 1e-8; break;
                        case 5: dmm_zeta =      dmm_zeta*smallstep; break;
                        case 6: stepsize =      stepsize*smallstep; break;
                    }    
                }
                if (strategy==0) {
                    // random point
                    switch (tune_param) {
                        case 0: dmm_alpha =     rand_alpha(generator); break;
                        case 1: dmm_beta =      rand_beta(generator); break;
                        case 2: dmm_gamma =     rand_gamma(generator); break;
                        case 3: dmm_delta =     rand_delta(generator); break;
                        case 4: dmm_epsilon =   rand_epsilon(generator); if (dmm_epsilon>=dmm_gamma) dmm_epsilon = dmm_gamma - 1e-8; break;
                        case 5: dmm_zeta =      rand_zeta(generator); break;
                        case 6: stepsize =      rand_stepsize(generator); break;
                    }    
                }
                
                printf(TEXT_SILVER);
                switch (tune_param) {
                    case 0: printf("c STRATEGY %d CHANGED α TO %.15f \n", strategy, dmm_alpha); break;
                    case 1: printf("c STRATEGY %d CHANGED β TO %.15f \n", strategy, dmm_beta); break;
                    case 2: printf("c STRATEGY %d CHANGED γ TO %.15f \n", strategy, dmm_gamma); break;
                    case 3: printf("c STRATEGY %d CHANGED ε TO %.15f \n", strategy, dmm_delta); break;
                    case 4: printf("c STRATEGY %d CHANGED δ TO %.15f \n", strategy, dmm_epsilon); break;
                    case 5: printf("c STRATEGY %d CHANGED ζ TO %.15f \n", strategy, dmm_zeta); break;
                    case 6: printf("c STRATEGY %d CHANGED stepsize TO %.15f \n", strategy, stepsize); break;
                }
                printf(TEXT_DEFAULT);
            }    
        }
        
        // run GPU integration:
        tune_loc = solveGPU();
        
        tune_change_params = true;

        // did we improve? (we also accept similar loc as this deviates the parameter space a bit)
        if (tune_loc<tune_global) {
            if (tune_loc<tune_global) tune_no_improvement_since = 0;
            tune_global = tune_loc;
            printf(TEXT_YELLOW);
            printf("c LOC=%d ACCEPTED.\n",tune_loc);
            printf(TEXT_DEFAULT);
            //tune_change_params = false;
        } else {
            printf(TEXT_SILVER);
            printf("c LOC=%d (GLOBAL=%d) - REJECTED\n",tune_loc, tune_global);
            printf(TEXT_DEFAULT);
            // reset parameters:
            for (int p=0; p<PC; p++) {
                if (rem_tune_param_idx[p]==1) {
                    switch (p) {
                        case 0: dmm_alpha = rem_tune_param[0]; break;
                        case 1: dmm_beta = rem_tune_param[1]; break;
                        case 2: dmm_gamma = rem_tune_param[2]; break;
                        case 3: dmm_delta = rem_tune_param[3]; break;
                        case 4: dmm_epsilon = rem_tune_param[4]; break;
                        case 5: dmm_zeta = rem_tune_param[5]; break;
                        case 6: stepsize = rem_tune_param[6]; break;
                    }
                }
            }
        }
        
        // adjust walltime_timeout? only if not continous...
        if (tune_no_improvement_since>=50) {//25) {
            tune_no_improvement_since = 0;
            walltime_timeout += 1.0;
            walltime_abs_timeout += 1.0;
            tune_change_params = false; // we increased the time, so we want to see if it gets better with the current (best) params...
            printf(TEXT_CYAN);
            printf("c UPDATED walltime_timout TO %.2fs, walltime_abs_timeout TO %.2fs\n",walltime_timeout,walltime_abs_timeout);
            printf(TEXT_DEFAULT);
        }

        tune_iteration++;
        tune_no_improvement_since++;

        // switchfraction:
        double _switchvars = n * switchfraction;
        int switchvars = abs(_switchvars);
        for (int i=0; i<switchvars; i++) {
            int _var = rand_choosevar(generator);
            if (_var<n) initial_assignments[_var] = rand_v(generator);
            if (_var>=n && _var < n+m) initial_assignments[_var] = rand_Xs(generator);
            if (_var>=n+m) initial_assignments[_var] = rand_Xl(generator);
        }
        printf("c SWITCHFRACTION CHANGED %d VARS\n",switchvars);

    }
    printf("c FINISHED. SOLUTION FOUND.\n");
    
}

//-------------------------------------------------------------------------------------------------------------------
// keyboard runtime menu
void  INThandler(int sig)
{
    signal(SIGINT, INThandler);

}


//-------------------------------------------------------------------------------------------------------------------
//parse command line options
int scan_opt(int argc, char **argv, const char *opt) {
    char c;
    while ((c = getopt (argc, argv, opt)) != -1)
        switch (c) {
            case 't': tune=atoi(optarg); break;
            case 'w': THREAD_COUNT=atoi(optarg); break;
            case 'q': quiet=atoi(optarg); break;
            case 'o': INTEGRATION_MODE=atoi(optarg); break;
            case 'i': strcpy(input_filename, optarg); break;
            case 'd': init_dt=atof(optarg); break;
            case 'x': rk_errorrate_1=atof(optarg); break;
            case 'y': rk_errorrate_2=atof(optarg); break;
            case 'z': maxsteps=atof(optarg); break;
            case 's': stepsize=atof(optarg); break;
            case 'a': timeout   =atof(optarg); break;
            
            case 'e': load_partable=atoi(optarg); break;
            case 'f': adaptive_max   =atof(optarg); break;
            case 'g': tune_mode_params   =atoi(optarg); break;

            case 'b': dmm_beta=atof(optarg); break;
            case 'c': dmm_alpha=atof(optarg); break;
            case 'n': dmm_gamma=atof(optarg); break;
            case 'h': dmm_delta=atof(optarg); break;
            case 'j': dmm_epsilon=atof(optarg); break;
            case 'k': dmm_zeta=atof(optarg); break;
            case 'l': seed=atoi(optarg); break;
            case 'm': xl_max=atoi(optarg); break;

            case 'p': heuristics=atoi(optarg); break;
            case 'u': walltime_timeout   =atof(optarg); break;

            case 'v': switchfraction   =atof(optarg); break;
            case 'r': tune_mode   =atoi(optarg); break;
            
            default: return(-1);
        }
    return(0);
}

//-------------------------------------------------------------------------------------------------------------------
void printDevProp(cudaDeviceProp devProp)
{
    printf(TEXT_GREEN);
    printf("c [GPU] Major revision number:         %d\n",  devProp.major);
    printf("c [GPU] Minor revision number:         %d\n",  devProp.minor);
    printf("c [GPU] Name:                          %s\n",  devProp.name);
    printf("c [GPU] Total global memory:           %u\n",  devProp.totalGlobalMem);
    printf("c [GPU] Total shared memory per block: %u\n",  devProp.sharedMemPerBlock);
    printf("c [GPU] Total registers per block:     %d\n",  devProp.regsPerBlock);
    printf("c [GPU] Warp size:                     %d\n",  devProp.warpSize);
    printf("c [GPU] Maximum memory pitch:          %u\n",  devProp.memPitch);
    printf("c [GPU] Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
    for (int i = 0; i < 3; ++i)
    printf("c [GPU] Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
    for (int i = 0; i < 3; ++i)
    printf("c [GPU] Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
    printf("c [GPU] Clock rate:                    %d\n",  devProp.clockRate);
    printf("c [GPU] Total constant memory:         %u\n",  devProp.totalConstMem);
    printf("c [GPU] Texture alignment:             %u\n",  devProp.textureAlignment);
    printf("c [GPU] Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
    printf("c [GPU] Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
    printf("c [GPU] Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
    printf(TEXT_DEFAULT);
    return;
}

//-------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
    int i, j;
    char buffer[32];
    global_best_thread = -1;

    char *syntax =
    "c          GENERAL OPTIONS:\n"
    "c          -i file      : INPUT FILE (CNF)\n"
    "c          -o [1;2;3;4] : ODE INTEGRATION: 1=CONSTANT, 2=CUSTOM ADAPTIVE, 3=RUNGE KUTTA, 4=IMPLICIT (DEFAULT:3)\n"
    "c          -q [0;1]     : QUIET MODE; 0=OFF 1=ON (DEFAULT:0)\n"
    "c          -u [double]  : WALL TIME TIMEOUT AFTER x s\n"
    "c          -w [int]     : NUMBER OF PARALLEL THREADS\n"
    "c          -e [0;1]     : LOAD partable.txt\n"
    "c\n"
    "c          -c [double]  : ALPHA (GROWTH RATE FOR LONG TERM MEMORY Xl)\n"
    "c          -b [double]  : BETA (GROWTH RATE FOR SHORT TERM MEMORY Xs)\n"
    "c          -n [double]  : GAMMA (RESTRICTION FOR CM IN Xs)\n"
    "c          -h [double]  : DELTA (GROWTH RATE FOR CM IN Xl)\n"
    "c          -j [double]  : EPSILON (REMOVE SPURIOUS SOLUTION Xm=0)\n"
    "c          -k [double]  : ZETA (REDUCTION FACTOR OF RIGIDITY G)\n"
    "c          -l [int]     : RANDOM SEED\n"
    "c          -m [int]     : MAX VALUE FOR Xl\n"
    "c\n"
    "c          -p [0;1]     : ALPHA HEURISTICS; 0=OFF 1=ON (DEFAULT:1)\n"
    "c\n"
    "c          TUNING OPTIONS:\n"
    "c          -t [0;1]     : TUNE CIRCUIT; 0=OFF 1=ON (DEFAULT:0)\n"
    "c          -g [0;1;2]   : PARAMS TO TUNE: 0:α..δ 1:ODE 2:ALL\n"
    "c          -r [0;1]     : TUNING MODE: 0=AWAYS FROM START; 1=CONTINOUS FROM BEST; 2=CONTINOUS OVERALL BEST\n"
    "c          -v [double]  : SWITCH FRACTION (DEFAULT 0.0001)\n"
    "c\n"
    "c          ODE - CONSTANT OPTIONS:\n"
    "c          -s [double]  : STEP SIZE\n"
    "c          -a [double]  : TIME OUT (INTEGRATION TIME)\n"
    "c\n"
    "c          ODE - RUNGE KUTTA OPTIONS:\n"
    "c          -x [double]  : ERROR RATE 1\n"
    "c          -y [double]  : ERROR RATE 2\n"
    "c          -d [double]  : INITIAL dt VALUE\n"
    "c          -z [double]  : MAXIMUM INTEGRATION STEPS\n"
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

    if (!quiet) printf("c ------------------------------------------------------------------------------------\n");
    if (!quiet) printf("c TURINGX SAT-SOLVER GPU VERSION                                               (C)2022\n");
    if (!quiet) printf("c ------------------------------------------------------------------------------------\n");
    if (!quiet) printf("c INSTANCE  : %s\n", input_filename);

    /// load CNF header:
    FILE *file = fopen(input_filename, "r");
    if (strcmp(buffer, "c") == 0) {
        while (strcmp(buffer, "\n") != 0) {
            fscanf(file, "%s", buffer);
        }
    }
    while (strcmp(buffer, "p") != 0) {
        fscanf(file, "%s", buffer);
    }
    fscanf(file, " cnf %i %i", &n, &m);

    if (coupledsystem) m = (int)(m * 2); 
    
    if (!quiet) printf("c VARIABLES : %'d\n", n);
    if (!quiet) printf("c CLAUSES   : %'d\n", m);
    if (!quiet) printf("c RATIO     : %lf\n", (double) m / n);

    xl_max = xl_max * m;
    if (xl_max<=0) xl_max = INT_MAX;

    /// reserve  memory - needs to be done before anything else:
    cls = (int *) calloc((size_t) m*max_lits_system, sizeof(int));
    for (int i=0; i<m*max_lits_system; i++) cls[i] = 0; 
    numOccurrenceT = (int *) calloc((size_t) n+1, sizeof(int));
    clauseSizes = (int *) calloc((size_t) m, sizeof(int));
        
    /// read CNF: /////////////////////////////////////////
    int lit; int lit_coupled;
    for (i = 0; i < m; i++) {
        j = 0; 
        do {
            fscanf(file, "%s", buffer);
            if (strcmp(buffer, "c") == 0) {
                continue;
            }
            lit = atoi(buffer);
            if (lit!=0) cls[i*max_lits_system+j] = lit;

            if (coupledsystem) {
                if (lit!=0) cls[((int)(m/2)+i)*max_lits_system+j] = lit;    
            }

            // increase number of Occurence of the variable, max number of occurences
            if (lit!=0) {
                numOccurrenceT[abs(lit)]++;
                if (numOccurrenceT[abs(lit)]>maxNumOccurences) {maxNumOccurences=numOccurrenceT[abs(lit)];}
                clauseSizes[i] = j+1;
            }

            if (coupledsystem) {
                if (lit!=0) {
                    numOccurrenceT[abs(lit)]++;
                    clauseSizes[(int)(m/2)+i] = j+1;
                }    
            }

            j++;
        } while (strcmp(buffer, "0") != 0);
        j--;
        if (j > max_lits_system) {
            printf("c ERROR: CLAUSE %d HAS MORE THAN %d LITERALS.\n",i,max_lits_system);
            return EXIT_FAILURE;
        }
    }
    
    if (!quiet) printf("c MAX VARIABLE OCCURENCE: %'d\n", maxNumOccurences);

    if (!quiet) printf("c FIRST 10 CLAUSES:\n");
    for (i = 0; i < 11; i++) {
        if (!quiet) printf("c CLAUSE %i: ",i);
        for (j = 0; j < clauseSizes[i]; j++) {if (!quiet) printf(" %d",cls[i*max_lits_system+j]);}
        if (!quiet) printf(" (%d)",clauseSizes[i]);
        if (!quiet) printf("\n");
    }

    //build occurence array: [var][cls...] 
    occurrence = (int *) calloc((size_t) (n+1)*maxNumOccurences, sizeof(int));
    occurenceCounter = (int *) calloc((size_t) n+1, sizeof(int));
    
    for (i=0; i<m; i++) {
        for (j = 0; j < clauseSizes[i]; j++) {
            lit = abs(cls[i*max_lits_system+j]);
            occurrence[lit*maxNumOccurences+occurenceCounter[lit]] = i;
            occurenceCounter[lit]++;
        }
    }

    /// initialize arrays:
    v_best = (TFloat *) calloc((size_t) THREAD_COUNT*(n+m*2), sizeof(TFloat));
    loc_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    global_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    global_all_runs_thread = (int *) calloc((size_t) THREAD_COUNT, sizeof(int));
    time_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    time_thread_actual = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    walltime_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    t_begin_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    t_end_thread = (double *) calloc((size_t) THREAD_COUNT, sizeof(double));
    initial_assignments = (TFloat *) calloc((size_t) THREAD_COUNT*(n+m*2), sizeof(TFloat));
    thread_params = (double *) calloc((size_t) THREAD_COUNT*PARAM_COUNT, sizeof(double));
    partable = (double *) calloc((size_t) PARAM_COUNT*4, sizeof(double));
    defaults = (double *) calloc((size_t) PARAM_COUNT*128, sizeof(double));

    h_v = (TFloat *) calloc((size_t) n, sizeof(TFloat));

    gpuErrorCheck(cudaMalloc(&d_x, (n+m*2)*sizeof(TFloat) ),true);
    gpuErrorCheck(cudaMalloc(&d_x_tmp, (n+m*2)*sizeof(TFloat) ),true);
    gpuErrorCheck(cudaMalloc(&d_dxdt, (n+m*2)*sizeof(TFloat) ),true); 
    gpuErrorCheck(cudaMalloc(&d_cls, m*max_lits_system*sizeof(int) ),true);
    gpuErrorCheck(cudaMemcpy(d_cls, cls, m*max_lits_system*sizeof(int), cudaMemcpyHostToDevice),true);
    gpuErrorCheck(cudaMalloc(&d_clauseSizes, m*sizeof(int) ),true);
    gpuErrorCheck(cudaMemcpy(d_clauseSizes, clauseSizes, m*sizeof(int), cudaMemcpyHostToDevice),true);
    cudaMalloc(&d_loc, sizeof(int));
    cudaMalloc(&d_varchanges, sizeof(int));
    cudaMalloc(&d_energy, sizeof(float));

    if (!quiet) printf("c ARRAYS INITIALISED\n");

    /// set _all_runs vars: ---------------------------------------------------------------------
    global_all_runs = m;
    for (int i=0; i<THREAD_COUNT; i++) global_all_runs_thread[i] = m;

    /// load solution: --------------------------------------------------------------------------
    if (loadsolution) {
        FILE *fs;
        fs = fopen(SOLUTION_FILE, "r");
        if (fs == NULL) {
            fprintf(stdout, "c solution.txt NOT PROVIDED. USING RANDOM ASSIGNMENTS.\n");
        } else {
            double ivar;
            for (int i=0; i<n;i++) {
                fscanf(fs,"%lf,",&ivar);
                v_best[i] = ivar;
                //if (i<5) printf("SOLUTION_FILE p%d=%.5f, ivar=%.5f ",i,v_best[i],ivar);
            }
            fclose(fs);
            printf("c SOLUTION solution.txt LOADED.\n");
        }
    } else {
        printf("c SOLUTION solution.txt NOT LOADED (SEED!=0)\n");
    }

    /// load partable: --------------------------------------------------------------------------
    if (load_partable) {
        FILE *fp;
        fp = fopen(PARTABLE_FILE, "r");
        if (fp == NULL) {
            fprintf(stdout, "c partable.txt NOT PROVIDED. USING DEFAULT SETTINGS.\n");
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
            printf("c PARTABLE LOADED %d PARAMS AND %d DEFAULTS.\n",PARAM_COUNT,num_defaults);
            for (int i=0; i<num_defaults; i++) {
                printf("c DEFAULTS SET %d: ",i);
                for (int j=0; j<PARAM_COUNT; j++) printf("%.32f ",defaults[i*PARAM_COUNT+j]);
                printf("\n");
            }
            partable_loaded = true;
        }
    } else {
        fprintf(stdout, "c partable.txt NOT LOADED.\n");
    }

    /// OUPUT SETTINGS: ---------------------------------------------------------------------------
    if (!quiet) printf(TEXT_CYAN);
    if (!quiet) printf("c SETTINGS:\n");
    if (!quiet) printf("c #THREADS        : %d\n",THREAD_COUNT);
    if (!quiet) printf("c SWITCHFRACTION  : %.17f\n",switchfraction);
    if (INTEGRATION_MODE==ODE_RUNGEKUTTA) {
        if (!quiet) printf("c INTGRATION MODE : ADAPTIVE RUNGE KUTTA\n");
        if (!quiet) printf("c ERROR RATE 1    : %.17f\n",rk_errorrate_1);
        if (!quiet) printf("c ERROR RATE 2    : %.17f\n",rk_errorrate_2);
        if (!quiet) printf("c INITAL dt       : %.17f\n",init_dt);
        if (!quiet) printf("c MAX STEPS       : %.0f\n",maxsteps);
    }
    if (INTEGRATION_MODE==ODE_CUSTOM_ADAPTIVE) {
        if (!quiet) printf("c INTGRATION MODE : CUSTOM ADAPTIVE\n");
        if (!quiet) printf("c MIN STEPSIZE    : %.10f\n",adaptive_min);
        if (!quiet) printf("c MAX STEPSIZE    : %.2f\n",adaptive_max);
        if (!quiet) printf("c MAX ERROR       : %.2f\n",step_error);
        if (!quiet) printf("c TIMEOUT         : %.0f\n",timeout);
    }
    if (INTEGRATION_MODE==ODE_CONSTANT) {
        if (!quiet) printf("c INTGRATION MODE : CONSTANT\n");
        if (!quiet) printf("c STEPSIZE        : %.15f\n",stepsize);
        if (!quiet) printf("c TIMEOUT         : %.0f\n",timeout);
    }
    if (!quiet) printf("c TUNE CIRCUIT    : %d\n",tune);
    if (!quiet) printf("c ALPHA HEURISTICS: %d\n",heuristics);
    if (!quiet) printf("c ALPHA           : %.17f\n",dmm_alpha);
    if (!quiet) printf("c BETA            : %.17f\n",dmm_beta);
    if (!quiet) printf("c GAMMA           : %.17f\n",dmm_gamma);
    if (!quiet) printf("c DELTA           : %.17f\n",dmm_delta);
    if (!quiet) printf("c EPSILON         : %.17f\n",dmm_epsilon);
    if (!quiet) printf("c ZETA            : %.17f\n",dmm_zeta);
    if (!quiet) printf("c SEED            : %d\n",seed);
    if (!quiet) printf("c XL_MAX          : %.d\n",xl_max);

    if (!quiet) printf(TEXT_DEFAULT);

    /// init thread specific parameters //////////////////////
    
    /* RNG */
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int>  rand_v(-1.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xs(0.0, 1.0);
    std::uniform_real_distribution<double>  rand_Xl(1.0, 10.0); //xl_max); //10.00
    
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
        thread_params[i*PARAM_COUNT+8]    = init_dt;
        
        //initial assignments:
        int set_seed = seed + (i*4096);
        for (int j=0; j<n; j++) initial_assignments[i*(n+m*2)+j] = rand_v(generator);
        for (int j=n; j<n+m; j++) initial_assignments[i*(n+m*2)+j] = 0.0; //rand_Xs(generator);
        for (int j=n+m; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = 1.0;//rand_Xl(generator);
        //special case: seed 0: all V zeros:
        if (seed==0 && i==0) for (int j=0; j<n+m*2; j++) initial_assignments[i*(n+m*2)+j] = 0.0; //v_best[i];// -1.00;
        //solution loaded?
        if (loadsolution) {
            for (int j=0; j<n; j++) {
                initial_assignments[i*(n+m*2)+j] = v_best[j];
            }
        }
    }

    /// start solver (GPU) ///////////////////////////////////////////////
    printf("c SOLVING WITH GPU....\n");
    // Number of CUDA devices
    int devCount;
    cudaGetDeviceCount(&devCount);
    printf("c [GPU] There are %d CUDA devices.\n", devCount);

    // Iterate through devices:
    for (int i = 0; i < devCount; ++i)
    {
        // Get device properties
        printf("c [GPU] CUDA Device #%d:\n", i);
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, i);
        printDevProp(devProp);
    }

    // start ODE with GPU:
    t_rem = 0.0; // in tuning mode, we continue from last lwoest loc and also at that time (t_rem)
    if (tune) {
        tuneGPU();  
    } else {
        solveGPU();
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

    // free device memory:
    cudaFree(d_x);
    cudaFree(d_x_tmp);
    cudaFree(d_dxdt);
    cudaFree(d_cls);
    cudaFree(d_clauseSizes);
    cudaFree(d_loc);
    cudaFree(d_energy);
    
    return EXIT_SUCCESS;
}


