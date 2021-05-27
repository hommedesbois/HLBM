#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <time.h>

#include <sys/time.h>
#include <sched.h>

#include "ludwig.h" // lattice parameters

// constants model D3Q19
const double w[NPOP] = {4./9., // center
		1./9., 1./9., 1./9., 1./9., // face
		1./36., 1./36., 1./36., 1./36.}; //edge


/* discrete velocity
 *
 * index :            0   1   2   3   4   5   6   7   8
                      |   |   |   |   |   |   |   |   | */
const int ex[NPOP] = {0,  0,  0,  1, -1,  1, -1,  1, -1};
const int ey[NPOP] = {0,  1, -1,  0,  0,  1,  1, -1, -1};
/*                    |   |   |   |   |   |   |   |   | */
const int finv[NPOP]={0,  2,  1,  4,  3,  8,  7,  6,  5};

const double cc[NPOP]={0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5}; 

const double H2xx[NPOP] = {-1.0/3.0, -1.0/3.0, -1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0}; 
const double H2xy[NPOP] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0};
const double H2yy[NPOP] = {-1.0/3.0, 2.0/3.0, 2.0/3.0, -1.0/3.0, -1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0}; 

const double H3xxy[NPOP] = {0.0, -1.0/3.0, 1.0/3.0, 0.0, 0.0, 2.0/3.0, 2.0/3.0, -2.0/3.0, -2.0/3.0};
const double H3xyy[NPOP] = {0.0, 0.0, 0.0, -1.0/3.0, 1.0/3.0, 2.0/3.0, -2.0/3.0, 2.0/3.0, -2.0/3.0};

const double H4xxyy[NPOP] = {1.0/9.0, -2.0/9.0, -2.0/9.0, -2.0/9.0, -2.0/9.0, 4.0/9.0, 4.0/9.0, 4.0/9.0, 4.0/9.0};


/*     _______________
 *    |               |
 *    |      NS       |
 *    |_______________|
 *    |               |
 *    |               |
 *    |      LBM      |
 *    |               |
 *    |_______________|
 *    |               |
 *    |      NS       |
 *    |_______________|
 */ 

/*** The following code couples the solution of the standard LBM on a D2Q9 lattice
 * to the solution of an isothermal Navier-Stokes model as descibed in the sketch above. 
 * The coulping is achieved in a multi-grid approach. Once the information is copied from NS to LBM,
 * the entire LBM solution is copied onto the NS domain.
 ***/

const double Csound = 343.2; // sound speed
double Mach; 
double U0;
double V0; // Votex is advected across the interface with M=0.2
double Es = 0.01; // vortex strength
double CFL;

double rho0 = 1.0; // reference density
double Pref; // pressure reference
double P0; // pressure parameter of solution -> p_solver


double viscosity = 0.003432; // viscosity

double pos = 0.5; // controls relative position of the vortex in y direction
double DX; // dimensional quadrant size
double Dt; // non-dimensional time step

double omega_g, tau, tau_g; // relaxation parameters

double Cs2 = 1./3.;
double invCs2 = 3.0;
double invCs4 = 9.0;
double invCs6 = 27.0;
double invCs8 = 81.0; 

int current_slot = 0, other_slot = 1; // two-field configuration

 /** Initialization
  */
void InitializeFluidNS(double **, int);
void InitializeFluidLB(double *, double **, int);

void GuardCells(double *, double *);

/**  Propagate
  */
// LBM
void Streaming(double *);
// NS
void FiniteVol(double **);
void Flux(double **, double ***, int);
void FiniteVolHeun(double **, double ***);

void FiniteVolRK2_1(double **, double ***, double ***, double ***);
void FiniteVolRK2_2(double **, double ***, double ***, double ***, double *);

void FiniteVolRK3_1(double **, double ***, double ***, double ***);
void FiniteVolRK3_2(double **, double ***, double ***, double ***);
void FiniteVolRK3_3(double **, double ***, double ***, double ***, double *);

void FiniteVolRK4_1(double **, double ***, double ***, double ***);
void FiniteVolRK4_23(double **, double ***, double ***, double ***, int);
void FiniteVolRK4_4(double **, double ***, double ***, double ***, double *);

/**  Collision
  */
void ComputeFcolFromCentMomCollision(double *, double **);
void HRRCollision(double *, double **, double **);
void RBGKCollision(double *, double **, int);

/**  Moments
  */

void ComputeMacroFromF(double *, double **, double **);
void ComputeGrad(double **, double **, double **);


/** Boundary Conditions
 */

void ComputeFcol_BC(double *);
void ComputeMacro_BC(double **);
void ComputeMacroLB_BC(double **);
void ComputeMacroP_BC(double **);
void GuardRK(double ***, int);
void GuardRKP(double ***, int);


void Closure(double **);
void ClosureRK4(double ***, int);

void ComputeMassNS(double **);
void ComputeMassLB(double *);

/** Coupling
 */  
void ComputeTransferLBM2FV(double *, double **);

/** Galileian Correction
 */
void ComputeCorrection(double **, double **);

/** Sponge Zone
  */
void ComputeSigma(double *, int);

/**  Output functions VTK
  */
void DumpMacroNS(double **, int);
void DumpMacroLB(double **, int);

/** Filter
  */
void FilterX7_NS(double **, int);
void FilterY7_NS(double **, int);

void FilterX7_LB(double **, int);
void FilterY7_LB(double **, int);


int main (int argc, const char * argv[]){
  clock_t start_t, end_t, total_t;
  int i, b, k;

  Pref = (Csound*Csound)/1.4;
  P0 = Csound * Csound;

	DX = 0.01;
	Dt = 1.;

  double *sigma = (double *) malloc (xMax_NS*yMax_NS*sizeof(double));

  double *fcol;

  double **macro_NS, ***flux, ***RK, ***QRK, ***macro_Surf_Int;
  double **macro_LB, **grad_LB;

  fcol = (double *) malloc (N_fcol * sizeof(double));

  macro_NS = (double **) malloc (2*xMaxp_NS*yMaxp_NS*sizeof(double *));
  for(i=0; i<(2 * xMaxp_NS * yMaxp_NS); i++){
    macro_NS[i] = (double *)malloc(sizeof(double)*10);
        }
  macro_LB = (double **) malloc (2*xMaxp_LB*yMaxp_LB*sizeof(double *));
  for(i=0; i<(2 * xMaxp_LB * yMaxp_LB); i++){
    macro_LB[i] = (double *)malloc(sizeof(double)*10);
        }
  grad_LB = (double **) malloc (2*xMaxp_LB*yMaxp_LB*sizeof(double *));
  for(i=0; i<(2 * xMaxp_LB * yMaxp_LB); i++){
    grad_LB[i] = (double *)malloc(sizeof(double)*6);
        }

  flux = (double ***) malloc (xMaxp_NS*yMaxp_NS*sizeof(double **));
  for(i=0; i<(xMaxp_NS * yMaxp_NS); i++){
    flux[i] = (double **)malloc(sizeof(double *)*6);
    for(k=0; k<6; k++){
      flux[i][k] = (double *)malloc(sizeof(double)*4);
    }
  }

  QRK = (double ***) malloc (4*sizeof(double **));
  for(b=0; b<4; b++){
    QRK[b] = (double **) malloc(xMaxp_NS*yMaxp_NS*sizeof(double *));
    for(i=0; i<(xMaxp_NS*yMaxp_NS); i++){
      QRK[b][i] = (double *)malloc(sizeof(double)*3);
    }
  }

  RK = (double ***) malloc (4*sizeof(double **));
  for(b=0; b<4; b++){
    RK[b] = (double **)malloc(xMaxp_NS * yMaxp_NS*sizeof(double *));
    for(i=0; i<(xMaxp_NS * yMaxp_NS); i++){
      RK[b][i] = (double *)malloc(sizeof(double)*7);
    }
  }

  macro_Surf_Int = (double ***) malloc (2*sizeof(double **));
  for(b=0; b<2; b++){
    macro_Surf_Int[b] = (double **)malloc(xMaxp_NS * sizeof(double *));
    for(i=0; i<xMaxp_NS; i++){
      macro_Surf_Int[b][i] = (double *)malloc(sizeof(double)*7);
    }
  }
    /** Simulation parameters
	 */

  enum{shearlayer=0, vortex, pulse};

  // set testcase
  int testcase = shearlayer;
  U0 = 102.96;
  V0 = 0.0;

	Mach = sqrt(U0*U0 + V0*V0) / Csound;
	CFL = (Mach + 1.)/sqrt(3.)*Dt; // CFL number with respect to the macroscopic velocity
  fprintf(stdout, "####### SIMULATION PARAMETERS ######\n");
  fprintf(stdout, "Cs = %3.2f m/s \n", Csound);
  fprintf(stdout, "U_mag = %3.2f m/s \n", sqrt(U0*U0 + V0*V0));
	fprintf(stdout, "Mach = %3.2f \n", Mach);
	fprintf(stdout, "CFL_NS = %3.2f \n", CFL);
  fprintf(stdout, "Dx = %g m \n", DX);
  fprintf(stdout, "Dt = %g s \n", (DX/(sqrt(3)*Csound)));

  fprintf(stdout, "##### END SIMULATION PARAMETERS ####\n");


  // non-dimensional relaxation time

  tau = sqrt(3)*viscosity/(Csound * DX); 
  fprintf(stdout, "tau = %g \n", tau);

  tau_g = 0.5 + sqrt(3)*viscosity/(Csound * DX);
  fprintf(stdout, "tau_g = %g \n", tau_g);

  omega_g = 1. / tau_g;
  fprintf(stdout, "omega_g = %g \n", omega_g);

  
  // Initialize simulation in fluid domain (NS) and resolution domain 1 (LB)
  InitializeFluidNS(macro_NS, testcase); // -> NS
  InitializeFluidLB(fcol, macro_LB, testcase); // -> LBM

	int s;
	int index = 0;

	system("mkdir -p vtk");
  system("rm vtk/*.*");

  double cpu = 0.0;

	for(s=0; s<NSTEPS; s++){
    start_t = clock();

    /**  PROPAGATION
      */
    Streaming(fcol); // LBM

    switch(time_marching){ // can be changed in header file (ludwig.h)

      case 0: // Four-stage Runge Kutta scheme for FV-NS

        /*
         * STEP 1
         */

        FiniteVolRK4_1(macro_NS, RK, QRK, macro_Surf_Int);
        GuardRK(RK, 1); // treament periodic boundaries rho, u
        ClosureRK4(RK, 1); // Computation of S_ij with updated velocities
        GuardRKP(RK, 1); // treament periodic boundaries P

        /*
         * STEP 2
         */

        FiniteVolRK4_23(macro_NS, RK, QRK, macro_Surf_Int, 2);
        GuardRK(RK, 2);
        ClosureRK4(RK, 2);
        GuardRKP(RK, 2);

        /*
         * STEP 3
         */

        FiniteVolRK4_23(macro_NS, RK, QRK, macro_Surf_Int, 3);
        GuardRK(RK, 3);
        ClosureRK4(RK, 3);
        GuardRKP(RK, 3);

        /*
         * STEP 4
         */

        FiniteVolRK4_4(macro_NS, RK, QRK, macro_Surf_Int, sigma);
        ComputeMacro_BC(macro_NS);
        Closure(macro_NS);
        ComputeMacroP_BC(macro_NS);
      break;

      case 1: // Three-stage Runge Kutta scheme for FV-NS

        /*
         * STEP 1
         */

        FiniteVolRK3_1(macro_NS, RK, QRK, macro_Surf_Int);
        GuardRK(RK, 1); // treament periodic boundaries rho, u
        ClosureRK4(RK, 1); // Computation of S_ij with updated velocities
        GuardRKP(RK, 1); // treament periodic boundaries P

        /*
         * STEP 2
         */

        FiniteVolRK3_2(macro_NS, RK, QRK, macro_Surf_Int);
        GuardRK(RK, 2);
        ClosureRK4(RK, 2);
        GuardRKP(RK, 2);

        /*
         * STEP 3
         */

        FiniteVolRK3_3(macro_NS, RK, QRK, macro_Surf_Int, sigma);
        ComputeMacro_BC(macro_NS);
        Closure(macro_NS);
        ComputeMacroP_BC(macro_NS);
      break;

      case 2: // Two-stage Runge Kutta scheme for FV-NS

        /*
         * STEP 1
         */

        FiniteVolRK2_1(macro_NS, RK, QRK, macro_Surf_Int);
        GuardRK(RK, 1); // treament periodic boundaries rho, u
        ClosureRK4(RK, 1); // Computation of S_ij with updated velocities
        GuardRKP(RK, 1); // treament periodic boundaries P

        /*
         * STEP 2
         */

        FiniteVolRK2_2(macro_NS, RK, QRK, macro_Surf_Int, sigma);
        ComputeMacro_BC(macro_NS);
        Closure(macro_NS);
        ComputeMacroP_BC(macro_NS);
      break;

      case 3: // Heun Predictor Corrector

        /*
         * STEP 1
         */

        FiniteVol(macro_NS);
        ComputeMacro_BC(macro_NS);
        Closure(macro_NS);
        ComputeMacroP_BC(macro_NS);

        Flux(macro_NS, flux, other_slot); // calculate fluxes at intermediate time-step

        /*
         * STEP 2
         */

        FiniteVolHeun(macro_NS, flux);
        ComputeMacro_BC(macro_NS);
        Closure(macro_NS);
        ComputeMacroP_BC(macro_NS);
      break;
    }

    /**  END PROPAGATION
      */
    other_slot = current_slot;
    current_slot = 1 - current_slot;

    /**  COLLISION
      */
    ComputeMacroFromF(fcol, macro_NS, macro_LB); // Transfer NS -> LBM taken into account here
    ComputeMacroLB_BC(macro_LB);
    ComputeGrad(macro_LB, macro_NS, grad_LB);
    //FilterX7_LB(macro_LB, current_slot); FilterY7_LB(macro_LB, current_slot); // Filter LBM solution
    //ComputeFcolFromCentMomCollision(fcol, macro_LB);
    HRRCollision(fcol, macro_LB, grad_LB);
    ComputeFcol_BC(fcol);

    /**  END COLLISION
      */
    ComputeTransferLBM2FV(fcol, macro_NS); // Copying LBM macro field onto NS field: Only one output required.

    // Filter NS solution
    FilterX7_NS(macro_NS, current_slot); FilterY7_NS(macro_NS, current_slot);


    end_t = clock();
    total_t = (double) (end_t-start_t)/CLOCKS_PER_SEC;
    cpu += total_t;

    if(!(s%period)){
      //time
      fprintf(stdout, "time-step : %d\n", s);

      double physical_time = (s * DX)/(sqrt(3)*Csound);
      fprintf(stdout, "simulated time: %f s \n", physical_time);

      ComputeMassNS(macro_NS);
      ComputeMassLB(fcol);

      DumpMacroNS(macro_NS, index);
      DumpMacroLB(macro_LB, index);
      //SigmaASCII(sigma);

      index++;
    }

  }
  return 0;
}

void ComputeCorrection(double **macro, double **grad_LB){
    int i, j, l, inode, ilink, xlink, ylink;
    double rho, u, v;

    int off = current_slot * xMaxp_LB*yMaxp_LB;
    for(i=1; i<=xMax_LB; i++)
        for(j=2; j<=(yMax_LB-1); j++){

            double u3x=0.0, v3y=0.0;

            for(l=0; l<NPOP; l++){

                xlink = i+ex[l];
                ylink = j+ey[l];
                if(xlink<1) xlink=xMax_LB;
                if(xlink>xMax_LB) xlink=1;
// 16.03 check because doamin is not periodic in y anymore 
                //if(ylink<1) ylink=yMax_rootgrid; 
                //if(ylink>yMax_rootgrid) ylink=1;


                ilink = off + IDM_LB(xlink, ylink);

                rho = macro[ilink][RHO];
                u = macro[ilink][RHOUX] / rho;
                v = macro[ilink][RHOUY] / rho;

                u3x += (cc[l]/4.0) * ex[l] * (rho*u*u*u);
                v3y += (cc[l]/4.0) * ey[l] * (rho*v*v*v);
                }

            inode = off + IDM_LB(i, j);

            grad_LB[inode][UX3X] = u3x;
            grad_LB[inode][UY3Y] = v3y;
      }
}