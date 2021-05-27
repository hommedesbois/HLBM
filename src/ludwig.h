/*
 *  ludwig.h
 *
 *  Created by Emmanuel Leveque and Tobias Horstmann
 *  Copyright 2016 CNRS. All rights reserved.
 *
 */
#define NSTEPS 1001
#define period 100

#define alpha 0.2 // filter coefficient

/* 
 * Time-marching schemes
 * 0 -> 4 stage Runge-Kutta
 * 1 -> 3 stage Runge-Kutta
 * 2 -> 2 stage Runge-Kutta
 * 3 -> Heun (2 stage RK)
 */

#define time_marching 0
#define QUICK 0 // QUICK or central difference

#define sigma_hrr 0.98 // coefficient for hybrid recursive regularization 
#define sigma_corr 0.0 // Mach number correction. Leave 0 for now.


enum {E = 0, W, N, S};
enum {RHO=0, RHOUX, RHOUY, PXX, PXY, PYX, PYY};
enum {PIXX=3, PIXY, PIYX, PIYY};
enum {UXX=0, UXY, UYX, UYY, UNEXT};
enum {UX3X = UNEXT , UY3Y};

#define max(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })

/*
 *
 *		Y
 *		^
 *		|
 *		|
 *		Z ----> X
 *
 *  periodic BC in all directions for root grid
 */

/*
 * CELL-CENTRED DATA POINTS
 */

 // NS
#define xMax_NS 100
#define yMax_NS 100

#define xMaxp_NS (xMax_NS + 2) // padded dimensions for periodic BC
#define yMaxp_NS (yMax_NS + 2)

// LB
#define xMax_LB 100
#define yMax_LB 50

#define xMaxp_LB (xMax_LB + 2) // padded dimensions for periodic BC and transition overlay
#define yMaxp_LB (yMax_LB + 2)

/*
 * MESH
 */

#define NX_NS (xMax_NS + 1)
#define NY_NS (yMax_NS + 1)

#define NX_LB (xMax_LB + 1)
#define NY_LB (yMax_LB + 1)


#define shift 25 // vertical displacement between NS and LB

#define i_c (i+1)/2   // shifted by 1
#define i_f i*2-1   // shifted by 1

#define NPOP 9 // model D2Q9 by default

#define N_fcol xMaxp_LB * yMaxp_LB * NPOP * 2

#define IDF(i,j,l)  (NPOP*((j)+yMaxp_LB*(i))+(l))
#define IDM_NS(i,j)     (i) * yMaxp_NS + (j)
#define IDM_LB(i,j)    (i) * yMaxp_LB + (j)

#define u2 ux*ux+uy*uy

// shared variables
extern const double w[NPOP];
extern const int ex[NPOP];
extern const int ey[NPOP];
extern const int finv[NPOP];
extern const double H2xx[NPOP], H2xy[NPOP], H2yy[NPOP];
extern const double H3xxy[NPOP], H3xyy[NPOP];
extern const double H4xxyy[NPOP];

extern int current_slot, other_slot;

extern double rho0; //  density parameter
extern const double Csound; // sound speed
extern double U0, V0; //  velocity parameter
extern double Es;

extern double Pref; // pressure reference
extern double P0; // pressure parameter
extern double viscosity; // viscosity
extern double Mach;

extern double Cs2;
extern double invCs2, invCs4, invCs6, invCs8;

extern double Length;
extern double pos;
extern double DX;
extern double Dt;

extern double omega_g, tau, tau_g;

