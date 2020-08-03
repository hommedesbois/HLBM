 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include <math.h>
#include "ludwig.h"

void InterfaceFlux(double **Macro, double ***RK, double ***Surf, int step){
  int m, i, s;
  int id_0, id_1;
  int off;

  off = current_slot * xMaxp_NS * yMaxp_NS;
  //fprintf(stdout, "step %d\n", step);
    switch(step){
    case 0:
        for(m=0; m<=6; m++) // macroskopic variable
            for(i=1; i<=xMax_NS; i++)
                for(s=0; s<=1; s++){

                if(s==0){
                    id_0 = off + IDM_NS(i,0);
                    id_1 = off + IDM_NS(i,1);
                }
                else if(s==1){
                    id_0 = off + IDM_NS(i,yMax_NS);
                    id_1 = off + IDM_NS(i,(yMax_NS+1));
                }
            Surf[s][i][m] = 0.5 * (Macro[id_0][m] + Macro[id_1][m]);
            //fprintf(stdout, "surf %g s %d i %d m %d \n", Surf[s][i][0], s, i, m);
        }
    break;

    repeat:
    case 1:
        for(m=0; m<=6; m++) // macroscopic variable
            for(i=1; i<=xMax_NS; i++)
                for(s=0; s<=1; s++){

                if(s==0){
                    id_0 =  IDM_NS(i,0);
                    id_1 =  IDM_NS(i,1);
                }
                else if(s==1){
                    id_0 = IDM_NS(i,yMax_NS);
                    id_1 = IDM_NS(i,(yMax_NS+1));
                }
            Surf[s][i][m] = 0.5 * (RK[step][id_0][m] + RK[step][id_1][m]);
        }

    break;
    case 2:
        goto repeat;
        break;
    case 3:
        goto repeat;
        break;

    }
}

void ComputeTransferLBM2FV(double *cell, double **macro_NS){
int i, j, j_rd1, l;
int id_fd, id_rd1_f;
int off_NS, off_LB;
double u, v;
double Pixx, Pixy, Piyy;
double Pxx_eq, Pxy_eq, Pyy_eq;
double Pxx_neq, Pxy_neq, Pyy_neq;

int interspace = 0;
interspace ++;

int ymin = 0;
int ymax = NY_NS;  // padded domain. cells for computation go from 1 to 100
int start_index_LB_on_NS   =  ymin + shift + interspace;
int end_index_LB_on_NS = ymax - (shift + interspace);


off_NS = current_slot * xMaxp_NS * yMaxp_NS;
off_LB =  current_slot * xMaxp_LB * yMaxp_LB;

for(i=0; i<xMaxp_NS; i++)
    for(j=start_index_LB_on_NS; j<=end_index_LB_on_NS; j++){

        id_fd = off_NS + IDM_NS(i,j);

        double rho = 0.0, rho_u = 0.0, rho_v = 0.0;
        Pixx = 0.0, Pixy = 0.0, Piyy=0; 

        j_rd1 = j - shift;


        for(l=0; l<NPOP; l++){

            id_rd1_f = off_LB * NPOP + IDF(i,j_rd1,l);

            // get macros from LB
            rho     += cell[id_rd1_f];
            rho_u   += cell[id_rd1_f] * ex[l];
            rho_v   += cell[id_rd1_f] * ey[l];

            Pixx    += cell[id_rd1_f]*ex[l]*ex[l];
            Pixy    += cell[id_rd1_f]*ex[l]*ey[l];
            Piyy    += cell[id_rd1_f]*ey[l]*ey[l]; 
        }
        u = rho_u/rho;
        v = rho_v/rho;

        Pxx_eq = rho * (1./3. + (u * u));
        Pxy_eq = rho * u * v;
        Pyy_eq = rho * (1./3. + (v * v));

        // Rescaling
        Pxx_neq = (Pixx - Pxx_eq) * (tau/tau_g);
        Pxy_neq = (Pixy - Pxy_eq) * (tau/tau_g);
        Pyy_neq = (Piyy - Pyy_eq) * (tau/tau_g);

        macro_NS[id_fd][RHO] = rho ;
        macro_NS[id_fd][RHOUX] = rho_u ;
        macro_NS[id_fd][RHOUY] = rho_v ;
        macro_NS[id_fd][PXX] = Pxx_eq + Pxx_neq ;
        macro_NS[id_fd][PXY] = Pxy_eq + Pxx_neq ;
        macro_NS[id_fd][PYX] = Pxy_eq + Pxy_neq ;
        macro_NS[id_fd][PYY] = Pyy_eq + Pyy_neq ;
        }
    }
