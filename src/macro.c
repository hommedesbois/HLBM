 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include <math.h>
#include "ludwig.h"

/* Calculation of the macroscopic variables */

void ComputeMacroFromF(double *cell, double **macro_NS, double **macro_LB){

int i, j, l, j_fd, id_rd1_f, id_rd1_m, id_fd;
double feq_1, feq_2, feq_3, feq_4, feq_5, feq_6, feq_7, feq_8;
double rho, u, v, ux, uy, rho_u, rho_v;
double Pxx, Pxy, Pyy;
double Pxx_eq_t0, Pxy_eq_t0, Pyy_eq_t0;
double Pxx_eq_t1, Pxy_eq_t1, Pyy_eq_t1;

int off_NS_t1 = current_slot * xMaxp_NS * yMaxp_NS;
int off_NS_t0 = other_slot * xMaxp_NS * yMaxp_NS;
int off_LB_t1 = current_slot * xMaxp_LB * yMaxp_LB;
int off_LB_t0 = other_slot * xMaxp_LB * yMaxp_LB;

for(i=1; i<=xMax_LB; i++)
  for(j=0; j<=NY_LB; j++){

     j_fd       = j + shift;
     id_fd      = off_NS_t1    + IDM_NS(i, j_fd);
     id_rd1_m   = off_LB_t1   + IDM_LB(i,j);

       if(j==0){ // macroscopic data from fluid simulation domain (NS)

          
          /*
          *  COMPUTATION  OF P^eq at previous time-step
          */
          
          // indices that point to valid data in the LB domain
          int id_1 = off_LB_t0 + IDM_LB(i, 1);
          int id_5 = off_LB_t0 + IDM_LB(i+1, 1);
          int id_6 = off_LB_t0 + IDM_LB(i-1, 1);

          // indices that point to data not availbe in the LB domain and therefore taken from NS
          int id_2 = off_NS_t0 + IDM_NS(i, (shift - 1));
          int id_3 = off_NS_t0 + IDM_NS(i+1, shift);
          int id_4 = off_NS_t0 + IDM_NS(i-1, shift);
          int id_7 = off_NS_t0 + IDM_NS(i+1, (shift-1));
          int id_8 = off_NS_t0 + IDM_NS(i-1, (shift-1));

          
          rho = macro_LB[id_1][0], ux = macro_LB[id_1][1] / rho, uy = macro_LB[id_1][2] / rho;
          // feq_1 = w[1] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*uy + 4.5*uy*uy + 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));
          feq_2 = w[2] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*uy + 4.5*uy*uy - 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_2][0], ux = macro_NS[id_2][1] / rho, uy = macro_NS[id_2][2] / rho;
          // feq_2 = w[2] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*uy + 4.5*uy*uy - 0.5*uy*(-9 * uy*uy + 9 * (ux*ux + uy*uy)));
          feq_1 = w[1] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*uy + 4.5*uy*uy + 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_3][0], ux = macro_NS[id_3][1]/rho, uy = macro_NS[id_3][2]/rho;
          // feq_3 = w[3] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*ux + 4.5*ux*ux + 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));
          feq_4 = w[4] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*ux + 4.5*ux*ux - 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_4][0], ux = macro_NS[id_4][1]/rho, uy = macro_NS[id_4][2]/rho;
          // feq_4 = w[4] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*ux + 4.5*ux*ux - 0.5*ux*(-9 * ux*ux + 9 * (ux*ux + uy*uy)));
          feq_3 = w[3] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*ux + 4.5*ux*ux + 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));

          rho = macro_LB[id_5][0], ux = macro_LB[id_5][1] / rho, uy = macro_LB[id_5][2] / rho;
          // feq_5 = w[5] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) + 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));
          feq_8 = w[8] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) - 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_LB[id_6][0], ux = macro_LB[id_6][1] / rho, uy = macro_LB[id_6][2] / rho;
          // feq_6 = w[6] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux - uy) - 4.5*(ux - uy)*(ux - uy) - 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) + 9 * (ux*ux + uy*uy)));
          feq_7 = w[7] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) + 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_7][0], ux = macro_NS[id_7][1] / rho, uy = macro_NS[id_7][2] / rho;
          // feq_7 = w[7] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) + 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));
          feq_6 = w[6] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) - 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_8][0], ux = macro_NS[id_8][1] / rho, uy = macro_NS[id_8][2] / rho;
          //feq_8 = w[8] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux + uy) - 4.5*(ux + uy)*(ux + uy) - 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) + 9 * (ux*ux + uy*uy)));
          feq_5 = w[5] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) + 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));


          Pxx_eq_t0 = feq_3 + feq_4 + feq_5 + feq_6 + feq_7 + feq_8;
          Pxy_eq_t0 = feq_5 - feq_6 - feq_7 + feq_8;
          Pyy_eq_t0 = feq_1 + feq_2 + feq_5 + feq_6 + feq_7 + feq_8;

         /*
          *  COMPUTATION OF MACRO at current time step
          */

          rho         = macro_NS[id_fd][0];
          rho_u       = macro_NS[id_fd][1];
          rho_v       = macro_NS[id_fd][2];

          u = rho_u / rho;
          v = rho_v / rho;

          Pxx_eq_t1 = rho * (1./3. + (u * u));
          Pxy_eq_t1 = rho * u * v;
          Pyy_eq_t1 = rho * (1./3. + (v * v));

          //   P reconstructed from Peq at two different time steps. Emmanuel's formula.
          Pxx = Pxx_eq_t1 - tau_g * (Pxx_eq_t1 - Pxx_eq_t0);
          Pxy = Pxy_eq_t1 - tau_g * (Pxy_eq_t1 - Pxy_eq_t0);
          Pyy = Pyy_eq_t1 - tau_g * (Pyy_eq_t1 - Pyy_eq_t0);
          /*/

          // only equilibrium part taken into account
          Pxx = Pxx_eq_t1;
          Pxy = Pxy_eq_t1;
          Pyy = Pyy_eq_t1;
          //*/

          //fprintf(stdout, " Pxx = %g Pxy  = %g Pyy = %g \n", Pxx, Pxy, Pyy);

       }else if(j==NY_LB){

          int trans_top = shift + yMax_LB + 1;
          //if(i==1)
          //fprintf(stdout, "trans_bottom %d trans_top %d \n", shift, trans_top);
          // indices that point to valid data in the LB domain
          int id_2 = off_LB_t0 + IDM_LB(i, yMax_LB);
          int id_7 = off_LB_t0 + IDM_LB(i+1, yMax_LB);
          int id_8 = off_LB_t0 + IDM_LB(i-1, yMax_LB);

          // indices that point to data not availbe in the LB domain and therefore taken from NS

          // Ich glaube das muss plus sein
          int id_1 = off_NS_t0 + IDM_NS(i, (trans_top + 1));
          int id_3 = off_NS_t0 + IDM_NS(i+1, trans_top);
          int id_4 = off_NS_t0 + IDM_NS(i-1, trans_top);
          int id_5 = off_NS_t0 + IDM_NS(i+1, (trans_top + 1));
          int id_6 = off_NS_t0 + IDM_NS(i-1, (trans_top + 1));


          rho = macro_NS[id_1][0], ux = macro_NS[id_1][1] / rho, uy = macro_NS[id_1][2] / rho;
          // feq_1 = w[1] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*uy + 4.5*uy*uy + 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));
          feq_2 = w[2] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*uy + 4.5*uy*uy - 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));

          rho = macro_LB[id_2][0], ux = macro_LB[id_2][1] / rho, uy = macro_LB[id_2][2] / rho;
          // feq_2 = w[2] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*uy + 4.5*uy*uy - 0.5*uy*(-9 * uy*uy + 9 * (ux*ux + uy*uy)));
          feq_1 = w[1] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*uy + 4.5*uy*uy + 0.5*uy*(9 * uy*uy - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_3][0], ux = macro_NS[id_3][1]/rho, uy = macro_NS[id_3][2]/rho;
          // feq_3 = w[3] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*ux + 4.5*ux*ux + 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));
          feq_4 = w[4] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*ux + 4.5*ux*ux - 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_4][0], ux = macro_NS[id_4][1]/rho, uy = macro_NS[id_4][2]/rho;
          // feq_4 = w[4] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*ux + 4.5*ux*ux - 0.5*ux*(-9 * ux*ux + 9 * (ux*ux + uy*uy)));
          feq_3 = w[3] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*ux + 4.5*ux*ux + 0.5*ux*(9 * ux*ux - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_5][0], ux = macro_NS[id_5][1] / rho, uy = macro_NS[id_5][2] / rho;
          // feq_5 = w[5] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) + 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));
          feq_8 = w[8] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) - 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_NS[id_6][0], ux = macro_NS[id_6][1] / rho, uy = macro_NS[id_6][2] / rho;
          // feq_6 = w[6] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux - uy) - 4.5*(ux - uy)*(ux - uy) - 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) + 9 * (ux*ux + uy*uy)));
          feq_7 = w[7] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) + 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_LB[id_7][0], ux = macro_LB[id_7][1] / rho, uy = macro_LB[id_7][2] / rho;
          // feq_7 = w[7] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) + 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));
          feq_6 = w[6] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux - uy) + 4.5*(ux - uy)*(ux - uy) - 0.5*(ux - uy)*(9 * (ux - uy)*(ux - uy) - 9 * (ux*ux + uy*uy)));

          rho = macro_LB[id_8][0], ux = macro_LB[id_8][1] / rho, uy = macro_LB[id_8][2] / rho;
          //feq_8 = w[8] * rho*(1.0 - 1.5*(ux*ux + uy*uy) - 3.0*(ux + uy) - 4.5*(ux + uy)*(ux + uy) - 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) + 9 * (ux*ux + uy*uy)));
          feq_5 = w[5] * rho*(1.0 - 1.5*(ux*ux + uy*uy) + 3.0*(ux + uy) + 4.5*(ux + uy)*(ux + uy) + 0.5*(ux + uy)*(9 * (ux + uy)*(ux + uy) - 9 * (ux*ux + uy*uy)));


          Pxx_eq_t0 = feq_3 + feq_4 + feq_5 + feq_6 + feq_7 + feq_8;
          Pxy_eq_t0 = feq_5 - feq_6 - feq_7 + feq_8;
          Pyy_eq_t0 = feq_1 + feq_2 + feq_5 + feq_6 + feq_7 + feq_8;


          rho         = macro_NS[id_fd][0];
          rho_u       = macro_NS[id_fd][1];
          rho_v       = macro_NS[id_fd][2];

          u = rho_u / rho;
          v = rho_v / rho;

          Pxx_eq_t1 = rho * (1./3. + (u * u));
          Pxy_eq_t1 = rho * u * v;
          Pyy_eq_t1 = rho * (1./3. + (v * v));

          //
          //P reconstructed from Peq at two different time steps. Emmanuel's formula.
          Pxx = Pxx_eq_t1 - tau_g * (Pxx_eq_t1 - Pxx_eq_t0);
          Pxy = Pxy_eq_t1 - tau_g * (Pxy_eq_t1 - Pxy_eq_t0);
          Pyy = Pyy_eq_t1 - tau_g * (Pyy_eq_t1 - Pyy_eq_t0);
          
          /* rho         = macro_NS[id_fd][0];
          rho_u       = macro_NS[id_fd][1];
          rho_v       = macro_NS[id_fd][2];

          u = rho_u / rho;
          v = rho_v / rho;

          // for now only equilibrium part!
          Pxx = rho * (1./3. + (u * u));
          Pxy = rho * u * v;
          Pyy = rho * (1./3. + (v * v)); */

       }else{ // moments of f
       rho = 0.0, rho_u = 0.0, rho_v = 0.0;
       Pxx = 0.0, Pxy = 0.0, Pyy = 0.0;
            for(l=0; l<NPOP; l++){

            id_rd1_f = off_LB_t1 * NPOP + IDF(i,j,l);

            rho         += cell[id_rd1_f];
            rho_u       += cell[id_rd1_f]*ex[l];
            rho_v       += cell[id_rd1_f]*ey[l];
            Pxx         += cell[id_rd1_f]*ex[l]*ex[l];
            Pxy         += cell[id_rd1_f]*ex[l]*ey[l];
            Pyy         += cell[id_rd1_f]*ey[l]*ey[l];

            }
        }
        macro_LB[id_rd1_m][RHO] = rho;
        macro_LB[id_rd1_m][RHOUX] = rho_u;
        macro_LB[id_rd1_m][RHOUY] = rho_v;
        macro_LB[id_rd1_m][PIXX] = Pxx;
        macro_LB[id_rd1_m][PIXY] = Pxy;
        macro_LB[id_rd1_m][PIYX] = Pxy;
        macro_LB[id_rd1_m][PIYY] = Pyy;
    } // end i,j  loop
} // end function


void ComputeGrad(double **macro, double **macro_NS, double **grad){

int i, j, id, idxp, idxm, idyp, idym;

int off_NS_t1 = current_slot * xMaxp_NS * yMaxp_NS;
// int off_NS_t0 = other_slot * xMaxp_NS * yMaxp_NS;
int off_LB_t1 = current_slot * xMaxp_LB * yMaxp_LB;
// int off_LB_t0 = other_slot * xMaxp_LB * yMaxp_LB;

int trans_top = shift + yMax_LB + 1;

for(i=1; i<=xMax_LB; i++)
  for(j=0; j<=NY_LB; j++){

      id = off_LB_t1 + IDM_LB(i,j);
      idxp = off_LB_t1 + IDM_LB(i+1,j);
      idxm = off_LB_t1 + IDM_LB(i-1,j);
      idyp = off_LB_t1 + IDM_LB(i,j+1);
      idym = off_LB_t1 + IDM_LB(i,j-1);

      // if(i==1) // brqcuhe ich erstmal nicht da BC für marco in x    
      //    idxm = off_LB_t1 + IDM_LB(xMax_LB,j);
      // if(i=xMax_LB)
      //    idxp = off_LB_t1 + IDM_LB(1,j);
      if(j==0)
         idym = off_NS_t1 + IDM_NS(i, (shift - 1));
      if(j==NY_LB)
         idyp = off_NS_t1 + IDM_NS(i, (trans_top + 1));


      grad[id][UXX] = (macro[idxp][1]/macro[idxp][0] - macro[idxm][1]/macro[idxm][0]) * 0.5;
      grad[id][UYX] = (macro[idxp][2]/macro[idxp][0] - macro[idxm][2]/macro[idxm][0]) * 0.5;
   
      if(j==0){
         grad[id][UXY] = (macro[idyp][1]/macro[idyp][0] - macro_NS[idym][1]/macro_NS[idym][0]) * 0.5;
         grad[id][UYY] = (macro[idyp][2]/macro[idyp][0] - macro_NS[idym][2]/macro_NS[idym][0]) * 0.5;
      }else if(j==NY_LB){
         grad[id][UXY] = (macro_NS[idyp][1]/macro_NS[idyp][0] - macro[idym][1]/macro[idym][0]) * 0.5;
         grad[id][UYY] = (macro_NS[idyp][2]/macro_NS[idyp][0] - macro[idym][2]/macro[idym][0]) * 0.5;
      }else{
         grad[id][UXY] = (macro[idyp][1]/macro[idyp][0] - macro[idym][1]/macro[idym][0]) * 0.5;
         grad[id][UYY] = (macro[idyp][2]/macro[idyp][0] - macro[idym][2]/macro[idym][0]) * 0.5;
      }
   }
}