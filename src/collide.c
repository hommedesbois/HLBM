 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include "ludwig.h"

void RBGKCollision(double *cell, double **macro, int coll){
int i, j, l, id, id_m, off, off_m;
double rho, rho_u, rho_v, u, v, feq, f_1;
double Pxx_eq, Pxy_eq, Pyy_eq;
double Pxx_neq, Pxy_neq, Pyy_neq;

off = current_slot * xMaxp_LB * yMaxp_LB * NPOP;
off_m = current_slot * xMaxp_LB * yMaxp_LB;

for(i=1; i<=xMax_LB; i++)
    for(j=0; j<yMaxp_LB; j++){

        id_m = off_m + IDM_LB(i,j);

        rho = macro[id_m][0];
        rho_u = macro[id_m][1];
        rho_v = macro[id_m][2];

        u = rho_u/rho;
        v = rho_v/rho;


        Pxx_eq = rho * (1./3. + (u * u));
        Pxy_eq = rho * u * v;
        Pyy_eq = rho * (1./3. + (v * v));

        Pxx_neq     = macro[id_m][3]-Pxx_eq;
        Pxy_neq     = macro[id_m][4]-Pxy_eq;
        Pyy_neq     = macro[id_m][6]-Pyy_eq;

        for(l=0; l<NPOP; l++){

            id = off + IDF(i,j,l);

            feq = w[l]*rho*(1.0 -
            1.5*(u*u + v*v) +
            3.0*(ex[l]*u + ey[l]*v) +
            4.5*(ex[l]*u + ey[l]*v) *    (ex[l]*u + ey[l]*v) +
            0.5*(ex[l]*u + ey[l]*v) * (9*(ex[l]*u + ey[l]*v) * (ex[l]*u + ey[l]*v) - 9*(u*u+v*v)));

        switch(coll){
            case 0:
            cell[id] = (1-omega_g)*cell[id] + omega_g * feq;
            break;

            case 1:
            f_1 = w[l]*0.5*invCs4*((ex[l]*ex[l]-Cs2)*Pxx_neq+2.0*ex[l]*ey[l]*Pxy_neq+(ey[l]*ey[l]-Cs2)*Pyy_neq);
            /*/
            cell[id] = feq[id] + f_1[id];
            /*/
            cell[id] = (1-omega_g)*f_1+feq;
            //*/
            break;
            }
        }
    }
}


void HRRCollision(double *cell, double **macro, double **grad, double **corr){
    int i, j, l, id, id_m, off, off_m;
    double rho, rho_u, rho_v, u, v;
    double feq, f_1, f00, f01, f02, f03, f04;
    double Pxx_eq, Pxy_eq, Pyy_eq;
    double Pxx_neq, Pxy_neq, Pyy_neq;
    double Sxx, Syy, Sxy;
    double a12xx, a12yy, a12xy;
    double a13xyy, a13xxy;
    double a14xxyy;
    double Scorr[NPOP];


    off = current_slot * xMaxp_LB * yMaxp_LB * NPOP;
    off_m = current_slot * xMaxp_LB * yMaxp_LB;


    for(i=1; i<=xMax_LB; i++)
        for(j=0; j<yMaxp_LB; j++){


            id_m = off_m + IDM_LB(i,j);

            rho = macro[id_m][RHO];
            rho_u = macro[id_m][RHOUX];
            rho_v = macro[id_m][RHOUY];

            u = rho_u/rho;
            v = rho_v/rho;

            Sxx = grad[id_m][UXX] - 0.5 * (grad[id_m][UXX] + grad[id_m][UYY]);
            Syy = grad[id_m][UYY] - 0.5 * (grad[id_m][UXX] + grad[id_m][UYY]);
            Sxy = 0.5 * (grad[id_m][UXY] + grad[id_m][UYX]);

            Pxx_eq = rho * (1./3. + (u * u));
            Pxy_eq = rho * u * v;
            Pyy_eq = rho * (1./3. + (v * v));

            Pxx_neq     = macro[id_m][PXX]-Pxx_eq;
            Pxy_neq     = macro[id_m][PXY]-Pxy_eq;
            Pyy_neq     = macro[id_m][PYY]-Pyy_eq;

            a12xx = Pxx_neq * sigma_hrr + (1.0-sigma_hrr)*(-2.0*tau_g*rho*Sxx*(1./3.)); 
            a12yy = Pyy_neq * sigma_hrr + (1.0-sigma_hrr)*(-2.0*tau_g*rho*Syy*(1./3.)); 
            a12xy = Pxy_neq * sigma_hrr + (1.0-sigma_hrr)*(-2.0*tau_g*rho*Sxy*(1./3.)); 

            a13xyy = 2.0*v*a12xy + u*a12yy;
            a13xxy = v*a12xx + 2.0*u*a12xy;
            
            a14xxyy = v*a13xxy + u*u*a12yy + 2*u*v*a12xy;

            for(l=0; l<NPOP; l++){

                id = off + IDF(i,j,l);

                f00 = 1.0;
                f01 = ex[l]*u + ey[l]*v;
                f02 = H2xx[l]*(u*u) + H2yy[l]*(v*v) + 2*H2xy[l]*u*v;
                f03 = H3xxy[l]*u*u*v + H3xyy[l]*u*v*v;
                f04 = H4xxyy[l]*u*u*v*v;
                
                //
                feq = w[l]*rho*(f00 + invCs2*f01 + 0.5*invCs4*f02 + 0.5*invCs6*f03 + 0.25*invCs8*f04);  
                
                /*/
                feq = w[l]*rho*(1.0   
                    + 3.0*(ex[l]*u + ey[l]*v)
                    + 4.5*(ex[l]*u + ey[l]*v) * (ex[l]*u + ey[l]*v) -1.5*(u*u + v*v) // f02
                    + 0.5*(ex[l]*u + ey[l]*v) * (9*(ex[l]*u + ey[l]*v) * (ex[l]*u + ey[l]*v) - 9*(u*u+v*v)) 
                    + 0.25*invCs8*f04); // f04
                //*/
                
                f_1 = w[l] * (0.5*invCs4*(H2xx[l]*a12xx + 2.0*H2xy[l]*a12xy + H2yy[l]*a12yy) 
                              + 0.5 * invCs6 * (H3xxy[l]*a13xxy + H3xyy[l]*a13xyy)
                              + 0.25 * invCs8 * H4xxyy[l]*a14xxyy); // fourth order regularization
           
                Scorr[l] = w[l]*0.5*invCs4*(H2xx[l]*grad[id_m][UX3X] + H2yy[l]*grad[id_m][UY3Y]);

                cell[id] = feq + (1-omega_g)*f_1 + sigma_corr * Scorr[l];
            }
        }
}

void ComputeFcolFromCentMomCollision(double *cell, double **macro){
  int id_0, id_1, id_2, id_3, id_4, id_5, id_6, id_7, id_8, idm;
  int i,j;
  int off_m, off_f;

  double k3_star, k4_star, k5_star, k6_star, k7_star, k8_star;
  double m3, m4, m5;
  double Rho, Ux, Uy;
  double Pxx, Pxy, Pyy;

  off_f = current_slot * xMaxp_LB * yMaxp_LB * NPOP;
  off_m = current_slot * xMaxp_LB * yMaxp_LB;

  for(i=1; i<=xMax_LB; i++)
    for(j=0; j<yMaxp_LB; j++){

      id_0 = off_f + IDF(i,j,0);
      id_1 = off_f + IDF(i,j,1);
      id_2 = off_f + IDF(i,j,2);
      id_3 = off_f + IDF(i,j,3);
      id_4 = off_f + IDF(i,j,4);
      id_5 = off_f + IDF(i,j,5);
      id_6 = off_f + IDF(i,j,6);
      id_7 = off_f + IDF(i,j,7);
      id_8 = off_f + IDF(i,j,8);

      idm = off_m + IDM_LB(i,j);

      Rho = macro[idm][0];
      Ux  = macro[idm][1]; Ux /=Rho;
      Uy  = macro[idm][2]; Uy /=Rho;
      Pxx = macro[idm][PIXX];
      Pxy = macro[idm][PIXY];
      Pyy = macro[idm][PIYY];

      m3  = Pxx + Pyy - Rho * (Ux*Ux+Uy*Uy);
      m4  = Pxx - Pyy - Rho * (Ux*Ux-Uy*Uy);
      m5  = Pxy - Rho*Ux*Uy;

      //RELAXATION OF THE MOMENTS

      m3 = m3 * (1. - omega_g) + omega_g*2.0*Rho/3.;
      m4 = m4 * (1. - omega_g);
      m5 = m5 * (1. - omega_g);

      k3_star = (2./3.)*Rho; // m3
      k4_star = m4;
      k5_star = m5;
      k6_star = 0.0;
      k7_star = 0.0;
      k8_star = -Rho*(9.*Ux*Ux*Uy*Uy-1.)/9.;

      cell[id_0] = k8_star + 2.*Ux*k7_star + 2.*Uy*k6_star - Rho*(- Ux*Ux*Uy*Uy + Ux*Ux + Uy*Uy - 1.) + k3_star*(Ux*Ux/2. + Uy*Uy/2. - 1.) - k4_star*(Ux*Ux/2. - Uy*Uy/2.) + 4.*Ux*Uy*k5_star;

      cell[id_1] = -k8_star/2. - k3_star*(Ux*Ux/4. + Uy*Uy/4. + Uy/4. - 1./4.) - k4_star*(- Ux*Ux/4. + Uy*Uy/4. + Uy/4. + 1./4.) - Ux*k7_star - k6_star*(Uy + 1./2.) - Ux*k5_star*(2.*Uy + 1.) - (Rho*Uy*(Ux*Ux - 1.)*(Uy + 1.))/2.;

      cell[id_2] = k3_star*(- Ux*Ux/4. - Uy*Uy/4. + Uy/4. + 1./4.) - k8_star/2. + k4_star*(Ux*Ux/4. - Uy*Uy/4. + Uy/4. - 1./4.) - Ux*k7_star - k6_star*(Uy - 1./2.) - Ux*k5_star*(2.*Uy - 1.) - (Rho*Uy*(Ux*Ux - 1.)*(Uy - 1.))/2.;

      cell[id_3] = k4_star*(Ux*Ux/4. + Ux/4. - Uy*Uy/4. + 1./4.) - k3_star*(Ux*Ux/4. + Ux/4. + Uy*Uy/4. - 1./4.) - k8_star/2. - Uy*k6_star - k7_star*(Ux + 1./2.) - Uy*k5_star*(2.*Ux + 1.) - (Rho*Ux*(Uy*Uy - 1.)*(Ux + 1.))/2.;

      cell[id_4] = k3_star*(- Ux*Ux/4. + Ux/4. - Uy*Uy/4. + 1./4.) - k8_star/2. - k4_star*(- Ux*Ux/4. + Ux/4. + Uy*Uy/4. - 1./4.) - Uy*k6_star - k7_star*(Ux - 1./2.) - Uy*k5_star*(2.*Ux - 1.) - (Rho*Ux*(Uy*Uy - 1.)*(Ux - 1.))/2.;

      cell[id_5] = k8_star/4. + k3_star*(Ux*Ux/8. + Ux/8. + Uy*Uy/8. + Uy/8.) - k4_star*(Ux*Ux/8. + Ux/8. - Uy*Uy/8. - Uy/8.) + k7_star*(Ux/2. + 1./4.) + k6_star*(Uy/2. + 1./4.) + (k5_star*(2.*Ux + 1.)*(2.*Uy + 1.))/4. + (Rho*Ux*Uy*(Ux + 1.)*(Uy + 1.))/4;

      cell[id_6] = k8_star/4. + k3_star*(Ux*Ux/8. - Ux/8. + Uy*Uy/8. + Uy/8.) + k4_star*(- Ux*Ux/8. + Ux/8. + Uy*Uy/8. + Uy/8.) + k7_star*(Ux/2. - 1./4.) + k6_star*(Uy/2. + 1./4.) + (k5_star*(2.*Ux - 1)*(2.*Uy + 1))/4. + (Rho*Ux*Uy*(Ux - 1.)*(Uy + 1.))/4.;

      cell[id_7] = k8_star/4. + k3_star*(Ux*Ux/8. + Ux/8. + Uy*Uy/8. - Uy/8.) - k4_star*(Ux*Ux/8. + Ux/8. - Uy*Uy/8. + Uy/8.) + k7_star*(Ux/2. + 1./4.) + k6_star*(Uy/2. - 1./4.) + (k5_star*(2.*Ux + 1.)*(2.*Uy - 1))/4. + (Rho*Ux*Uy*(Ux + 1.)*(Uy - 1.))/4.;

      cell[id_8] = k8_star/4. - k3_star*(- Ux*Ux/8. + Ux/8. - Uy*Uy/8. + Uy/8.) + k4_star*(- Ux*Ux/8. + Ux/8. + Uy*Uy/8. - Uy/8.) + k7_star*(Ux/2. - 1./4.) + k6_star*(Uy/2. - 1./4.) + (k5_star*(2.*Ux - 1.)*(2.*Uy - 1))/4. + (Rho*Ux*Uy*(Ux - 1.)*(Uy - 1.))/4.;
    }
}
