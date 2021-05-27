 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include <math.h>
#include "ludwig.h"


void Closure(double **macro){
int i,j;
int off, id, idxp1, idxm1, idyp1, idym1;
double rho, ux, uy;
double duxdx, duxdy, duydx, duydy;
double Sxx, Sxy, Syx, Syy;

off = other_slot  * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

        id      = off + IDM_NS(i,j);
        idxp1   = off + IDM_NS((i+1),j);
        idxm1   = off + IDM_NS((i-1),j);
        idyp1   = off + IDM_NS(i,(j+1));
        idym1   = off + IDM_NS(i,(j-1));

        rho = macro[id][RHO];
        ux  = macro[id][RHOUX]/rho;
        uy  = macro[id][RHOUY]/rho;

        duxdx = (macro[idxp1][1]/macro[idxp1][0] -  macro[idxm1][1]/macro[idxm1][0])/2.0;
        duxdy = (macro[idyp1][1]/macro[idyp1][0] -  macro[idym1][1]/macro[idym1][0])/2.0;
        duydx = (macro[idxp1][2]/macro[idxp1][0] -  macro[idxm1][2]/macro[idxm1][0])/2.0;
        duydy = (macro[idyp1][2]/macro[idyp1][0] -  macro[idym1][2]/macro[idym1][0])/2.0;

        Sxx = duxdx;
        Sxy = 0.5 * (duxdy + duydx);
        Syx = 0.5 * (duydx + duxdy);
        Syy = duydy;
switch(time_marching){
case 0:
        macro[id][3] = rho * (1./3. + (ux * ux)) - 2./3.*tau*rho*Sxx;
        macro[id][4] = rho * ux * uy - 2./3. * tau * rho* Sxy;
        macro[id][5] = rho * uy * ux - 2./3. * tau * rho* Syx;
        macro[id][6] = rho * (1./3. + (uy * uy)) - 2./3.*tau*rho*Syy;
        break;
case 1:
        macro[id][3] = rho * (1./3. + (ux * ux)) - 2./3.*tau*rho*Sxx;
        macro[id][4] = rho * ux * uy - 2./3. * tau * rho* Sxy;
        macro[id][5] = rho * uy * ux - 2./3. * tau * rho* Syx;
        macro[id][6] = rho * (1./3. + (uy * uy)) - 2./3.*tau*rho*Syy;
        break;
case 2:
        macro[id][3] = rho * (1./3. + (ux * ux)) - 2./3.*tau*rho*Sxx;
        macro[id][4] = rho * ux * uy - 2./3. * tau * rho* Sxy;
        macro[id][5] = rho * uy * ux - 2./3. * tau * rho* Syx;
        macro[id][6] = rho * (1./3. + (uy * uy)) - 2./3.*tau*rho*Syy;
        break;
case 3:
        macro[id][3] = 1./3. * rho - 2./3.*tau*rho*Sxx;
        macro[id][4] = - 2./3. * tau * rho* Sxy;
        macro[id][5] = - 2./3. * tau * rho* Syx;
        macro[id][6] = 1./3. * rho - 2./3.*tau*rho*Syy;
        break;
        }
    }
}


void ClosureRK4(double ***RK, int step){
int i,j;
int id, idxp1, idxm1, idyp1, idym1;
double rho, ux, uy;
double duxdx, duxdy, duydx, duydy;
double Sxx, Sxy, Syx, Syy;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

        id      = IDM_NS(i,j);
        idxp1   = IDM_NS((i+1),j);
        idxm1   = IDM_NS((i-1),j);
        idyp1   = IDM_NS(i,(j+1));
        idym1   = IDM_NS(i,(j-1));

        rho = RK[step][id][0];
        ux  = RK[step][id][1]/rho;
        uy  = RK[step][id][2]/rho;


        duxdx = (RK[step][idxp1][1]/RK[step][idxp1][0] -  RK[step][idxm1][1]/RK[step][idxm1][0])/2.0;
        duxdy = (RK[step][idyp1][1]/RK[step][idyp1][0] -  RK[step][idym1][1]/RK[step][idym1][0])/2.0;
        duydx = (RK[step][idxp1][2]/RK[step][idxp1][0] -  RK[step][idxm1][2]/RK[step][idxm1][0])/2.0;
        duydy = (RK[step][idyp1][2]/RK[step][idyp1][0] -  RK[step][idym1][2]/RK[step][idym1][0])/2.0;

        Sxx = duxdx;
        Sxy = 0.5 * (duxdy + duydx);
        Syx = 0.5 * (duydx + duxdy);
        Syy = duydy;

        RK[step][id][3] = rho * (1./3. + (ux * ux)) - 2./3.*tau*rho*Sxx;
        RK[step][id][4] = rho * ux * uy - 2./3. * tau * rho* Sxy;
        RK[step][id][5] = rho * uy * ux - 2./3. * tau * rho* Syx;
        RK[step][id][6] = rho * (1./3. + (uy * uy)) - 2./3.*tau*rho*Syy;
    }
}
