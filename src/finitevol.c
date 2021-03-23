 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include <math.h>
#include "ludwig.h"


void FiniteVol(double **Macro){
int i, j, b, c, d, e, k;
int id, id_p1, id_m1, id_p2, id_old, id_new;
int off0, off1;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef_m1, coef, coef_p1, coef_p2;
double rho_u_surf[3][4];
double rho_u_u_surf[5][4];
double sigma_surf[5][4];

double u[3];

off0 = current_slot * xMaxp_NS * yMaxp_NS;
off1 = other_slot   * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    id_old = off0 + IDM_NS(i, j);
    id_new = off1 + IDM_NS(i, j);

    double sum_rho_u = 0.0, P_x = 0.0, P_y = 0.0;

    id = off0 + IDM_NS(i,j);

    for(b=0; b<4; b++){

        int i2 = i + 2*nx[b];
        int j2 = j + 2*ny[b];

        if(j2 < 0)
            j2 = yMax_NS-1;

        if(j2 > (yMax_NS+1))
            j2 = 2;

        if(i2 < 0)
            i2 = xMax_NS-1;

        if(i2 > (xMax_NS+1))
            i2 = 2;

        id_p1 = off0 + IDM_NS((i + nx[b]), (j + ny[b]));
        id_m1 = off0 + IDM_NS((i - nx[b]), (j - ny[b]));
        id_p2 = off0 + IDM_NS(i2, j2);

        coef_m1 = 0.0;
        coef    = 0.5;
        coef_p1 = 0.5;
        coef_p2 = 0.0;

        for(c=1; c<=2; c++){

            u[c] = 0.5 * (Macro[id][c]/Macro[id][0] + Macro[id_p1][c]/Macro[id_p1][0]); // velocity at the interface
            double u_dot_n = u[1] * nx[b] + u[2] * ny[b];

            if(u_dot_n > 0.0){
                //
                coef_m1 = -1./8.;
                coef    =  6./8.;
                coef_p1 =  3./8.;
                coef_p2 =  0.0;
                /*/
                coef_m1 = 0.0;
                coef    = 1.0;
                coef_p1 = 0.0;
                coef_p2 = 0.0;
                //*/
                }
            else if(u_dot_n < 0.0){
                //
                coef_m1 =  0.0;
                coef    =  3./8.;
                coef_p1 =  6./8.;
                coef_p2 = -1./8.;
                /*/
                coef_m1 =  0.0;
                coef    =  0.0;
                coef_p1 =  1.0;
                coef_p2 =  0.0;
                //*/
                }
            rho_u_surf[c][b] = coef_m1 * u[c] * Macro[id_m1][(0)] + coef * u[c] * Macro[id][0] + coef_p1 * u[c] * Macro[id_p1][0] + coef_p2 * u[c] * Macro[id_p2][0];

            for(d=1; d<=2; d++){
                e = (c-1)*2 + d;
                k = e+2;
                rho_u_u_surf[e][b] = coef_m1 * u[c] * Macro[id_m1][d] + coef * u[c] * Macro[id][d] + coef_p1 * u[c] * Macro[id_p1][d] + coef_p2 * u[c] * Macro[id_p2][d];
                sigma_surf[e][b] = 0.5 * Macro[id][k] + 0.5 * Macro[id_p1][k];
                    }
                }
            sum_rho_u += nx[b]*rho_u_surf[1][b] + ny[b]*rho_u_surf[2][b];
            P_x += nx[b] * (rho_u_u_surf[1][b] + sigma_surf[1][b]) + ny[b] * (rho_u_u_surf[2][b] + sigma_surf[2][b]);
            P_y += nx[b] * (rho_u_u_surf[3][b] + sigma_surf[3][b]) + ny[b] * (rho_u_u_surf[4][b] + sigma_surf[4][b]);


            }
    Macro[id_new][0] = Macro[id_old][0] - Dt * sum_rho_u;
    Macro[id_new][1] = Macro[id_old][1] - Dt * P_x;
    Macro[id_new][2] = Macro[id_old][2] - Dt * P_y;
    }
}

void Flux(double **Macro, double ***flux, int which_slot){
int i, j, b, c, d, e, k, i2, j2;
int id, id_flux, id_p1, id_m1, id_p2;
int off;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef_m1, coef, coef_p1, coef_p2;
double rho_u_surf[3][4];
double rho_u_u_surf[5][4];
double sigma_surf[5][4];


double u[3];


off = which_slot * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    id = off + IDM_NS(i,j);
    id_flux = IDM_NS(i,j);

    for(b=0; b<4; b++){

        i2 = i + 2*nx[b];
        j2 = j + 2*ny[b];

        if(j2 < 0)
            j2 = yMax_NS-1;

        if(j2 > (yMax_NS+1))
            j2 = 2;

        if(i2 < 0)
            i2 = xMax_NS-1;

        if(i2 > (xMax_NS+1))
            i2 = 2;

        id_p1 = off + IDM_NS((i+nx[b]), (j+ny[b]));
        id_m1 = off + IDM_NS((i-nx[b]), (j-ny[b]));
        id_p2 = off + IDM_NS(i2, j2);

        for(c=1; c<=2; c++){

            u[c] = 0.5 * (Macro[id][c]/Macro[id][0] + Macro[id_p1][c]/Macro[id_p1][0]); // velocity at the interface
            double u_dot_n = u[1] * nx[b] + u[2] * ny[b];

            if(u_dot_n > 0.0){
                /*/
                coef_m1 = -1./8.;
                coef    =  6./8.;
                coef_p1 =  3./8.;
                coef_p2 =  0.0;
                /*/
                coef_m1 = 0.0;
                coef    = 0.5;
                coef_p1 = 0.5;
                coef_p2 = 0.0;
                //*/
                }
            else if(u_dot_n < 0.0){
                /*/
                coef_m1 =  0.0;
                coef    =  3./8.;
                coef_p1 =  6./8.;
                coef_p2 = -1./8.;
                /*/
                coef_m1 =  0.0;
                coef    =  0.5;
                coef_p1 =  0.5;
                coef_p2 =  0.0;
                //*/
                }
            rho_u_surf[c][b] = coef_m1 * u[c] * Macro[id_m1][(0)] + coef * u[c] * Macro[id][0] + coef_p1 * u[c] * Macro[id_p1][0] + coef_p2 * u[c] * Macro[id_p2][0];

            for(d=1; d<=2; d++){
                e = (c-1)*2 + d;
                k = 2 + e;
                rho_u_u_surf[e][b] = coef_m1 * u[c] * Macro[id_m1][(d)] + coef * u[c] * Macro[id][d] + coef_p1 * u[c] * Macro[id_p1][d] + coef_p2 * u[c] * Macro[id_p2][d];
                sigma_surf[e][b] = 0.5 * Macro[id][k] + 0.5 * Macro[id_p1][k];
                }
            }
        flux[id_flux][0][b] = rho_u_surf[1][b];
        flux[id_flux][1][b] = rho_u_surf[2][b];
        flux[id_flux][2][b] = rho_u_u_surf[1][b] + sigma_surf[1][b];
        flux[id_flux][3][b] = rho_u_u_surf[2][b] + sigma_surf[2][b];
        flux[id_flux][4][b] = rho_u_u_surf[3][b] + sigma_surf[3][b];
        flux[id_flux][5][b] = rho_u_u_surf[4][b] + sigma_surf[4][b];

//        for(k=1; k<=6; k++){
//
//        double macro = Macro[id_node][k];
//        double rm_dot_n = macro * n[b];
//        rm_dot_n = 0.0;
//
//            if((rm_dot_n > 0.0) && (k < 3)) {
//                //
//                coef_m1 = -1./8.;
//                coef    =  6./8.;
//                coef_p1 =  3./8.;
//                coef_p2 =  0.0;
//                /*/
//                coef_m1 = 0.0;
//                coef    = 1.0;
//                coef_p1 = 0.0;
//                coef_p2 = 0.0;
//                //*/
//                }
//            else if((rm_dot_n < 0.0) && (k < 3)){
//                //
//                coef_m1 =  0.0;
//                coef    =  3./8.;
//                coef_p1 =  6./8.;
//                coef_p2 = -1./8.;
//                /*/
//                coef_m1 =  0.0;
//                coef    =  0.0;
//                coef_p1 =  1.0;
//                coef_p2 =  0.0;
//                //*/
//                }
//            else{
//                coef_m1 = 0.0;
//                coef    = 0.5;
//                coef_p1 = 0.5;
//                coef_p2 = 0.0;
//                }

//            flux[id_flux][(k-1)][b] = coef_m1 * Macro[id_m1][(k)] + coef * Macro[id_node][k] + coef_p1 * Macro[id_p1][k] + coef_p2 * Macro[id_p2][k];
            //flux[id_flux][(k-1)][b] = coef * Macro[id_node][k] + coef_p1 * Macro[id_p1][k];
            }
        }
    }

void FiniteVolHeun(double **Macro, double ***flux){
int i, j, b, c, d, e, k;
int id, id_p1, id_m1, id_p2, id_new, id_flux;

int off0, off1;

/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef_m1, coef, coef_p1, coef_p2;
double rho_u_surf[3][4];
double rho_u_u_surf[5][4];
double sigma_surf[5][4];

double u[3];

off0 = current_slot * xMaxp_NS * yMaxp_NS;
off1 = other_slot   * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    id = off0 + IDM_NS(i,j);
    id_new = off1 + IDM_NS(i,j);
    id_flux = IDM_NS(i,j);


    double rho_u = 0.0, P_x = 0.0, P_y = 0.0;


    for(b=0; b<4; b++){

        int i2 = i + 2*nx[b];
        int j2 = j + 2*ny[b];

        if(j2 < 0)
            j2 = yMax_NS-1;

        if(j2 > (yMax_NS+1))
            j2 = 2;

        if(i2 < 0)
            i2 = xMax_NS-1;

        if(i2 > (xMax_NS+1))
            i2 = 2;

        id_p1 = off0 + IDM_NS((i + nx[b]), (j + ny[b]));
        id_m1 = off0 + IDM_NS((i - nx[b]), (j - ny[b]));
        id_p2 = off0 + IDM_NS(i2, j2);

        coef_m1 = 0.0;
        coef    = 0.5;
        coef_p1 = 0.5;
        coef_p2 = 0.0;

        for(c=1; c<=2; c++){

            u[c] = 0.5 * (Macro[id][c]/Macro[id][0] + Macro[id_p1][c]/Macro[id_p1][0]); // velocity at the interface
            double u_dot_n = u[1] * nx[b] + u[2] * ny[b];

            if(u_dot_n > 0.0){
                //
                coef_m1 = -1./8.;
                coef    =  6./8.;
                coef_p1 =  3./8.;
                coef_p2 =  0.0;
                /*/
                coef_m1 = 0.0;
                coef    = 1.0;
                coef_p1 = 0.0;
                coef_p2 = 0.0;
                //*/
                }
            else if(u_dot_n < 0.0){
                //
                coef_m1 =  0.0;
                coef    =  3./8.;
                coef_p1 =  6./8.;
                coef_p2 = -1./8.;
                /*/
                coef_m1 =  0.0;
                coef    =  0.0;
                coef_p1 =  1.0;
                coef_p2 =  0.0;
                //*/
                }
            rho_u_surf[c][b] = coef_m1 * u[c] * Macro[id_m1][(0)] + coef * u[c] * Macro[id][0] + coef_p1 * u[c] * Macro[id_p1][0] + coef_p2 * u[c] * Macro[id_p2][0];

            for(d=1; d<=2; d++){
                e = (c-1)*2 + d;
                k=e+2;
                rho_u_u_surf[e][b] = coef_m1 * u[c] * Macro[id_m1][(d)] + coef * u[c] * Macro[id][d] + coef_p1 * u[c] * Macro[id_p1][d] + coef_p2 * u[c] * Macro[id_p2][d];
                sigma_surf[e][b] = 0.5 * Macro[id][k] + 0.5 * Macro[id_p1][k];
                }
            }
            rho_u += 0.5 * nx[b]* (rho_u_surf[1][b] + flux[id_flux][0][b]) + 0.5 * ny[b]*(rho_u_surf[2][b] + flux[id_flux][1][b]);
            P_x += 0.5 * nx[b] * ((rho_u_u_surf[1][b] + sigma_surf[1][b]) + flux[id_flux][2][b]) + 0.5 * ny[b]*((rho_u_u_surf[2][b]+ sigma_surf[2][b]) + flux[id_flux][3][b]);
            P_y += 0.5 * nx[b] * ((rho_u_u_surf[3][b] + sigma_surf[3][b]) + flux[id_flux][4][b]) + 0.5 * ny[b]*((rho_u_u_surf[4][b]+ sigma_surf[4][b]) + flux[id_flux][5][b]);
            }

    Macro[id_new][0] = Macro[id][0] - Dt * rho_u;
    Macro[id_new][1] = Macro[id][1] - Dt * P_x;
    Macro[id_new][2] = Macro[id][2] - Dt * P_y;
    }
}

void FiniteVolRK4_1(double **Macro, double ***RK, double ***QRK, double ***Surf){
int i, j, k, b;
int id, id_rk, id_p1;
int off;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off = current_slot * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id = off + IDM_NS(i,j);
    id_rk = IDM_NS(i,j);

    for(b=0; b<4; b++){
        for(k=RHOUX; k<=PYY; k++){

            id_p1 = off + IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

		/*
		***** rechange if periodicity removed ***
            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = Surf[0][i][k];
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = Surf[1][i][k];
            else
            	*/
	    Macro_Surf[(k-1)][b] = coef * Macro[id][k] + coef_p1 * Macro[id_p1][k]; // Interpolation of marco on surface

            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b]; // slope
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[0][id_rk][0] = sum_rho_u;
    QRK[0][id_rk][1] = sum_Px;
    QRK[0][id_rk][2] = sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */

    RK[1][id_rk][0] = Macro[id][RHO] - 0.5 * Dt * QRK[0][id_rk][0];
    RK[1][id_rk][1] = Macro[id][RHOUX] - 0.5 * Dt * QRK[0][id_rk][1];
    RK[1][id_rk][2] = Macro[id][RHOUY] - 0.5 * Dt * QRK[0][id_rk][2];
    }
}

void FiniteVolRK4_23(double **Macro, double ***RK, double ***QRK, double ***Surf, int step){
int i, j, k, b;
int id, id_rk, id_p1;
int off;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off = current_slot * xMaxp_NS * yMaxp_NS;
step = step -1;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id = off + IDM_NS(i,j);
    id_rk = IDM_NS(i,j);

    for(b=0; b<4; b++){
        for(k=RHOUX; k<=PYY; k++){

            id_p1 = IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = coef * RK[step][id_rk][k] + coef_p1 * RK[step][IDM_NS((i + nx[b]), yMax_NS)][k]; // Surf[0][i][k]; --> change back if periodicity removed
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = coef * RK[step][id_rk][k] + coef_p1 * RK[step][IDM_NS((i + nx[b]), 1)][k];       // Surf[1][i][k]; 
            else
            Macro_Surf[(k-1)][b] = coef * RK[step][id_rk][k] + coef_p1 * RK[step][id_p1][k];
            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b];
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[step][id_rk][0] = sum_rho_u;
    QRK[step][id_rk][1] = sum_Px;
    QRK[step][id_rk][2] = sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */
    int sp1 = step+1;
    RK[sp1][id_rk][0] = Macro[id][RHO] - 1./(3.-step) * Dt * QRK[step][id_rk][0];
    RK[sp1][id_rk][1] = Macro[id][RHOUX] - 1./(3.-step) * Dt * QRK[step][id_rk][1];
    RK[sp1][id_rk][2] = Macro[id][RHOUY] - 1./(3.-step) * Dt * QRK[step][id_rk][2];
    }
}

void FiniteVolRK4_4(double **Macro, double ***RK, double ***QRK, double ***Surf, double *s){

int i, j, k, b;
int id_rk, id_p1, id_s, id_old, id_new;
int off0, off1;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off0 = current_slot * xMaxp_NS * yMaxp_NS;
off1 = other_slot * xMaxp_NS * yMaxp_NS;


for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id_rk = IDM_NS(i,j);
    id_old = off0 + IDM_NS(i,j);
    id_new = off1 + IDM_NS(i,j);
    id_s = (j-1)*xMax_NS + (i-1);

    for(b=0; b<4; b++){
        for(k=1; k<=6; k++){

            id_p1 = IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = coef * RK[3][id_rk][k] + coef_p1 * RK[3][IDM_NS((i + nx[b]), yMax_NS)][k]; // Surf[0][i][k];
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = coef * RK[3][id_rk][k] + coef_p1 * RK[3][IDM_NS((i + nx[b]), 1)][k]; // Surf[1][i][k];
            else
            Macro_Surf[(k-1)][b] = coef * RK[3][id_rk][k] + coef_p1 * RK[3][id_p1][k];
            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b];
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[3][id_rk][0] = sum_rho_u;
    QRK[3][id_rk][1] = sum_Px;
    QRK[3][id_rk][2] = sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */
//    if(j>xMax_NS){
//    Macro[id_new][0] = Macro[id_old][0] - Dt/6. * (QRK[0][id_rk][0] + 2*QRK[1][id_rk][0] + 2*QRK[2][id_rk][0] + QRK[3][id_rk][0]);
//    Macro[id_new][1] = Macro[id_old][1] - Dt/6. * (QRK[0][id_rk][1] + 2*QRK[1][id_rk][1] + 2*QRK[2][id_rk][1] + QRK[3][id_rk][1]);
//    Macro[id_new][2] = Macro[id_old][2] - Dt/6. * (QRK[0][id_rk][2] + 2*QRK[1][id_rk][2] + 2*QRK[2][id_rk][2] + QRK[3][id_rk][2]);

    s[id_s] = 0.0; // sponge zones deactivated

    Macro[id_new][RHO] = Macro[id_old][RHO] - Dt/6. * (QRK[0][id_rk][0] + 2*QRK[1][id_rk][0] + 2*QRK[2][id_rk][0] + QRK[3][id_rk][0]) - Dt * s[id_s]* (Macro[id_old][RHO]-rho0);
    Macro[id_new][RHOUX] = Macro[id_old][RHOUX] - Dt/6. * (QRK[0][id_rk][1] + 2*QRK[1][id_rk][1] + 2*QRK[2][id_rk][1] + QRK[3][id_rk][1]) - Dt * s[id_s]* (Macro[id_old][RHOUX]-rho0*U0/(sqrt(3)*Csound));
    Macro[id_new][RHOUY] = Macro[id_old][RHOUY] - Dt/6. * (QRK[0][id_rk][2] + 2*QRK[1][id_rk][2] + 2*QRK[2][id_rk][2] + QRK[3][id_rk][2]) - Dt * s[id_s]* (Macro[id_old][RHOUY]-rho0*V0/(sqrt(3)*Csound));

    }
}


void FiniteVolRK3_1(double **Macro, double ***RK, double ***QRK, double ***Surf){
int i, j, k, b;
int id, id_rk, id_p1;
int off;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off = current_slot * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id = off + IDM_NS(i,j);
    id_rk = IDM_NS(i,j);

    for(b=0; b<4; b++){
        for(k=RHOUX; k<=PYY; k++){

            id_p1 = off + IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

		/*
		***** rechange if periodicity removed ***
            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = Surf[0][i][k];
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = Surf[1][i][k];
            else
            	*/
	    Macro_Surf[(k-1)][b] = coef * Macro[id][k] + coef_p1 * Macro[id_p1][k]; // Interpolation of marco on surface

            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b]; // slope
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[0][id_rk][0] = -sum_rho_u;
    QRK[0][id_rk][1] = -sum_Px;
    QRK[0][id_rk][2] = -sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */

    RK[1][id_rk][0] = Macro[id][RHO] + Dt * QRK[0][id_rk][0];
    RK[1][id_rk][1] = Macro[id][RHOUX] + Dt * QRK[0][id_rk][1];
    RK[1][id_rk][2] = Macro[id][RHOUY] + Dt * QRK[0][id_rk][2];
    }
}


void FiniteVolRK3_2(double **Macro, double ***RK, double ***QRK, double ***Surf){
int i, j, k, b;
int id, id_rk, id_p1;
int off;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off = current_slot * xMaxp_NS * yMaxp_NS;

for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id = off + IDM_NS(i,j);
    id_rk = IDM_NS(i,j);

    for(b=0; b<4; b++){
        for(k=RHOUX; k<=PYY; k++){

            id_p1 = IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = coef * RK[1][id_rk][k] + coef_p1 * RK[1][IDM_NS((i + nx[b]), yMax_NS)][k]; // Surf[0][i][k]; --> change back if periodicity removed
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = coef * RK[1][id_rk][k] + coef_p1 * RK[1][IDM_NS((i + nx[b]), 1)][k];       // Surf[1][i][k]; 
            else
            Macro_Surf[(k-1)][b] = coef * RK[1][id_rk][k] + coef_p1 * RK[1][id_p1][k];
            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b];
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[1][id_rk][0] = -sum_rho_u;
    QRK[1][id_rk][1] = -sum_Px;
    QRK[1][id_rk][2] = -sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */
    RK[2][id_rk][0] = 0.75 * Macro[id][RHO] +   0.25 * RK[1][id_rk][0] + 0.25 * Dt * QRK[1][id_rk][0];
    RK[2][id_rk][1] = 0.75 * Macro[id][RHOUX] + 0.25 * RK[1][id_rk][1] + 0.25 * Dt * QRK[1][id_rk][1];
    RK[2][id_rk][2] = 0.75 * Macro[id][RHOUY] + 0.25 * RK[1][id_rk][2] + 0.25 * Dt * QRK[1][id_rk][2];
    }
}

void FiniteVolRK3_3(double **Macro, double ***RK, double ***QRK, double ***Surf, double *s){

int i, j, k, b;
int id_rk, id_p1, id_s, id_old, id_new;
int off0, off1;
/*
 *  Surface normals */
/*                 E   W  N   S
                   |   |  |   |            */
const int nx[4] = {1, -1, 0,  0};
const int ny[4] = {0,  0, 1, -1};

double coef, coef_p1;
double Macro_Surf[6][4];

off0 = current_slot * xMaxp_NS * yMaxp_NS;
off1 = other_slot * xMaxp_NS * yMaxp_NS;


for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

    double sum_rho_u = 0.0, sum_Px = 0.0, sum_Py = 0.0;

    id_rk = IDM_NS(i,j);
    id_old = off0 + IDM_NS(i,j);
    id_new = off1 + IDM_NS(i,j);
    id_s = (j-1)*xMax_NS + (i-1);

    for(b=0; b<4; b++){
        for(k=1; k<=6; k++){

            id_p1 = IDM_NS((i + nx[b]), (j + ny[b]));

            coef    = 0.5;
            coef_p1 = 0.5;

            if(j==1 && b==3)
            Macro_Surf[(k-1)][b] = coef * RK[2][id_rk][k] + coef_p1 * RK[2][IDM_NS((i + nx[b]), yMax_NS)][k]; // Surf[0][i][k];
            else if(j==yMax_NS && b==2)
            Macro_Surf[(k-1)][b] = coef * RK[2][id_rk][k] + coef_p1 * RK[2][IDM_NS((i + nx[b]), 1)][k]; // Surf[1][i][k];
            else
            Macro_Surf[(k-1)][b] = coef * RK[2][id_rk][k] + coef_p1 * RK[2][id_p1][k];
            }

        sum_rho_u += nx[b]*Macro_Surf[0][b] + ny[b]*Macro_Surf[1][b];
        sum_Px    += nx[b]*Macro_Surf[2][b] + ny[b]*Macro_Surf[3][b];
        sum_Py    += nx[b]*Macro_Surf[4][b] + ny[b]*Macro_Surf[5][b];
        }

    QRK[2][id_rk][0] = -sum_rho_u;
    QRK[2][id_rk][1] = -sum_Px;
    QRK[2][id_rk][2] = -sum_Py;

    /*
    RK[0] = Macro
    RK[0] = Macro
    RK[0] = Macro
    */
//    if(j>xMax_NS){
//    Macro[id_new][0] = Macro[id_old][0] - Dt/6. * (QRK[0][id_rk][0] + 2*QRK[1][id_rk][0] + 2*QRK[2][id_rk][0] + QRK[3][id_rk][0]);
//    Macro[id_new][1] = Macro[id_old][1] - Dt/6. * (QRK[0][id_rk][1] + 2*QRK[1][id_rk][1] + 2*QRK[2][id_rk][1] + QRK[3][id_rk][1]);
//    Macro[id_new][2] = Macro[id_old][2] - Dt/6. * (QRK[0][id_rk][2] + 2*QRK[1][id_rk][2] + 2*QRK[2][id_rk][2] + QRK[3][id_rk][2]);

    s[id_s] = 0.0; // sponge zones deactivated

    Macro[id_new][RHO] =   1./3. * Macro[id_old][RHO]   + 2./3.*RK[2][id_rk][0] + 2./3.*Dt*QRK[2][id_rk][0] - Dt * s[id_s]* (Macro[id_old][RHO]-rho0);
    Macro[id_new][RHOUX] = 1./3. * Macro[id_old][RHOUX] + 2./3.*RK[2][id_rk][1] + 2./3.*Dt*QRK[2][id_rk][1] - Dt * s[id_s]* (Macro[id_old][RHOUX]-rho0*U0/(sqrt(3)*Csound));
    Macro[id_new][RHOUY] = 1./3. * Macro[id_old][RHOUY] + 2./3.*RK[2][id_rk][2] + 2./3.*Dt*QRK[2][id_rk][2] - Dt * s[id_s]* (Macro[id_old][RHOUY]-rho0*V0/(sqrt(3)*Csound));

    }
}
