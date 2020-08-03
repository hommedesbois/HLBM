/* T. Horstmann
* Department of Aeroacoustics
* ONERA
* Chatillon
* France
*/

#include <stdio.h>
#include <math.h>
#include "ludwig.h"

void FilterX7_NS(double **Macro, int which_slot){
  int i, j, id;
  int idm3x, idm2x, idm1x, idp1x, idp2x, idp3x;
  int off;

  double Macro_filt[7];

  off = which_slot*xMaxp_NS*yMaxp_NS;

  for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

        id = off + IDM_NS(i,j);

        idm3x = off + IDM_NS(i-3,j);
        idm2x = off + IDM_NS(i-2,j);
        idm1x = off + IDM_NS(i-1,j);
        idp1x = off + IDM_NS(i+1,j);
        idp2x = off + IDM_NS(i+2,j);
        idp3x = off + IDM_NS(i+3,j);

        if(i==1){
            idm3x = off + IDM_NS(xMax_NS-2,j);
            idm2x = off + IDM_NS(xMax_NS-1,j);
            idm1x = off + IDM_NS(xMax_NS,j);
        }
        else if(i==2){
            idm3x = off + IDM_NS(xMax_NS-1,j);
            idm2x = off + IDM_NS(xMax_NS,j);
        }
        else if(i==3){
            idm3x = off + IDM_NS(xMax_NS,j);
        }
        else if(i==xMax_NS){
            idp1x = off + IDM_NS(1,j);
            idp2x = off + IDM_NS(2,j);
            idp3x = off + IDM_NS(3,j);
        }
        else if(i==(xMax_NS-1)){
            idp2x = off + IDM_NS(1,j);
            idp3x = off + IDM_NS(2,j);
        }
        else if(i==(xMax_NS-2)){
            idp3x = off + IDM_NS(1,j);
        }

        Macro_filt[0] = 1./64.*(Macro[idm3x][0] +  Macro[idp3x][0]) - 3./32.*(Macro[idm2x][0]+Macro[idp2x][0]) + 15./64.*(Macro[idm1x][0]+Macro[idp1x][0]) + 11./16.*Macro[id][0];
        Macro_filt[1] = 1./64.*(Macro[idm3x][1] +  Macro[idp3x][1]) - 3./32.*(Macro[idm2x][1]+Macro[idp2x][1]) + 15./64.*(Macro[idm1x][1]+Macro[idp1x][1]) + 11./16.*Macro[id][1];
        Macro_filt[2] = 1./64.*(Macro[idm3x][2] +  Macro[idp3x][2]) - 3./32.*(Macro[idm2x][2]+Macro[idp2x][2]) + 15./64.*(Macro[idm1x][2]+Macro[idp1x][2]) + 11./16.*Macro[id][2];
        Macro_filt[3] = 1./64.*(Macro[idm3x][3] +  Macro[idp3x][3]) - 3./32.*(Macro[idm2x][3]+Macro[idp2x][3]) + 15./64.*(Macro[idm1x][3]+Macro[idp1x][3]) + 11./16.*Macro[id][3];
        Macro_filt[4] = 1./64.*(Macro[idm3x][4] +  Macro[idp3x][4]) - 3./32.*(Macro[idm2x][4]+Macro[idp2x][4]) + 15./64.*(Macro[idm1x][4]+Macro[idp1x][4]) + 11./16.*Macro[id][4];
        Macro_filt[5] = 1./64.*(Macro[idm3x][5] +  Macro[idp3x][5]) - 3./32.*(Macro[idm2x][5]+Macro[idp2x][5]) + 15./64.*(Macro[idm1x][5]+Macro[idp1x][5]) + 11./16.*Macro[id][5];
        Macro_filt[6] = 1./64.*(Macro[idm3x][6] +  Macro[idp3x][6]) - 3./32.*(Macro[idm2x][6]+Macro[idp2x][6]) + 15./64.*(Macro[idm1x][6]+Macro[idp1x][6]) + 11./16.*Macro[id][6];
        

        Macro[id][0]=alpha * Macro_filt[0] + (1-alpha)*Macro[id][0];
        Macro[id][1]=alpha * Macro_filt[1] + (1-alpha)*Macro[id][1];
        Macro[id][2]=alpha * Macro_filt[2] + (1-alpha)*Macro[id][2];
        Macro[id][3]=alpha * Macro_filt[3] + (1-alpha)*Macro[id][3];
        Macro[id][4]=alpha * Macro_filt[4] + (1-alpha)*Macro[id][4];
        Macro[id][5]=alpha * Macro_filt[5] + (1-alpha)*Macro[id][5];
        Macro[id][6]=alpha * Macro_filt[6] + (1-alpha)*Macro[id][6];
        }
}

void FilterY7_NS(double **Macro, int which_slot){
  int i, j, id;
  int idm3y, idm2y, idm1y, idp1y, idp2y, idp3y;
  int off;

  double Macro_filt[7];


  off = which_slot*xMaxp_NS*yMaxp_NS;

  for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=yMax_NS; j++){

        id = off + IDM_NS(i,j);

        idm3y = off + IDM_NS(i,j-3);
        idm2y = off + IDM_NS(i,j-2);
        idm1y = off + IDM_NS(i,j-1);
        idp1y = off + IDM_NS(i,j+1);
        idp2y = off + IDM_NS(i,j+2);
        idp3y = off + IDM_NS(i,j+3);


        if(j==1){
            idm3y = off + IDM_NS(i,yMax_NS-2);
            idm2y = off + IDM_NS(i,yMax_NS-1);
            idm1y = off + IDM_NS(i,yMax_NS);
        }
        else if(j==2){
            idm3y = off + IDM_NS(i,yMax_NS-1);
            idm2y = off + IDM_NS(i,yMax_NS);
        }
        else if(j==3){
            idm3y = off + IDM_NS(i,yMax_NS);
        }
        else if(j==yMax_NS){
            idp1y = off + IDM_NS(i,1);
            idp2y = off + IDM_NS(i,2);
            idp3y = off + IDM_NS(i,3);
        }
        else if(j==(yMax_NS-1)){
            idp2y = off + IDM_NS(i,1);
            idp3y = off + IDM_NS(i,2);
        }
        else if(j==(yMax_NS-2)){
            idp3y = off + IDM_NS(i,1);
        }

        Macro_filt[0] = 1./64.*(Macro[idm3y][0] +  Macro[idp3y][0]) - 3./32.*(Macro[idm2y][0]+Macro[idp2y][0]) + 15./64.*(Macro[idm1y][0] + Macro[idp1y][0]) + 11./16. * Macro[id][0];
        Macro_filt[1] = 1./64.*(Macro[idm3y][1] +  Macro[idp3y][1]) - 3./32.*(Macro[idm2y][1]+Macro[idp2y][1]) + 15./64.*(Macro[idm1y][1] + Macro[idp1y][1]) + 11./16. * Macro[id][1];
        Macro_filt[2] = 1./64.*(Macro[idm3y][2] +  Macro[idp3y][2]) - 3./32.*(Macro[idm2y][2]+Macro[idp2y][2]) + 15./64.*(Macro[idm1y][2] + Macro[idp1y][2]) + 11./16. * Macro[id][2];
        Macro_filt[3] = 1./64.*(Macro[idm3y][3] +  Macro[idp3y][3]) - 3./32.*(Macro[idm2y][3]+Macro[idp2y][3]) + 15./64.*(Macro[idm1y][3] + Macro[idp1y][3]) + 11./16. * Macro[id][3];
        Macro_filt[4] = 1./64.*(Macro[idm3y][4] +  Macro[idp3y][4]) - 3./32.*(Macro[idm2y][4]+Macro[idp2y][4]) + 15./64.*(Macro[idm1y][4] + Macro[idp1y][4]) + 11./16. * Macro[id][4];
        Macro_filt[5] = 1./64.*(Macro[idm3y][5] +  Macro[idp3y][5]) - 3./32.*(Macro[idm2y][5]+Macro[idp2y][5]) + 15./64.*(Macro[idm1y][5] + Macro[idp1y][5]) + 11./16. * Macro[id][5];
        Macro_filt[6] = 1./64.*(Macro[idm3y][6] +  Macro[idp3y][6]) - 3./32.*(Macro[idm2y][6]+Macro[idp2y][6]) + 15./64.*(Macro[idm1y][6] + Macro[idp1y][6]) + 11./16. * Macro[id][6];

        Macro[id][0]=alpha * Macro_filt[0] + (1-alpha)*Macro[id][0];
        Macro[id][1]=alpha * Macro_filt[1] + (1-alpha)*Macro[id][1];
        Macro[id][2]=alpha * Macro_filt[2] + (1-alpha)*Macro[id][2];
        Macro[id][3]=alpha * Macro_filt[3] + (1-alpha)*Macro[id][3];
        Macro[id][4]=alpha * Macro_filt[4] + (1-alpha)*Macro[id][4];
        Macro[id][5]=alpha * Macro_filt[5] + (1-alpha)*Macro[id][5];
        Macro[id][6]=alpha * Macro_filt[6] + (1-alpha)*Macro[id][6];
        }
}

void FilterX7_LB(double **Macro, int which_slot){
  int i, j, id;
  int idm3x, idm2x, idm1x, idp1x, idp2x, idp3x;
  int off;

  double Macro_filt[7];

  off = which_slot*xMaxp_LB*yMaxp_LB;

  for(i=1; i<=xMax_LB; i++)
    for(j=3; j<(yMax_LB-1); j++){

        id = off + IDM_LB(i,j);

        idm3x = off + IDM_LB(i-3,j);
        idm2x = off + IDM_LB(i-2,j);
        idm1x = off + IDM_LB(i-1,j);
        idp1x = off + IDM_LB(i+1,j);
        idp2x = off + IDM_LB(i+2,j);
        idp3x = off + IDM_LB(i+3,j);

        if(i==1){
            idm3x = off + IDM_LB(xMax_LB-2,j);
            idm2x = off + IDM_LB(xMax_LB-1,j);
            idm1x = off + IDM_LB(xMax_LB,j);
        }
        else if(i==2){
            idm3x = off + IDM_LB(xMax_LB-1,j);
            idm2x = off + IDM_LB(xMax_LB,j);
        }
        else if(i==3){
            idm3x = off + IDM_LB(xMax_LB,j);
        }
        else if(i==xMax_LB){
            idp1x = off + IDM_LB(1,j);
            idp2x = off + IDM_LB(2,j);
            idp3x = off + IDM_LB(3,j);
        }
        else if(i==(xMax_LB-1)){
            idp2x = off + IDM_LB(1,j);
            idp3x = off + IDM_LB(2,j);
        }
        else if(i==(xMax_LB-2)){
            idp3x = off + IDM_LB(1,j);
        }

        if(j==1 || j==yMax_LB){
        Macro_filt[0] = 1./4.*(Macro[idm1x][0] + Macro[idp1x][0]) + 2./4.*Macro[id][0];
        Macro_filt[1] = 1./4.*(Macro[idm1x][1] + Macro[idp1x][1]) + 2./4.*Macro[id][1];
        Macro_filt[2] = 1./4.*(Macro[idm1x][2] + Macro[idp1x][2]) + 2./4.*Macro[id][2];
        Macro_filt[3] = 1./4.*(Macro[idm1x][3] + Macro[idp1x][3]) + 2./4.*Macro[id][3];
        Macro_filt[4] = 1./4.*(Macro[idm1x][4] + Macro[idp1x][4]) + 2./4.*Macro[id][4];
        Macro_filt[5] = 1./4.*(Macro[idm1x][5] + Macro[idp1x][5]) + 2./4.*Macro[id][5];
        Macro_filt[6] = 1./4.*(Macro[idm1x][6] + Macro[idp1x][6]) + 2./4.*Macro[id][6];
        }

        else if(j==2 || j==(yMax_LB-1)){
        Macro_filt[0] = -1./16.*( Macro[idm2x][0] + Macro[idp2x][0]) + 2./8.*(Macro[idm1x][0] + Macro[idp1x][0]) + 5./8.*Macro[id][0];
        Macro_filt[1] = -1./16.*( Macro[idm2x][1] + Macro[idp2x][1]) + 2./8.*(Macro[idm1x][1] + Macro[idp1x][1]) + 5./8.*Macro[id][1];
        Macro_filt[2] = -1./16.*( Macro[idm2x][2] + Macro[idp2x][2]) + 2./8.*(Macro[idm1x][2] + Macro[idp1x][2]) + 5./8.*Macro[id][2];
        Macro_filt[3] = -1./16.*( Macro[idm2x][3] + Macro[idp2x][3]) + 2./8.*(Macro[idm1x][3] + Macro[idp1x][3]) + 5./8.*Macro[id][3];
        Macro_filt[4] = -1./16.*( Macro[idm2x][4] + Macro[idp2x][4]) + 2./8.*(Macro[idm1x][4] + Macro[idp1x][4]) + 5./8.*Macro[id][4];
        Macro_filt[5] = -1./16.*( Macro[idm2x][5] + Macro[idp2x][5]) + 2./8.*(Macro[idm1x][5] + Macro[idp1x][5]) + 5./8.*Macro[id][5];
        Macro_filt[6] = -1./16.*( Macro[idm2x][6] + Macro[idp2x][6]) + 2./8.*(Macro[idm1x][6] + Macro[idp1x][6]) + 5./8.*Macro[id][6];
        }

        else{
        Macro_filt[0] = 1./64.*(Macro[idm3x][0] +  Macro[idp3x][0]) - 3./32.*(Macro[idm2x][0]+Macro[idp2x][0]) + 15./64.*(Macro[idm1x][0]+Macro[idp1x][0]) + 11./16.*Macro[id][0];
        Macro_filt[1] = 1./64.*(Macro[idm3x][1] +  Macro[idp3x][1]) - 3./32.*(Macro[idm2x][1]+Macro[idp2x][1]) + 15./64.*(Macro[idm1x][1]+Macro[idp1x][1]) + 11./16.*Macro[id][1];
        Macro_filt[2] = 1./64.*(Macro[idm3x][2] +  Macro[idp3x][2]) - 3./32.*(Macro[idm2x][2]+Macro[idp2x][2]) + 15./64.*(Macro[idm1x][2]+Macro[idp1x][2]) + 11./16.*Macro[id][2];
        Macro_filt[3] = 1./64.*(Macro[idm3x][3] +  Macro[idp3x][3]) - 3./32.*(Macro[idm2x][3]+Macro[idp2x][3]) + 15./64.*(Macro[idm1x][3]+Macro[idp1x][3]) + 11./16.*Macro[id][3];
        Macro_filt[4] = 1./64.*(Macro[idm3x][4] +  Macro[idp3x][4]) - 3./32.*(Macro[idm2x][4]+Macro[idp2x][4]) + 15./64.*(Macro[idm1x][4]+Macro[idp1x][4]) + 11./16.*Macro[id][4];
        Macro_filt[5] = 1./64.*(Macro[idm3x][5] +  Macro[idp3x][5]) - 3./32.*(Macro[idm2x][5]+Macro[idp2x][5]) + 15./64.*(Macro[idm1x][5]+Macro[idp1x][5]) + 11./16.*Macro[id][5];
        Macro_filt[6] = 1./64.*(Macro[idm3x][6] +  Macro[idp3x][6]) - 3./32.*(Macro[idm2x][6]+Macro[idp2x][6]) + 15./64.*(Macro[idm1x][6]+Macro[idp1x][6]) + 11./16.*Macro[id][6];
        }

        Macro[id][0]=alpha * Macro_filt[0] + (1-alpha)*Macro[id][0];
        Macro[id][1]=alpha * Macro_filt[1] + (1-alpha)*Macro[id][1];
        Macro[id][2]=alpha * Macro_filt[2] + (1-alpha)*Macro[id][2];
        Macro[id][3]=alpha * Macro_filt[3] + (1-alpha)*Macro[id][3];
        Macro[id][4]=alpha * Macro_filt[4] + (1-alpha)*Macro[id][4];
        Macro[id][5]=alpha * Macro_filt[5] + (1-alpha)*Macro[id][5];
        Macro[id][6]=alpha * Macro_filt[6] + (1-alpha)*Macro[id][6];
        }
}


void FilterY7_LB(double **Macro, int which_slot){
  int i, j, id;
  int idm3y, idm2y, idm1y, idp1y, idp2y, idp3y;
  int off;

  double Macro_filt[7];


  off = which_slot*xMaxp_LB*yMaxp_LB;

  for(i=1; i<=xMax_LB; i++)
    for(j=3; j<(yMax_LB-1); j++){

        id = off + IDM_LB(i,j);

        idm3y = off + IDM_LB(i,j-3);
        idm2y = off + IDM_LB(i,j-2);
        idm1y = off + IDM_LB(i,j-1);
        idp1y = off + IDM_LB(i,j+1);
        idp2y = off + IDM_LB(i,j+2);
        idp3y = off + IDM_LB(i,j+3);


        if(j==1 || j==yMax_LB){
        Macro_filt[0] = 1./4.*(Macro[idm1y][0] + Macro[idp1y][0]) + 2./4. * Macro[id][0];
        Macro_filt[1] = 1./4.*(Macro[idm1y][1] + Macro[idp1y][1]) + 2./4. * Macro[id][1];
        Macro_filt[2] = 1./4.*(Macro[idm1y][2] + Macro[idp1y][2]) + 2./4. * Macro[id][2];
        Macro_filt[3] = 1./4.*(Macro[idm1y][3] + Macro[idp1y][3]) + 2./4. * Macro[id][3];
        Macro_filt[4] = 1./4.*(Macro[idm1y][4] + Macro[idp1y][4]) + 2./4. * Macro[id][4];
        Macro_filt[5] = 1./4.*(Macro[idm1y][5] + Macro[idp1y][5]) + 2./4. * Macro[id][5];
        Macro_filt[6] = 1./4.*(Macro[idm1y][6] + Macro[idp1y][6]) + 2./4. * Macro[id][6];

        }
        else if(j==2 || j==(yMax_LB-1)){
        Macro_filt[0] = -1./16.*( Macro[idm2y][0] + Macro[idp2y][0]) + 2./8.*(Macro[idm1y][0] + Macro[idp1y][0]) + 5./8.*Macro[id][0];
        Macro_filt[1] = -1./16.*( Macro[idm2y][1] + Macro[idp2y][1]) + 2./8.*(Macro[idm1y][1] + Macro[idp1y][1]) + 5./8.*Macro[id][1];
        Macro_filt[2] = -1./16.*( Macro[idm2y][2] + Macro[idp2y][2]) + 2./8.*(Macro[idm1y][2] + Macro[idp1y][2]) + 5./8.*Macro[id][2];
        Macro_filt[3] = -1./16.*( Macro[idm2y][3] + Macro[idp2y][3]) + 2./8.*(Macro[idm1y][3] + Macro[idp1y][3]) + 5./8.*Macro[id][3];
        Macro_filt[4] = -1./16.*( Macro[idm2y][4] + Macro[idp2y][4]) + 2./8.*(Macro[idm1y][4] + Macro[idp1y][4]) + 5./8.*Macro[id][4];
        Macro_filt[5] = -1./16.*( Macro[idm2y][5] + Macro[idp2y][5]) + 2./8.*(Macro[idm1y][5] + Macro[idp1y][5]) + 5./8.*Macro[id][5];
        Macro_filt[6] = -1./16.*( Macro[idm2y][6] + Macro[idp2y][6]) + 2./8.*(Macro[idm1y][6] + Macro[idp1y][6]) + 5./8.*Macro[id][6];
        }
        else{
        Macro_filt[0] = 1./64.*(Macro[idm3y][0] +  Macro[idp3y][0]) - 3./32.*(Macro[idm2y][0]+Macro[idp2y][0]) + 15./64.*(Macro[idm1y][0] + Macro[idp1y][0]) + 11./16. * Macro[id][0];
        Macro_filt[1] = 1./64.*(Macro[idm3y][1] +  Macro[idp3y][1]) - 3./32.*(Macro[idm2y][1]+Macro[idp2y][1]) + 15./64.*(Macro[idm1y][1] + Macro[idp1y][1]) + 11./16. * Macro[id][1];
        Macro_filt[2] = 1./64.*(Macro[idm3y][2] +  Macro[idp3y][2]) - 3./32.*(Macro[idm2y][2]+Macro[idp2y][2]) + 15./64.*(Macro[idm1y][2] + Macro[idp1y][2]) + 11./16. * Macro[id][2];
        Macro_filt[3] = 1./64.*(Macro[idm3y][3] +  Macro[idp3y][3]) - 3./32.*(Macro[idm2y][3]+Macro[idp2y][3]) + 15./64.*(Macro[idm1y][3] + Macro[idp1y][3]) + 11./16. * Macro[id][3];
        Macro_filt[4] = 1./64.*(Macro[idm3y][4] +  Macro[idp3y][4]) - 3./32.*(Macro[idm2y][4]+Macro[idp2y][4]) + 15./64.*(Macro[idm1y][4] + Macro[idp1y][4]) + 11./16. * Macro[id][4];
        Macro_filt[5] = 1./64.*(Macro[idm3y][5] +  Macro[idp3y][5]) - 3./32.*(Macro[idm2y][5]+Macro[idp2y][5]) + 15./64.*(Macro[idm1y][5] + Macro[idp1y][5]) + 11./16. * Macro[id][5];
        Macro_filt[6] = 1./64.*(Macro[idm3y][6] +  Macro[idp3y][6]) - 3./32.*(Macro[idm2y][6]+Macro[idp2y][6]) + 15./64.*(Macro[idm1y][6] + Macro[idp1y][6]) + 11./16. * Macro[id][6];
        }

        Macro[id][0]=alpha * Macro_filt[0] + (1-alpha)*Macro[id][0];
        Macro[id][1]=alpha * Macro_filt[1] + (1-alpha)*Macro[id][1];
        Macro[id][2]=alpha * Macro_filt[2] + (1-alpha)*Macro[id][2];
        Macro[id][3]=alpha * Macro_filt[3] + (1-alpha)*Macro[id][3];
        Macro[id][4]=alpha * Macro_filt[4] + (1-alpha)*Macro[id][4];
        Macro[id][5]=alpha * Macro_filt[5] + (1-alpha)*Macro[id][5];
        Macro[id][6]=alpha * Macro_filt[6] + (1-alpha)*Macro[id][6];
        }
}

void ComputeMassNS(double **Macro_p){
	int i, j, id, off;
	off = current_slot * xMaxp_NS * yMaxp_NS;

	double MASS = 0.0;

	for(i=1; i<=xMax_NS; i++)
		for(j=1; j<=yMax_NS; j++){
            id = off + IDM_NS(i,j);
            MASS += Macro_p[id][0];
				}
	fprintf(stdout, "mass NS : %9.12f\n", MASS);
}

void ComputeMassLB(double *cell){
	int i, j, l, id, off;
	off = current_slot * xMaxp_LB * yMaxp_LB;

	double MASS = 0.0;

	for(i=1; i<=xMax_LB; i++)
		for(j=1; j<=yMax_LB; j++)
            for(l=0; l<NPOP; l++){
            id = off * NPOP + IDF(i,j,l);
            MASS += cell[id];
				}
	fprintf(stdout, "mass LB : %9.12f\n", MASS);
}
