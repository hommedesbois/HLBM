/* T. Horstmann
* Department of Aeroacoustics
* ONERA
* Chatillon
* France
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ludwig.h"

void ComputeSigma(double *sigma, int it){
int i,j;
double s;
double w = 0.4;
int Nw_0 = w/DX + 1;
int Nw_max = xMax_NS - Nw_0 + 1;

//fprintf(stdout, "Nw_0 %d, Nw_max %d, width %g, DX %g\n", Nw_0, Nw_max, w);


for(i=1; i<=xMax_NS; i++)
    for(j=1; j<=xMax_NS; j++){

        int id = (j-1)*xMax_NS + (i-1);

        if(j<=Nw_0)
            s = 3125./256. * ((1-j)*DX * (((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX))/(-w*-w*-w*-w*-w));
//          else if (j>=Nw_max && j<=xMax_NS)
//          s = 3125./256. * ((xMax_NS-j)*DX * (((j - Nw_max)*DX) * ((j - Nw_max)*DX) * ((j-Nw_max)*DX) * ((j-Nw_max)*DX)) / (w*w*w*w*w));
        else
            s = 0.0;

        sigma[id] = s;
    }

for(j=1; j<=xMax_NS; j++)
    for(i=j; i<=(xMax_NS+1-j); i++){

        int id = (i-1)*xMax_NS + (j-1);

        if(j<=Nw_0)
            s = 3125./256. * ((1-j)*DX * (((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX))/(-w*-w*-w*-w*-w));
        else if (j>=Nw_max && j<=xMax_NS)
            s = 3125./256. * ((xMax_NS-j)*DX * (((j - Nw_max)*DX) * ((j - Nw_max)*DX) * ((j-Nw_max)*DX) * ((j-Nw_max)*DX)) / (w*w*w*w*w));
        else
            s = 0.0;

        sigma[id] = s;
    }

for(j=1; j<=xMax_NS; j++)
    for(i=(xMax_NS-j+1); i<=(xMax_NS-(xMax_NS-j)); i++){

        int id = (i-1)*xMax_NS + (j-1);

        if(j<=Nw_0)
            s = 3125./256. * ((1-j)*DX * (((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX) * ((j - Nw_0)*DX))/(-w*-w*-w*-w*-w));
        else if (j>=Nw_max && j<=xMax_NS)
            s = 3125./256. * ((xMax_NS-j)*DX * (((j - Nw_max)*DX) * ((j - Nw_max)*DX) * ((j-Nw_max)*DX) * ((j-Nw_max)*DX)) / (w*w*w*w*w));
        else
            s = 0.0;

        sigma[id] = s;
    }

for(j=Nw_0; j<=yMax_NS; j++)
    for(i=1; i<=xMax_NS; i++){

        int id = (j-1)*xMax_NS + (i-1);

        if(i<=Nw_0)
            s = 3125./256. * ((1-i)*DX * (((i - Nw_0)*DX) * ((i - Nw_0)*DX) * ((i - Nw_0)*DX) * ((i - Nw_0)*DX))/(-w*-w*-w*-w*-w));
        else if (i>=Nw_max && i<=xMax_NS)
            s = 3125./256. * ((xMax_NS-i)*DX * (((i - Nw_max)*DX) * ((i - Nw_max)*DX) * ((i-Nw_max)*DX) * ((i-Nw_max)*DX)) / (w*w*w*w*w));
        else
            s = 0.0;

        sigma[id] = s;
        }
    }