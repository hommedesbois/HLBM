 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include "ludwig.h"


void Streaming(double *cell){
	int i, j, l, inv, idx0, idx1, off0, off1;

	off0 = current_slot * NPOP * yMaxp_LB * xMaxp_LB;
	off1 = other_slot * NPOP * yMaxp_LB * xMaxp_LB;

	for(i=1; i<=xMax_LB; i++)
		for(j=1; j<=yMax_LB; j++)

				for(l=0; l<NPOP; l++){
					inv = finv[l];

					idx0 = off0 + IDF(i+ex[inv],j+ey[inv],l);
					idx1 = off1 + IDF(i,j,l);

					cell[idx1] = cell[idx0];
				}
}
