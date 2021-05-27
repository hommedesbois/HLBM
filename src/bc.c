 /* T. Horstmann
    AT-TRA
    Deutsches Zentrum für Luft- und Raumfahrt (DLR)
    Berlin
	Germany
 */

#include <stdio.h>
#include <math.h>
#include "ludwig.h"


void ComputeFcol_BC(double *cell){
	int j,l;

	/*
	 *		Y
	 *		^
	 *		|
	 *		|
	 *		Z ----> X
	 */

	int i_in, i_out;
	int idx_in, idx_out, offset;

	offset = current_slot*NPOP*xMaxp_LB*yMaxp_LB;

	// faces
	i_out = 1;
	i_in = xMax_LB + 1;

	for(j=0; j<yMaxp_LB; j++)
			for(l=0; l<NPOP; l++){
				idx_out =  offset + IDF(i_out,j,l);
				idx_in = offset + IDF(i_in,j,l);
				cell[idx_in] = cell[idx_out];
			}

	i_out = xMax_LB;
	i_in = 0;

	for(j=0; j<yMaxp_LB; j++)
			for(l=0; l<NPOP; l++){
				idx_out =  offset + IDF(i_out,j,l);
				idx_in = offset + IDF(i_in,j,l);
				cell[idx_in] = cell[idx_out];
			}
}

void ComputeMacroLB_BC(double **Macro){
	int j,k;

	int i_in, i_out;
	int idx_in, idx_out, offset;

	offset = current_slot*xMaxp_LB*yMaxp_LB;

	// faces
	i_out = 1;
	i_in = xMax_LB + 1;

	for(j=0; j<=NY_LB; j++){			// change to 0 and <yMaxp when periodicity only in one direction
        idx_out =  offset + IDM_LB(i_out, j);
        idx_in = offset + IDM_LB(i_in, j);
        for(k=0; k<=2; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
            }
        }

	i_out = xMax_LB;
	i_in = 0;

	for(j=0; j<=NY_LB; j++){
        idx_out =  offset + IDM_LB(i_out, j);
        idx_in = offset + IDM_LB(i_in, j);
        for(k=0; k<=2; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
			}
    }
}

void ComputeMacro_BC(double **Macro){
	int i,j,k;

	int i_in, i_out, j_in, j_out;
	int idx_in, idx_out, offset;

	offset = other_slot*xMaxp_NS*yMaxp_NS;

	// faces
	i_out = 1;
	i_in = xMax_NS + 1;

	for(j=1; j<=yMax_NS; j++){			// change to 0 and <yMaxp when periodicity only in one direction
        idx_out =  offset + IDM_NS(i_out, j);
        idx_in = offset + IDM_NS(i_in, j);
        for(k=0; k<=2; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
            }
        }

	i_out = xMax_NS;
	i_in = 0;

	for(j=1; j<=yMax_NS; j++){
        idx_out =  offset + IDM_NS(i_out, j);
        idx_in = offset + IDM_NS(i_in, j);
        for(k=0; k<=2; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
			}
        }
	j_out = 1;
	j_in = yMax_NS + 1;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  offset + IDM_NS(i, j_out);
	idx_in = offset + IDM_NS(i, j_in);
	for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
	}

	j_out = yMax_NS;
	j_in = 0;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  offset + IDM_NS(i, j_out);
	idx_in = offset + IDM_NS(i, j_in);
	for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
	}

	//edges
    	i_out = 1;
    	j_out = 1;
    	i_in = xMax_NS + 1;
    	j_in = yMax_NS + 1;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}


    	i_out = xMax_NS;
    	j_out = yMax_NS;
    	i_in = 0;
    	j_in = 0;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}

    	i_out = xMax_NS;
    	j_out = 1;
    	i_in  = 0;
    	j_in  = yMax_NS + 1;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);

        for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
            }

    	i_out = 1;
    	j_out = yMax_NS;
    	i_in = xMax_NS + 1;
    	j_in = 0;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
}

void ComputeMacroP_BC(double **Macro){
	int i,j,k;

	/*
	 *		Y
	 *		^
	 *		|
	 *		|
	 *		Z ----> X
	 */

	int i_in, i_out, j_in, j_out;
	int idx_in, idx_out, offset;

	offset = other_slot*xMaxp_NS*yMaxp_NS;

	// faces
	i_out = 1;
	i_in = xMax_NS + 1;

	for(j=1; j<=yMax_NS; j++){  			// change to 0 and <yMaxp when periodicity only in one direction
        idx_out =  offset + IDM_NS(i_out, j);
        idx_in = offset + IDM_NS(i_in, j);
        for(k=3; k<=6; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
            }
        }

	i_out = xMax_NS;
	i_in = 0;

	for(j=1; j<=yMax_NS; j++){
        idx_out =  offset + IDM_NS(i_out, j);
        idx_in = offset + IDM_NS(i_in, j);
        for(k=3; k<=6; k++){
            Macro[idx_in][k] = Macro[idx_out][k];
			}
	}

	j_out = 1;
	j_in = yMax_NS + 1;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  offset + IDM_NS(i,j_out);
	idx_in = offset + IDM_NS(i,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
	}

	j_out = yMax_NS;
	j_in = 0;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  offset + IDM_NS(i,j_out);
	idx_in = offset + IDM_NS(i,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
	}

	//edges
    	i_out = 1;
    	j_out = 1;
    	i_in = xMax_NS + 1;
    	j_in = yMax_NS + 1;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}


    	i_out = xMax_NS;
    	j_out = yMax_NS;
    	i_in = 0;
    	j_in = 0;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}

    	i_out = xMax_NS;
    	j_out = 1;
    	i_in = 0;
    	j_in = yMax_NS + 1;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}

    	i_out = 1;
    	j_out = yMax_NS;
    	i_in = xMax_NS + 1;
    	j_in = 0;

    	idx_out = offset + IDM_NS(i_out, j_out);
    	idx_in = offset + IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            Macro[idx_in][k] = Macro[idx_out][k];
		}
}


void GuardRK(double ***RK, int step){
	int i,j,k;

	int i_in, i_out, j_in, j_out;
	int idx_in, idx_out;

	// faces
	i_out = 1;
	i_in = xMax_NS + 1;

	for(j=1; j<=yMax_NS; j++){   // change to 0 and <yMaxp when periodicity only in one direction
        idx_out =  IDM_NS(i_out, j);
        idx_in = IDM_NS(i_in, j);
        for(k=0; k<=2; k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
            }
        }

	i_out = xMax_NS;
	i_in = 0;

	for(j=1; j<=yMax_NS; j++){
        idx_out = IDM_NS(i_out, j);
        idx_in = IDM_NS(i_in, j);
        for(k=0; k<=2; k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
        }
    }
	j_out = 1;
	j_in = yMax_NS + 1;

	for(i=1; i<=xMax_NS; i++){
        idx_out = IDM_NS(i,j_out);
	idx_in = IDM_NS(i,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}
	}

	j_out = yMax_NS;
	j_in = 0;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  IDM_NS(i,j_out);
	idx_in = IDM_NS(i,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}
	}

	//edges
    	i_out = 1;
    	j_out = 1;
    	i_in = xMax_NS + 1;
    	j_in = yMax_NS + 1;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}


    	i_out = xMax_NS;
    	j_out = yMax_NS;
    	i_in = 0;
    	j_in = 0;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

    	i_out = xMax_NS;
    	j_out = 1;
    	i_in = 0;
    	j_in = yMax_NS + 1;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

    	i_out = 1;
    	j_out = yMax_NS;
    	i_in = xMax_NS + 1;
    	j_in = 0;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=0; k<=2;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

}

void GuardRKP(double ***RK, int step){
	int i,j,k;

	/*
	 *		Y
	 *		^
	 *		|
	 *		|
	 *		Z ----> X
	 */

	int i_in, i_out, j_in, j_out;
	int idx_in, idx_out;

	// faces
	i_out = 1;
	i_in = xMax_NS + 1;

	for(j=1; j<=yMax_NS; j++){ 	// change to 0 and <yMaxp when periodicity only in one direction
        idx_out = IDM_NS(i_out, j);
        idx_in = IDM_NS(i_in, j);
        for(k=3; k<=6; k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
            }
        }

	i_out = xMax_NS;
	i_in = 0;

	for(j=1; j<=yMax_NS; j++){
        idx_out = IDM_NS(i_out, j);
        idx_in = IDM_NS(i_in, j);
        for(k=3; k<=6; k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
			}
        }

	j_out = 1;
	j_in = yMax_NS + 1;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  IDM_NS(i,j_out);
	idx_in = IDM_NS(i,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}
	}

	j_out = yMax_NS;
	j_in = 0;

	for(i=1; i<=xMax_NS; i++){
        idx_out =  IDM_NS(i,j_out);
	idx_in = IDM_NS(i,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}
	}

	//edges
    	i_out = 1;
    	j_out = 1;
    	i_in = xMax_NS + 1;
    	j_in = yMax_NS + 1;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}


    	i_out = xMax_NS;
    	j_out = yMax_NS;
    	i_in = 0;
    	j_in = 0;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

    	i_out = xMax_NS;
    	j_out = 1;
    	i_in = 0;
    	j_in = yMax_NS + 1;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

    	i_out = 1;
    	j_out = yMax_NS;
    	i_in = xMax_NS + 1;
    	j_in = 0;

    	idx_out = IDM_NS(i_out, j_out);
    	idx_in = IDM_NS(i_in,j_in);
	for(k=3; k<=6;k++){
            RK[step][idx_in][k] = RK[step][idx_out][k];
		}

}