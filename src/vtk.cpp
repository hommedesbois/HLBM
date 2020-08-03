#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

using namespace std;
#include "ludwig.h"

#ifndef __cplusplus
#error C++ compiler is required
#endif

extern "C" void DumpMacroNS(double **macro, int index){

	 int off = current_slot * xMaxp_NS * yMaxp_NS;

	stringstream output_filename;
	output_filename << "vtk/FluidSimulationDomain_t" << index << ".vtk";
	ofstream output_file;

	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS "<< NX_NS << " " << NY_NS << " " << 1 << endl;

	output_file << "X_COORDINATES " << NX_NS << " float\n";
	for(int i=0; i< NX_NS; i++)
		output_file << (double) i*DX << " ";
	output_file << endl;

	output_file << "Y_COORDINATES " << NY_NS << " float\n";
	for(int j=0; j< NY_NS; j++)
		output_file << (double) j*DX << " ";
	output_file << endl;

    output_file << "Z_COORDINATES " << 1 << " float\n";
    output_file << 0 << endl;


	output_file << "CELL_DATA " << xMax_NS * yMax_NS << endl;
	output_file << "SCALARS density float 1 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=yMax_NS; Y++)
			for(int X=1; X<= xMax_NS; X++){

				int idx = off + IDM_NS(X,Y);
				output_file << (macro[idx][0] - rho0) << endl;
			}
    output_file << "SCALARS velocity float 3 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=yMax_NS; Y++)
			for(int X=1; X<= xMax_NS; X++){

				int idx = off + IDM_NS(X,Y);
				double ux = (macro[idx][1]/macro[idx][0])*Csound*sqrt(3);
				double uy = (macro[idx][2]/macro[idx][0])*Csound*sqrt(3);
				double uz = 0;
				output_file << ux << " " << uy << " " << uz << endl;
			}

	output_file.close();
}

extern "C" void DumpMacroLB(double **macro, int index){

    int off = current_slot * xMaxp_LB * yMaxp_LB;

	stringstream output_filename;
	output_filename << "vtk/lbm_t" << index << ".vtk";
	ofstream output_file;


	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS "<< NX_LB << " " << NY_LB << " " << 1 << endl;

	output_file << "X_COORDINATES " << NX_LB << " float\n";
	for(int i=0; i<NX_LB; i++)
		output_file << (double) i*DX << " ";
	output_file << endl;

	output_file << "Y_COORDINATES " << NY_LB << " float\n";
	for(int j=0; j<NY_LB; j++)
		output_file << (double) (j+shift)*DX << " ";
	output_file << endl;

    output_file << "Z_COORDINATES " << 1 << " float\n";
    output_file << 0.01 << endl;


	output_file << "CELL_DATA " << xMax_LB * yMax_LB << endl;
	output_file << "SCALARS density float 1 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=yMax_LB; Y++)
			for(int X=1; X<= xMax_LB; X++){

				int idx = off + IDM_LB(X,Y);
				output_file << (macro[idx][0] - rho0) << endl;
			}
    output_file << "SCALARS velocity float 3 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=yMax_LB; Y++)
			for(int X=1; X<= xMax_LB; X++){

				int idx = off + IDM_LB(X,Y);
				double ux = (macro[idx][1]/macro[idx][0])*Csound*sqrt(3);
				double uy = (macro[idx][2]/macro[idx][0])*Csound*sqrt(3);
				double uz = 0.0;
				output_file << ux << " " << uy << " " << uz << endl;
			}

	output_file.close();
}
