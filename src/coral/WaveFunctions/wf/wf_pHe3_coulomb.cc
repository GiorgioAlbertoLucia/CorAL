// <<BEGIN-copyright>>
// 
//                 The GNU General Public License (GPL) Version 2, June 1991
// 
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence 
// Livermore National Laboratory. Written by Ron Soltz (soltz1@llnl.gov), David A. Brown 
// (dbrown@bnl.gov) and Scott Pratt (pratts@pa.msu.edu).
// 
// CODE-CODE-643336 All rights reserved. 
// 
// This file is part of CorAL, Version: 1.17.
// 
// Please see the file LICENSE.TXT in the main directory of this source code distribution.
// 
// This program is free software; you can redistribute it and/or modify it under the terms of 
// the GNU General Public License (as published by the Free Software Foundation) version 2, 
// dated June 1991.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the terms and conditions of the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with this program; 
// if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
// MA 02111-1307 USA
// 
// <<END-copyright>>
#ifndef __INCLUDE_PHE3_COULOMB_CC__
#define __INCLUDE_PHE3_COULOMB_CC__
#include "wavefunction.h"
using namespace std;

CWaveFunction_pHe3_Coulomb::CWaveFunction_pHe3_Coulomb(string parsfilename) {
    
    m1 = MPROTON;
    m2 = MHE3;
    ParsInit(parsfilename);
    q1q2 = (COULOMB == 0) ? 0. : 2.; // Positive charge interaction for proton-He3
    nchannels = 2;
    
    InitArrays(); // Initializes arrays if needed, inherited from CWaveFunction
    InitWaves();  // Assumes this sets up plane waves or other necessary initializations

    //nchannels = 0;
    channelweight[0] = 0.25;
    channelweight[1] = 0.75;
    //for(int ichannel = 0; ichannel < nchannels; ichannel++){
    //    for(int iq = 0; iq < nqmax; iq++){
    //        double q = qarray[iq];
    //        Wepsilon[ichannel][iq] = ddeltadq[ichannel][iq] 
    //            - GetIW(ell[ichannel], epsilon, q, q1q2, eta[iq], delta[ichannel][iq])
    //	        + GetIW(ell[ichannel], epsilon, q, q1q2, eta[iq], 0.0);
    //        Wepsilon[ichannel][iq] = 3.0 * Wepsilon[ichannel][iq] * pow(HBARC/epsilon,3)/(2.0*q*q);
    //
    //    }
    //}

    printf("Proton-Helium3 Coulomb wf initialized\n");
}

double CWaveFunction_pHe3_Coulomb::CalcPsiSquared(int iq, double r, double ctheta) {

    double psisquared = 0.0;
    double q;
    const double ROOT2 = sqrt(2.0);
    complex<double> psi, psia, psib;
    
    if (iq >= nqmax) {
        psisquared = 1.0;
    } else {
        q = qarray[iq];
        psia = planewave[iq]->planewave(r, ctheta);
        psib = planewave[iq]->planewave(r, -ctheta);
        psi = (psia + psib) / ROOT2;
        psisquared = real(psi * conj(psi));
        psisquared = real(psia * conj(psia));
    }

    // psisquared *= RelativisticCorrection(r, iq); // check wht this does
    return psisquared;
   
    /*
    double psisquared;
	complex<double> psi0;

	if(iq>=nqmax){
		psisquared=1.0;
	}
	else{
		psi0=planewave[iq]->planewave(r,ctheta);
		psisquared=real(psi0*conj(psi0));
		//psisquared*=RelativisticCorrection(r,iq);
	}
	return psisquared;
    */
}

/*
CWaveFunction_pHe3_Coulomb::~CWaveFunction_pHe3_Coulomb() {
    // Cleanup if necessary
}
*/

#endif
