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
#include "parametermap.h"
#include "cheezyparser.h"
#include "random.h"
#include "oscar.h"
#include "oscar_source_generator_1d.h"
#include "oscar_correlation_generator_1d.h"
#include "oscar_source_generator_3dcart.h"
#include "oscar_source_generator_3dsphr.h"
#include "oscar_correlation_generator_3dcart.h"
#include "message.h"
#include "basisfunc_imager1d.h"
#include "expander.h"
#include "uncoupled_imager3d.h"
#include "corr3d_ylm.h"
#include "tnt_array1d.h"
#include "wavefunction.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <memory>

struct OutputHistograms {
    CHistogram3d* source;
    CHistogram1d* CF;
    CHistogram1d* kstarNum;
    CHistogram1d* kstarDen;
    CHistogram1d* kstarX;
    CHistogram1d* kstarY;
    CHistogram1d* kstarZ;

    OutputHistograms() {
        source = new CHistogram3d();
        CF = new CHistogram1d();
        kstarNum = new CHistogram1d();
        kstarDen = new CHistogram1d();
        kstarX = new CHistogram1d();
        kstarY = new CHistogram1d();
        kstarZ = new CHistogram1d();
    }

    ~OutputHistograms() {
        std::cout << "Deleting output histograms" << std::endl;
        delete source;
        delete CF;
        delete kstarNum;
        delete kstarDen;
        delete kstarX;
        delete kstarY;
        delete kstarZ;
    }
};

enum class optionCF {
    oneD,
    oneDParallel,
    threeD
};

//-------------------------------------------------------------------------
// Function declarations
//-------------------------------------------------------------------------

void getHelp(void);
void computeCorrelation1d( parameterMap& inMap );
void displaySource( parameterMap& inMap );
CHistogram3d generateGaussianSource( parameterMap m );
void generateEvents1d( OutputHistograms& output1d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf );
void generateEvents1dParallel( OutputHistograms& output1d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf );
void generateEvents3d( OutputHistograms& output3d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf );

//-------------------------------------------------------------------------
// Main routines
//-------------------------------------------------------------------------

//! Print usage information & quit
void getHelp(void){
    cout<<"\nUsage: phemto <mode> [inputFile.dat]"<<endl;
    cout<<"    Valid modes are: --1d            : Generate 1d source and correlation"<<endl;
    cout<<"                     --source        : Generate 3d source distribution"<<endl;
    cout<<"                     -h, --help      : Prints this message, then quits"<<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** 1d implementation of the Proton HElium3 feMTOscopy (PHEMTO) using CorAL ***"<<endl;
    cout << "Version: 1.0"<<endl;
    cout << "Units: " << endl;
    cout << "    Length: fm"<<endl;
    cout << "    Energy, Mass, Momentum: MeV"<<endl;
    cout << endl;

    MESSAGE << CMessage::warning;

    bool skip_correlation=false;
    bool make_plots=false;

    // Parse command line
    if (argc==1) getHelp();
    string paramFile("");
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="--help") getHelp();
        if (sarg=="-h") getHelp();
        if (sarg=="--source") {
            cout << "Will output a 3d gaussian source distribution\n"<<endl;
        }
        else if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
        else paramFile = sarg;
    }

    // Read in the input parameters from the argument on the command line
    if (paramFile=="") {
        MESSAGE<<"No inputFile parameter file given!!"<<ENDM_FATAL;
        getHelp();
    }
    parameterMap inMap;
    ReadParsFromFile(inMap, paramFile);

    // Generate the correlation from the source, then put hte results in the map
    for (vector<string>::iterator mode=modeList.begin(); mode!=modeList.end(); ++mode) {
        if      (*mode == "--1d")           computeCorrelation1d( inMap );
        else if (*mode == "--source")       displaySource( inMap );
        else MESSAGE<< "Unknown mode: "<<*mode<<ENDM_WARN;
        cout << endl;
    }
    return 0;
}

//-------------------------------------------------------------------------
// Function definitions
//-------------------------------------------------------------------------

/**
 * Generate a Gaussian source and compute the correlation function for p-He3
*/
void computeCorrelation1d( parameterMap& inMap ){
    parameterMap sourceInMap = parameter::getMap(inMap, "source_settings");
    parameterMap kstarInMap = parameter::getMap(inMap, "kstar_settings");
    std::string wfMapFile = parameter::getS(inMap, "wf_map_file", "/Users/glucia/Projects/CorAL/examples/PHEMTO/input/wf_input.dat");

    // Generate the correlation function
    OutputHistograms qaHists;
    qaHists.kstarNum->Read(kstarInMap);
    qaHists.kstarDen->Read(kstarInMap);
    qaHists.CF->Read(kstarInMap);
    
    CWaveFunction_pHe3_Coulomb wf(wfMapFile);
    generateEvents1d(qaHists, sourceInMap, kstarInMap, wf);
    printf("Calculating correlation function\n");
    // Calculate CF by dividing num by den
    for (int ibin = 0; ibin < qaHists.kstarNum->ndata; ++ibin) {
        if (qaHists.kstarDen->data[ibin] > 0.0) {
            qaHists.CF->data[ibin] = qaHists.kstarNum->data[ibin] / qaHists.kstarDen->data[ibin];
            qaHists.CF->uncert[ibin] = sqrt(qaHists.CF->data[ibin] * (1./qaHists.kstarNum->data[ibin] + 1./qaHists.kstarDen->data[ibin]));
        }
    }
    printf("Done calculating correlation function\n");


    //computeCF( optionCF::oneD, sourceInMap, kstarInMap, wfMapFile, qaHists );

    cout << "Saving correlation function to file" << endl;
    parameterMap CFOutMap;
    qaHists.CF->Write(CFOutMap);
    WriteParsToFile(CFOutMap, parameter::getS(inMap, "correlation_file", "output_correlation_1d.dat"));

    // Write the kstar distributions to file
    cout << "Saving kstar distributions to file" << endl;
    parameterMap kstarNumOutMap, kstarDenOutMap;
    qaHists.kstarNum->Write(kstarNumOutMap);
    WriteParsToFile(kstarNumOutMap, parameter::getS(inMap, "kstarNum_file", "output_kstarNum_1d.dat"));

}

/**
 * Generate a Gaussian source for p-He3
*/
void displaySource( parameterMap& inMap ){
    
    parameterMap sourceInMap = parameter::getMap(inMap, "source_settings");
    
    // Generate the source distribution (example)
    CHistogram3d source = generateGaussianSource( sourceInMap );
    parameterMap sourceOutMap;
    source.Write(sourceOutMap);
    WriteParsToFile(sourceOutMap, parameter::getS(inMap, "source_file", "output_source_3d.dat"));

}

/**
 * Generate a Gaussian source distribution (3d model)
 * 
 * @param m parameterMap containing the source parameters.
 * Requiered parameters are RX, RY, RZ, nMcSamples
*/
CHistogram3d generateGaussianSource( parameterMap m ){
    
    cout << "\nGenerating source" << endl;

    // Initialize the source histogram (nbins, bin_width, offset)
    CHistogram3d source;
    source.Read(m);

    // Get the source parameters (gaussian source with input radius)
    double RX = parameter::getD(m,"RX", 10.0);
    double RY = parameter::getD(m,"RY", 10.0);
    double RZ = parameter::getD(m,"RZ", 10.0);
    double root2 = sqrt(2.0);

    CRandom randy(-12345);

    int nMcSamples = parameter::getI(m,"nMcSamples", 1000000);
    for (int ipart=0; ipart<nMcSamples; ++ipart){
        
        double x = root2*RX*randy.gauss();
        double y = root2*RY*randy.gauss();
        double z = root2*RZ*randy.gauss();
        
        int ibin = source.whatBin(x,y,z);
        if (ibin>=0 && ibin<source.ndata) {
            source.data[ibin] += 1.0;
            source.uncert[ibin] += 1.0;
        }
    }

    return source;
}

/**
 * Simulate the correlation function in 1d by sampling a flat relative momentum distribution
 * and a Gaussian relative position distribution (3d).
 * 
 * @param output1d OutputHistograms containing the output histograms.
 * @param nMcSamples int number of Monte Carlo samples to take.
 * @param KSTAR_BINS int number of bins for the correlation function.   
 * @param KSTAR_MAX double maximum value of the correlation function.
 * @param DELKSTAR double bin width of the correlation function.
 * @param SIGNED_KSTAR bool whether to use signed k* values.
 * @param RX double source radius in x direction.
 * @param RY double source radius in y direction.
 * @param RZ double source radius in z direction.
 * 
*/
void generateEvents1d( OutputHistograms& output1d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf ) {
    
    // Get the source parameters
    double RX = parameter::getD(sourceMap, "RX", 10.0);
    double RY = parameter::getD(sourceMap, "RY", 10.0);
    double RZ = parameter::getD(sourceMap, "RZ", 10.0);
    const double root2_RX = sqrt(2.0) * RX;
    const double root2_RY = sqrt(2.0) * RY;
    const double root2_RZ = sqrt(2.0) * RZ;

    // Correlation function parameters
    bool SIGNED_KSTAR = parameter::getB(kstarMap, "SIGNED_KSTAR", false);
    double DELKSTAR = parameter::getD(kstarMap, "KSTAR_DELTA", 1.0);
    int KSTAR_BINS = parameter::getI(kstarMap, "KSTAR_BINS", 200);
    double KSTAR_MAX = parameter::getD(kstarMap, "KSTAR_MAX", 200.0);

    int nMcSamples = parameter::getI(sourceMap, "nMcSamples", 1000000);
    printf("Number of Monte Carlo samples: %d\n", nMcSamples);
    
    CRandom randy(-12345);  // Different seed per thread

    cout << "Simulating 1d correlation function" << endl;
    for (int ik = 0; ik < KSTAR_BINS; ++ik) {
        cout << "\r" << static_cast<int>(100.0 * ik / KSTAR_BINS) << "%" << flush;
        double k_base = DELKSTAR * ik + DELKSTAR / 2.;

        for (int ipart = 0; ipart < nMcSamples; ++ipart) {
            double k = k_base;

            double x = root2_RX * randy.gauss();
            double y = root2_RY * randy.gauss();
            double z = root2_RZ * randy.gauss();
            double r = sqrt(x * x + y * y + z * z);
            double ctheta = 0.; // WIP: angle between k and r

            double weight = wf.CalcPsiSquared(k, r, ctheta);
            int ibinKstar = output1d.kstarNum->whatBin(k);
            if (ibinKstar >= 0 && ibinKstar < output1d.kstarNum->ndata) {
                output1d.kstarNum->data[ibinKstar] += weight;    // WIP: weight is given by the wave function
                output1d.kstarDen->data[ibinKstar] += 1.0;
            }
        }
    }
    cout << endl;
}

/**
 * Simulate the correlation function in 1d by sampling a flat relative momentum distribution
 * and a Gaussian relative position distribution (3d).
 * 
 * @param output1d OutputHistograms containing the output histograms.
 * @param nMcSamples int number of Monte Carlo samples to take.
 * @param KSTAR_BINS int number of bins for the correlation function.   
 * @param KSTAR_MAX double maximum value of the correlation function.
 * @param DELKSTAR double bin width of the correlation function.
 * @param SIGNED_KSTAR bool whether to use signed k* values.
 * @param RX double source radius in x direction.
 * @param RY double source radius in y direction.
 * @param RZ double source radius in z direction.
 * 
*/
void generateEvents1dParallel( OutputHistograms& output1d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf ) {
    /*
    // Get the source parameters
    double RX = parameter::getD(sourceMap, "RX", 10.0);
    double RY = parameter::getD(sourceMap, "RY", 10.0);
    double RZ = parameter::getD(sourceMap, "RZ", 10.0);
    const double root2_RX = sqrt(2.0) * RX;
    const double root2_RY = sqrt(2.0) * RY;
    const double root2_RZ = sqrt(2.0) * RZ;

    CRandom randy(-12345);  // single instance for seeding

    // Correlation function parameters
    bool SIGNED_KSTAR = parameter::getB(kstarMap, "SIGNED_KSTAR", false);
    double DELKSTAR = parameter::getD(kstarMap, "KSTAR_DELTA", 1.0);
    int KSTAR_BINS = parameter::getI(kstarMap, "KSTAR_BINS", 200);
    double KSTAR_MAX = parameter::getD(kstarMap, "KSTAR_MAX", 200.0);

    int nMcSamples = parameter::getI(sourceMap, "nMcSamples", 1000000);
    
    // Determine the number of threads to use
    const unsigned int MAX_THREADS = 6;
    unsigned int nThreads = min(std::thread::hardware_concurrency(), MAX_THREADS);
    cout << "Using " << nThreads << " threads" << endl;
    int samplesPerThread = nMcSamples / nThreads;

    // Mutex for safely updating the main kstarNum and kstarDen
    std::mutex mtx;

    // Function to handle Monte Carlo sampling in each thread
    auto monteCarloSampling = [&](int threadID) {
        auto localKstarNum = output1d.kstarNum;
        auto localKstarDen = output1d.kstarDen;
        CRandom randy(-12345 - threadID);  // Different seed per thread

        for (int ik = 0; ik < KSTAR_BINS; ++ik) {
            double k_base = DELKSTAR * ik + DELKSTAR / 2.;

            for (int ipart = 0; ipart < samplesPerThread; ++ipart) {
                //double k_sign = SIGNED_KSTAR ? 1.0 : pow(-1, static_cast<int>(randy.ran()));
                double k_sign = 1.;
                double k = k_base;

                if (k > KSTAR_MAX) continue;

                double x = root2_RX * randy.gauss();
                double y = root2_RY * randy.gauss();
                double z = root2_RZ * randy.gauss();
                double r_sq = x * x + y * y + z * z;

                double weight = 1.0;  // Modify as needed
                int ibinKstar = localKstarNum->whatBin(k);
                if (ibinKstar >= 0 && ibinKstar < localKstarNum->ndata) {
                    localKstarNum->data[ibinKstar] += weight;    // WIP: weight is given by the wave function
                    localKstarDen->data[ibinKstar] += 1.0;
                }
            }
        }

        // Lock and update main histograms safely
        std::lock_guard<std::mutex> lock(mtx);
        for (int i = 0; i < output1d.kstarNum->ndata; ++i) {
            output1d.kstarNum->data[i] += localKstarNum.data[i];
            output1d.kstarDen->data[i] += localKstarDen.data[i];
        }
    };

    // Launch threads
    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < nThreads; ++i) {
        threads.emplace_back(monteCarloSampling, i);
    }

    // Wait for threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    */
}

/**
 * Simulate the correlation function in 3d by sampling a flat relative momentum distribution
 * and a Gaussian relative position distribution (3d).
 * 
 * @param output3d OutputHistograms containing the output histograms.
 * @param nMcSamples int number of Monte Carlo samples to take.
 * @param KSTAR_BINS int number of bins for the correlation function.
 * @param KSTAR_MAX double maximum value of the correlation function.
 * @param DELKSTAR double bin width of the correlation function.
 * @param SIGNED_KSTAR bool whether to use signed k* values.
 * @param RX double source radius in x direction.
 * @param RY double source radius in y direction.
 * @param RZ double source radius in z direction.
 * 
*/
void generateEvents3d( OutputHistograms& output3d, parameterMap sourceMap, parameterMap kstarMap, CWaveFunction& wf ) {  
    /*
    // Get the source parameters
    double RX = parameter::getD(sourceMap, "RX", 10.0);
    double RY = parameter::getD(sourceMap, "RY", 10.0);
    double RZ = parameter::getD(sourceMap, "RZ", 10.0);
    const double root2_RX = sqrt(2.0) * RX;
    const double root2_RY = sqrt(2.0) * RY;
    const double root2_RZ = sqrt(2.0) * RZ;

    //CRandom randy(-12345);  // single instance for seeding

    // Correlation function parameters
    bool SIGNED_KSTAR = parameter::getB(kstarMap, "SIGNED_KSTAR", false);
    double DELKSTAR = parameter::getD(kstarMap, "KSTAR_DELTA", 1.0);
    int KSTAR_BINS = parameter::getI(kstarMap, "KSTAR_BINS", 200);
    double KSTAR_MAX = parameter::getD(kstarMap, "KSTAR_MAX", 200.0);

    int nMcSamples = parameter::getI(sourceMap, "nMcSamples", 1000000);
    
    // Determine the number of threads to use
    unsigned int nThreads = std::thread::hardware_concurrency();
    int samplesPerThread = nMcSamples / nThreads;

    // Mutex for safely updating the main kstarNum and kstarDen
    std::mutex mtx;

    // Function to handle Monte Carlo sampling in each thread
    auto monteCarloSampling = [&](int threadID) {
        auto localKstarNum = output3d.kstarNum;
        auto localKstarDen = output3d.kstarDen;
        auto localKstarX = output3d.kstarX;
        auto localKstarY = output3d.kstarY;
        auto localKstarZ = output3d.kstarZ;
        CRandom randy(-12345 - threadID);  // Different seed per thread

        for (int ikx = 0; ikx < KSTAR_BINS; ++ikx) {
            double kx_base = DELKSTAR * ikx;

            for (int iky = 0; iky < KSTAR_BINS; ++iky) {
                double ky_base = DELKSTAR * iky;

                for (int ikz = 0; ikz < KSTAR_BINS; ++ikz) {
                    double kz_base = DELKSTAR * ikz;

                    for (int ipart = 0; ipart < samplesPerThread; ++ipart) {
                        double kx_offset = SIGNED_KSTAR ? 1.0 : pow(-1, static_cast<int>(randy.ran()));
                        double ky_offset = SIGNED_KSTAR ? 1.0 : pow(-1, static_cast<int>(randy.ran()));
                        double kz_offset = SIGNED_KSTAR ? 1.0 : pow(-1, static_cast<int>(randy.ran()));

                        double kx = kx_offset * (kx_base + DELKSTAR * randy.ran());
                        double ky = ky_offset * (ky_base + DELKSTAR * randy.ran());
                        double kz = kz_offset * (kz_base + DELKSTAR * randy.ran());
                        double k_sq = kx * kx + ky * ky + kz * kz;

                        if (k_sq > (KSTAR_MAX * KSTAR_MAX)) continue;

                        double x = root2_RX * randy.gauss();
                        double y = root2_RY * randy.gauss();
                        double z = root2_RZ * randy.gauss();
                        double r_sq = x * x + y * y + z * z;

                        //double weight = 1.0;  // Modify as needed
                        double weight = wf.CalcPsiSquared(sqrt(k_sq), sqrt(r_sq), 0.0);

                        int ibinKstar = localKstarNum.whatBin(sqrt(k_sq));
                        int ibinKstarX = localKstarX.whatBin(kx);
                        int ibinKstarY = localKstarY.whatBin(ky);
                        int ibinKstarZ = localKstarZ.whatBin(kz);
                        if (ibinKstar >= 0 && ibinKstar < localKstarNum.ndata) {
                            localKstarNum.data[ibinKstar] += weight;    // WIP: weight is given by the wave function
                            localKstarDen.data[ibinKstar] += 1.0;
                            localKstarX.data[ibinKstarX] += weight;
                            localKstarY.data[ibinKstarY] += weight;
                            localKstarZ.data[ibinKstarZ] += weight;
                        }
                    }
                }
            }
        }

        // Lock and update main histograms safely
        std::lock_guard<std::mutex> lock(mtx);
        for (int i = 0; i < output3d.kstarNum.ndata; ++i) {
            output3d.kstarNum.data[i] += localKstarNum.data[i];
            output3d.kstarDen.data[i] += localKstarDen.data[i];
            output3d.kstarX.data[i] += localKstarX.data[i];
            output3d.kstarY.data[i] += localKstarY.data[i];
            output3d.kstarZ.data[i] += localKstarZ.data[i];
        }
    };

    // Launch threads
    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < nThreads; ++i) {
        threads.emplace_back(monteCarloSampling, i);
    }

    // Wait for threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    */
}
