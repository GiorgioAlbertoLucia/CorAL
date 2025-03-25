# PHEMTO: Proton-Helium-3 Femtoscopy  
A model for computing the correlation function of a proton–Helium-3 (p–He3) system in PbPb collisions.

### Added code
Here are listed the relevant files added/modified compared to the current version of CorAL
* include/wavefunction.h
* src/coral/Wavefunctions/wavefunction.h
* src/coral/Wavefunctions/wf/wf_pHe3_coulomb.cc
* src/phemto/phemto.cc

### Requirements  
- **ROOT** (only required for plotting)  

### Build Instructions  
To compile PHEMTO, navigate to the main **CorAL** directory and run:  
```sh
scons/scons.py
```

### Running PHEMTO  
Execute the binary from the `bin` directory:  
```sh
cd bin && ./phemto <mode> [inputFile.dat]
```
To see available options, run:  
```sh
./phemto --help
```

After running PHEMTO, visualize the results with:  
```sh
cd PHEMTO && python3 plot.py
```

Before running the simulation, ensure you update the configuration files in PHEMTO/input
