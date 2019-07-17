This work is released under GPLv3 license as in GPLv3.pdf and GPLv3.tex with copyright as in copyright.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
For proper functionality, two additional features must be obtained by users:
1) the generic uncertainty quantification and sparse pce construction is based on UQlab (https://www.uqlab.com/)
2) the dfba system was integrated using DFBAlab (https://yoric.mit.edu/software/dfbalab/how-obtain-dfbalab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

These should both be downloaded and put into the matlab path before this code can be run. 
Note that any dfba integrator can be used in place of DFBAlab. All relevant files should be placed in the "dfba_model" folder.

Once a dfba integrator has been installed, the user must modify "evalForwardModel.m" to properly interact with the integrator. This function should return the states S and time T profiles for all random samples X. The states at any given time t can be determined by interpolating the values in (S,T). 

The folder structure is as follows:
1) Examples -- contains two examples of constructing surrogates and solving a parameter estimation problem
2) Functions -- contains all functions needed to execute the scripts

The included user guide provides a walkthrough of the synthetic metabolic network problem implemented in Example2. Make sure the "Functions" folder is in your matlab path. The main_map.m script requires the non-smooth optimizer solvopt (https://imsc.uni-graz.at/kuntsevich/solvopt/download.html#Matlab). 

Note that the mat-file "nsPCE_tol=1e-3.mat" is included in Example1, which contains the ns-pce models for glucose, xylose, and biomass over 8 time points. This can be used to replicate results in the companion paper.  

When using this code, please cite the following paper.

@article{paulson19,
  title={Fast uncertainty quantification of dynamic flux balance analysis models using non-smooth polynomial chaos expansions},
  author={J. A. Paulson and M. Martin-Casas and A. Mesbah},
  journal={PLoS Computational Biology},
  pages={under review},
  year={2019}
}
