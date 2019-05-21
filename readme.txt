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

There are three user-specified functions:
1) qoi_func.m -- quantities of interest as function of dfba states
2) tdisc_func.m -- discontinuity time when using ns-pce
3) likelihood_func.m -- likelihood function for Bayesian estimation

There are two main codes:
1) main_pce.m -- builds the pce using either global or ns-pce, see the "User input" section for parameters that can be modified by the user
2) main_smc.m -- given a surrogate model (built and saved from "main_pce"), this script performs sequential monte carlo with resampling as well as forward uncertainty propagation to get predictive distributions. 

The mat-file "nsPCE_tol=1e-3.mat" is included, which contains the ns-pce models for glucose, xylose, and biomass over 8 time points. This can be used to replicate results in the companian paper.  

To cite the ns-pce method, please cite the following paper.

@article{paulson19,
  title={Fast uncertainty quantification of dynamic flux balance analysis models using non-smooth polynomial chaos expansions},
  author={J. A. Paulson and M. Martin-Casas and A. Mesbah},
  journal={PLoS Computational Biology},
  pages={under review},
  year={2019}
}
