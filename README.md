# Halo_Sparsity

This code computes the average halo sparsity for massive halos (>10<sup>13</sup> M<sub>Sun</sub> h<sup>-1</sup>) at a given redshift for a given cosmological model by solving Eq. (4) in Corasaniti et al. (2018) for two different overdensity cases (<s<sub>200,500</sub>> and <s<sub>500,1000</sub>>) assuming mass function parametrizations from Despali et al. (2016) or calibrated against the RayGal simulation halo catalogs. In the case of the Despali et al. mass functions the average sparsity includes the matched halo correction calibrated on the RayGal simulations. The input linear matter power spectrum for a given set of cosmological parameters is given by the Eisenstein & Hu (1998) formulae of the linear transfer function, while the linearly extrapolated spherical collapse threshold is given by the formula from Kitayama & Suto (1996). For more details see Corasaniti et al. (2018) and Corasaniti, Sereno & Ettori (in preparation).

# Installation & Running

Modify the Makefile to specify your Fortran compiler. Compile with command

> make

run the code calling the executable

> ./sparsity 

Enter the values of: Omega<sub>m</sub>, Omega<sub>b</sub> h<sup>2</sup>, h, sigma<sub>8</sub>, n<sub>s</sub>, w<sub>0</sub>, w<sub>a</sub> and z<sub>out</sub> the redshift output of the halo sparsity (values should be separated by a coma or entered one by one touching the return key after each entry) Namelist files for LCDM-WMAP7 and Planck cosmological parameters can be found in the NAMELIST folder, in such a case run as

> ./sparsity < NAMELIST/par_XXX.nml

The code output on the screen the value of the sparsity s<sub></sub>(z<sub>out</sub>)

The code can be easily modified to compute the sparsity over a discretized redshift interval z<sub>min</sub> < z < z<sub>max</sub>, and it can be hacked with a bit of work to include other mass function parametrizations (see corresponding functions in module MF in mf_commons.f90) and also extended to include sparsity definitions for other overdensity thresolds.
If you have any question or find bugs or having problem do not hesitate to contact Pier-Stefano Corasaniti (Pier-Stefano.Corasaniti _at_ obspm.fr)

# Acknowledgements

If you use the code please cite: 

Balmes et al., Mont. Not. Roy. Astron. Soc. 437, 2328 (2014), arXiv:1307.2922\
Corasaniti et al. (2017), Astrophys. J. 862, 40 (2018), arXiv:1711.00480\
Corasaniti, Sereno, Ettori, in preparation

## Author

* **Pier-Stefano Corasaniti** (LUTH, CNRS & Observatory of Paris) 
