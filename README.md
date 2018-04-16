# Halo_Sparsity

This code computes the average halo sparsity at a given redshift for a given cosmological model by solving Eq. (4) in Corasaniti et al. (2018) assuming a Sheth-Tormen parametrization of the halo mass function at M500c and M1000c masses which has been calibrated against the RayGal simulation halo catalogs. The input linear matter power spectrum for a given set of cosmological parameters is given by the Eisenstein & Hu (1998) formulae of the linear transfer function, while the linearly extrapolated spherical collapse threshold is given by the formula from Kitayama & Suto (1996). For more details see Corasaniti et al. (2018).

# Installation & Running

Modify the Makefile to specify your Fortran compiler. Compile with command

> make

run the code calling the executable

> ./sparsity 

Enter the values of: Omega_m, Omega_b h2, h, sigma8, n_2, w0, wa and the redshift output of the halo sparsity (values should be separated by a coma or entered one by one touching the return key after each entry) Namelist files for LCDM-WMAP7 and Planck cosmological parameters can be found in the NAMELIST folder, in such a case run as

> ./sparsity < NAMELIST/par_XXX.nml

The code can be easily modified to compute the sparsity over a discretized redshift interval z_min < z z_max, and it can be hacked with a bit of work to include other mass function parametrizations at M500c and M1000c (see corresponding functions in module MF in mf_commons.f90) and also extended to include sparsity definitions for other overdensity thresolds.
If you have any question or find bugs or having problem do not hesitate to contact Pier-Stefano Corasaniti (Pier-Stefano.Corasaniti _at_ obspm.fr)

# Acknowledgements

If you use the code please cite: 

Balmes et al., Mont. Not. Roy. Astron. Soc. 437, 2328 (2014), arXiv:1307.2922\
Corasaniti et al., arXiv:1711.00480

## Author

* **Pier-Stefano Corasaniti** (LUTH, CNRS & Observatory of Paris) 
