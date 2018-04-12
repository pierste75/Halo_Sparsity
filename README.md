# Halo_Sparsity

This code computes the average halo sparsity at a given redshift for a given cosmological model by solving Eq. (4) in Corasaniti et al. (2018) assuming a Sheth-Tormen parametrization of the halo mass function at M500c and M1000c masses which has been calibrated against the RayGal simulation halo catalogs. The input linear matter power spectrum for a given set of cosmological parameters is given by the Eisenstein & Hu (1998) formulae of the linear transfer function, while the linearly extrapolated spherical collapse threshold is given by the formula from Kitayama & Suto (1996). 

# Installation

Modify the Makefile to specify your Fortran compiler, after compiling the code run the executable, sparsity. Namelist files for LCDM-WMAP7 and Planck cosmological parameters can be found in the NAMELIST folder. The code can be easily modified to include other mass function parametrizations at M500c and M1000c, see corresponding functions in module MF in mf_commons.f90. It can also be easily extended to include sparsity definition for other overdensity thresolds.

# Acknowledgements

If you use the code please cite: 

Balmes et al., Mont. Not. Roy. Astron. Soc. 437, 2328 (2014), arXiv:1307.2922
Corasaniti et al., arXiv:1711.00480

## Author

* **Pier-Stefano Corasaniti** (LUTH, CNRS & Observatory of Paris) 
