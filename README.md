# Halo_Sparsity

This code computes the average halo sparsity at a given redshift for a given cosmological model by solving Eq. (8) in Balmes et al. (2014) assuming a Sheth-Tormen parametrization of the halo mass function at M500c and M1000c masses which has been calibrated against the RayGal simulation halo catalogs. The input linear matter power spectrum for a given set of cosmological parameters is given by the Eisenstein & Hu (1998) formulae of the linear transfer function, while the linearly extrapolated spherical collapse threshold is given by the formula from Kitayama & Suto (1996). 

# Installation

Modify the Makefile to specify your Fortran compiler, after compiling the code run the executable, sparse.
