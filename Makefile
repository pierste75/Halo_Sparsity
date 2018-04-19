f90comp = gfortran

sparsity: cosmo_commons.o tf_hueis.o pk_commons.o mf_commons.o sparsity_driver.o

	$(f90comp) -O3 -o sparsity cosmo_commons.o tf_hueis.o pk_commons.o mf_commons.o sparsity_driver.o

cosmo_commons.mod: cosmo_commons.o cosmo_commons.f90
	$(f90comp) -O3 -c cosmo_commons.f90
cosmo_commons.o: cosmo_commons.f90
	$(f90comp) -O3 -c cosmo_commons.f90
tf_hueis.o: tf_hueis.f
	$(f90comp) -O3 -c tf_hueis.f
pk_commons.mod: pk_commons.o pk_commons.f90
	$(f90comp) -O3 -c pk_commons.f90
pk_commons.o: pk_commons.f90
	$(f90comp) -O3 -c pk_commons.f90
mf_commons.mod: mf_commons.o mf_commons.f90
	$(f90comp) -O3 -c mf_commons.f90
mf_commons.o: mf_commons.f90
	$(f90comp) -O3 -c mf_commons.f90
sparsity_driver.o: cosmo_commons.mod pk_commons.mod mf_commons.mod sparsity_driver.f90
	$(f90comp) -O3 -c sparsity_driver.f90
clean:
	rm *.mod *.o sparsity

