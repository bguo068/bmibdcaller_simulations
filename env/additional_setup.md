## Setup `isorelate`

1. Activate the Conda environment: 
	`conda activate simulation`
2. Prep Conda environment for compilation: 
	`mamba install r-devtools`
3. install `isorelate` in R (select none for updates). Note a commit hash is
   used to control the version of `isorelate`. Within an R session of
   `simulation` environment, run:
	```
	devtools::install_github(
		"bahlolab/isoRelate@109ee470f1d1e9e4cd131045c23eea7047cff649"
	)
	```

## setup `hmmIBD`
1. Activate the simulation environment and change directory: 
	```
	conda activate simulation
	cd $CONDA_PREFIX
	```
2. Compile and install:
	```
	git clone https://github.com/glipsnort/hmmIBD.git
	cd hmmIBD/
	git checkout v2.0.4
	sed -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.667e-7/' -i hmmIBD.c
	$CXX -o hmmIBD -O3 -lm -Wall hmmIBD.c -fpermissive
	mv hmmIBD $CONDA_PREFIX/bin
	cd ..
	rm -rf hmmIBD
	```

## setup TPBWT

1. Activate the simulation environment and change directory: 
	```
	conda activate simulation
	cd $CONDA_PREFIX
	```
2. The `prep Conda environment for compilation` step is covered in Conda environment recipe
	`mamba install cython pandas numpy setuptools_scm`
3. git clone and install `phasedibd` into the Conda environment

	```
	git clone https://github.com/23andMe/phasedibd.git
	cd phasedibd
	python setup.py build_ext --inplace
	python setup.py install
	python tests/unit_tests.py
	cd ..
	rm -rf phasedibd
	```
