## Setup `isorelate`

1. Activate the Conda environment: 
	`conda activate bmibdcaller_simulations`
2. install `isorelate` in R (select none for updates). Note a commit hash is
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
	conda activate bmibdcaller_simulations
	cd $CONDA_PREFIX
	```
2. Compile and install:
	```
	git clone git@github.com:bguo068/hmmibd-rs.git
	cd hmmibd-rs/
	$CXX -o hmmIBDr -O3 -lm -Wall hmmIBD.c -fpermissive
	mv hmmIBDr $CONDA_PREFIX/bin
	cd ..
	rm -rf hmmIBD
	```

## setup TPBWT

1. Activate the simulation environment and change directory: 
	```
	conda activate bmibdcaller_simulations
	cd $CONDA_PREFIX
	```
2. In case, the following packages is not included in the environment
`bmibdcaller_simulations`, run:
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

## setup tskibd
```
	conda activate bmibdcaller_simulations
	git clone git@github.com:gbinux/tskibd.git
	cd tskibd
	git submodule update --init --recursive
	meson build
	ninja -C build tskibd
	cp  build/tskibd $CONDA_PREFIX/bin/
	cd ..
	rm -rf tskibd
```