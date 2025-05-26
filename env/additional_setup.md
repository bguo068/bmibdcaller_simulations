## Setup `isorelate`

1. Activate the Conda environment: 
	`conda activate bmibdcaller_simulations`
2. Prep Conda environment for compilation: 
	`mamba install r-devtools`
3. install `isorelate`  in R (select none for updates). Note a commit hash is
   used to control the version of `isorelate`. Within an R session of
   `bmibdcaller_simulations` environment, run:
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
2. Compile, patch and install `hmmIBD`:
	```
	git clone https://github.com/glipsnort/hmmIBD.git
	cmd hmmIBD
	git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
	# patch to allow specify recombination rate (-r) through command line
	patch hmmIBD.c PATH_PROJ_ROOT_DIR/env/hmmIBD.patch
	$CC -o hmmIBD -O3 -lm -Wall hmmIBD.c
	cd ..
	cp hmmIBD/hmmIBD bin/hmmIBD
	rm -rf hmmIBD
	```

## setup TPBWT

1. Activate the simulation environment and change directory: 
```bash
	conda activate bmibdcaller_simulations
	cd $CONDA_PREFIX
```
2. In case, the following packages is not included in the environment
`bmibdcaller_simulations`, run:
	`mamba install cython pandas numpy setuptools_scm`
3. git clone and install `phasedibd` into the Conda environment

```bash
	git clone https://github.com/23andMe/phasedibd.git
	cd phasedibd
	git checkout 9a7b9492b1024996cd5beadb863a7ed314d63077
	python setup.py build_ext --inplace
	python setup.py install
	python tests/unit_tests.py
	cd ..
	rm -rf phasedibd
```

## setup tskibd and tskibd-filter
```bash
	conda activate bmibdcaller_simulations
	git clone https://github.com/bguo068/tskibd.git
	cd tskibd
	git checkout v0.0.2
	git submodule update --init --recursive
	meson build
	ninja -C build tskibd
	cp  build/tskibd $CONDA_PREFIX/bin/
	# assuming cargo/rust is installed.
	cd tskibd-filter
	cargo build --release
	cp target/release/tskibd-filter $CONDA_PREFIX/bin/
	cd ..
	cd ..
	rm -rf tskibd
```

## setup ishare/ibdutils
Refer to https://github.com/bguo068/ishare?tab=readme-ov-file#requirements
```bash
	git clone https://github.com/bguo068/ishare.git
	cd ishare
	git checkout v0.1.11
	cargo build --release --bin ibdutils
	cd ..
	cp ishare/target/release/ibdutils bin/
	rm -rf ishare
```
