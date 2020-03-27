# pfp-data-structures
Prefix-Free Parsing data structures

### Download

To clone the repository, run:

> git clone https://github.com/maxrossi91/pfp-data-structures

### Compile

We use cmake to generate the Makefile. Create a build folder in the main pfp-data-structures folder and runc cmake and than make:

> mkdir build
> cd build; cmake ..
> make

#### Debug

> mkdir build
> cd build; cmake -DCMAKE_BUILD_TYPE=Debug ..
> make

### Run

After compiling, run

>  ./test/src/pfp_ds_test ../data/yeast.fasta


### Authors ###

Theoretical results:

* Christina Boucher
* Travis Gagie
* Jan Holub
* Giovanni Manzini

Implementation:

* Ond\v{r}ej Cvacho
* Massimiliano Rossi
