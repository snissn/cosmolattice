# CosmoLattice

## *A modern code for lattice simulations of scalar and gauge field dynamics in an expanding universe*
### Authors: Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg

---

### 🚀 **Custom Model: φ⁴ + Einstein–Cartan Torsion Extension**

This fork adds a new model, `lphi4Torsion`, which implements a φ⁴ scalar field with:

* An additional dynamical torsion field,
* Tunable mass (`mTorsion`) and quartic self-interaction (`lambdaT`) for torsion,
* Local Einstein–Cartan–inspired coupling between the torsion field and the main scalar field momentum, parameterized by `alpha_torsion`.

**Implementation notes:**

* The new model is located in `src/models/lphi4Torsion.h` and takes three scalar fields: `[phi0, phi1, Torsion]`.
* The local torsion–momentum coupling term is implemented as a patch to `src/include/CosmoInterface/evolvers/velocityverlet.h`
  (in `kickScalar()`, for field 0 only).
* See `src/models/parameter-files/lphi4Torsion.in` for an example input file with all torsion-related parameters.

#### **Building with MPI and HDF5 (recommended for large runs & output):**

```bash
cd build
cmake -DHDF5=ON -DMODEL=lphi4Torsion -DMPI=ON ../
bash fetchall.sh MyLibs    # Installs h5f and other dependencies if needed
make cosmolattice
```

#### **Running in Parallel**

You can use MPI to run on many CPU cores:

```bash
mpirun --oversubscribe -np 16 ./lphi4Torsion input=lphi4Torsion.in save_dir=./
```

> **Note:** Adjust `-np 16` to match the number of CPU cores you want to use.
> **The number of MPI processes (`-np`) should evenly divide your lattice size** (e.g., for `N=64`, use 1, 2, 4, 8, 16, 32, or 64).
> This ensures proper domain decomposition and optimal parallel performance.

**Key new input parameters:**

* `alpha_torsion`: Coupling strength (set to 0 to recover standard φ⁴ theory)
* `mTorsion`: Torsion mass parameter
* `lambdaT`: Torsion quartic self-coupling
* `initial_amplitudes`, `initial_momenta`: Now take 3 values: `[phi0, phi1, torsion]`

See the comments in `lphi4Torsion.h` for detailed documentation.

---

### Documentation

- Please visit the official webpage for CosmoLattice at [cosmolattice.net](https://cosmolattice.net).
- To learn how to install and execute the code as well as how it works :  <a href=https://arxiv.org/pdf/2102.01031.pdf target="_blank" rel="noopener noreferrer" > arXiv:2102.01031</a> .
- To learn about the underlying theoretical framework: <a href=https://arxiv.org/pdf/2006.15122.pdf target="_blank" rel="noopener noreferrer" > arXiv:2006.15122</a> .

### Basic installation

*Minimal requirements:* `CMake` version 3 or above, `g++` version 5 or above, `fftw3`.

```
git clone https://github.com/cosmolattice/cosmolattice.git
cd cosmolattice   
mkdir build                     
cd build                        
cmake -DMODEL=lphi4 ../
make cosmolattice
```

This will compile the ``lphi4`` model. To run it with the default input file, you can do

``
./lphi4 input=../src/models/parameter-files/lphi4.in
``

The above commands just represent a very brief guide for the installation and execution of CosmoLattice. 
For further information, see  Appendix A of the <a href=https://arxiv.org/pdf/2102.01031.pdf target="_blank" rel="noopener noreferrer" >user-manual</a>.
All options of CosmoLattice, as well as how to activate them and how to install the optional external 
libraries are explained at length there.

### Credits

CosmoLattice is freely available to anyone who wants to use or modify it. However, whenever 
using CosmoLattice in your research, no matter how much (or little) you modify the code, 
<b>please cite both <a href=https://arxiv.org/pdf/2006.15122.pdf target="_blank" rel="noopener noreferrer" > arXiv:2006.15122</a> 
and <a href=https://arxiv.org/pdf/2102.01031.pdf target="_blank" rel="noopener noreferrer" > arXiv:2102.01031</a> in your papers</b>. 
