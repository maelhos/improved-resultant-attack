
# Implementation of the attack described in "Improved Resultant Attack against Arithmetization-Oriented Primitives"

## Build instructions

### Cloning the repo and building NTL

Building this project requires *g++*, *GMP* and *libssl-dev* as well as a working [SageMath](https://github.com/sagemath/sage) installation for the Python scripts.

First clone the repo with the submodule :

```bash
git clone --recurse-submodules https://github.com/maelhos/polymerize.git
```

**Note that if you downloaded this repo as a zip or didn't use `--recurse-submodules` you must manually run `git clone https://github.com/libntl/ntl.git` in the main directory.**

Then, patch and build the NTL :

```bash
chmod +x patchNTL.sh
./patchNTL.sh
cd ntl/src
./configure
make clean && make && sudo make install
```

*Patching the NTL is required for all instances above a few rounds. If you get a "Polynomial too big for FFT" error, check that the following lines are correctly set in `ntl/include/FFT.h` :*

```c++
#if (36 <= NTL_FFTMaxRootBnd)
#define NTL_FFTMaxRoot (36)
```

Note that we use the NTL with static linking.
After building the NTL you should have it as a static library in `ntl/src/ntl.a`

### Building the attack binaries

- Anemoi :

```bash
make clean && make anemoi
```

- Griffin :

```bash
make clean && make griffin
```

- Rescue-Prime :

```bash
make clean && make rescue
```

## Configuration

### Attack configuration

#### Anemoi

For *Anemoi*, everything can be configured in **C++** as no round-skipping is involved :
You can change the $\alpha$ value and the seed for the generation of the constants in `anemoi/include/anemoi.h`.
The number of attacked rounds and the prime field size can be modified in `anemoi/src/main.cpp` as the **R** and **P** macros.

The Python scripts can be used to check for consistency with the C++ implementation.
Currently, a demo of an Anemoi11 CICO vector is checked in the `anemoi/python/anemoi.py` file.

#### Griffin

For *Griffin*, we provide a helper `griffin/python/griffin_helper.py`. If you want to change the parameters, you must do it at the top of this Python file, run it and paste the generated header in `griffin/include/griffin_mat.h` between the *// CONSTANTS :* and *// END CONSTANTS* tags before compiling.

#### Rescue

Same as *Griffin* but we also provide `rescue/python/rescue.py` to check the CICO afterward.
