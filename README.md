## Overview

This repo has a few simple test problems for measuring performance. The three test problems are:

1. Transient Thermal
2. Solid Dynamics
3. Definite Maxwell 

Of course, these test problems aren't actually going through the rigamarole of solving a realistic physics problem,
they just set up the kind of calculations that would appear in these simulations (i.e. "mass" + "stiffness" + "source" term).

The two main kernels of interest are the jacobian-vector product and the element jacobian calculation.

## Building

---

<details>
  <summary>libCEED (optional) </summary>
  <br>

  libCEED sources are included in this repo, and if you want to build with libCEED enabled you have to first `cd` into that directory and run
  
  ```sh
  make OPT='-O3 -ffp-contract=fast'
  ```
  or, if you're on an arm mac machine
  ```sh
  make OPT='-O3 -mcpu=apple-m1 -fvectorize -fslp-vectorize -ffp-contract=fast'
  ```

  Once that's done, configure the main CMake project with the additional argument `-DMFEM_USE_CEED=ON`
  <br>

</details>

---

Start by configuring CMake

```bash
cmake . -Bbuild -GNinja -DCMAKE_BUILD_TYPE=Release
```

Here are some other options that may be relevant

- `-DENABLE_CUDA=ON`
- `-DCMAKE_CUDA_CUDA_ARCHITECTURES=70`

Following that, build the project

```sh
cd build
ninja
```

and run one of the tests: `poisson`, `elasticity`, or `maxwell`. 
The arguments for specifying input meshes, refinement, and other parameters are explained by the executables.

> Note: some arguments are problem specific (e.g. `-pa` is supported by `poisson`, but not `maxwell`).
