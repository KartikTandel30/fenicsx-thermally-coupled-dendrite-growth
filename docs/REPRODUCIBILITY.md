# Reproducibility and execution guide

The repository documents an evolving FEniCSx/DOLFINx implementation. At the time this guide was prepared, exact versions of DOLFINx, PETSc, MPI, Basix and UFL used to generate the archived output files were **not recorded in a pinned environment file**. Accordingly, the commands below are **instructions for attempting a run in an already working, compatible environment**, not proof of an independently reproduced result.

## Dependencies

### Core numerical solver

- Python and NumPy
- DOLFINx/FEniCSx, UFL and Basix
- PETSc/petsc4py
- MPI/mpi4py
- An available PETSc LU solver, such as MUMPS or SuperLU_DIST, for the current configuration

### Selected development scripts / post-processing

- Matplotlib, PyVista, PyVistaQt, OpenPyXL, meshio (depending on the script)
- ParaView for interactive viewing of XDMF output

**Compatibility notice:** The solver imports `dolfinx.nls.petsc.NewtonSolver`, so verify against the exact DOLFINx API/version available in your environment. Package instructions and APIs may change between releases. Do **not** invent or assume a package version from the mere presence of a source import.

## Windows users

The GitHub setup and code management can be performed in **Git Bash** or **PowerShell**. Running DOLFINx itself is a separate requirement; a Linux or **WSL2 Ubuntu** environment with a compatible PETSc/MPI setup is a practical path for Windows. Git Bash by itself does not install or configure DOLFINx.

## Record environment information

Run these in the Python environment in which the solver will execute, and record the actual output in a future `docs/ENVIRONMENT.md`:

```bash
python --version
python -c "import dolfinx, basix, ufl, numpy, petsc4py, mpi4py; print('dolfinx:', dolfinx.__version__); print('basix:', basix.__version__); print('ufl:', ufl.__version__); print('numpy:', numpy.__version__); print('petsc4py:', petsc4py.__version__); print('mpi4py:', mpi4py.__version__)"
mpirun --version
```

Once the original working versions are known, add a pinned and tested environment specification, preferably from the solver's native installation method.

## Run the current root solver

From the repository root, in a **working FEniCSx environment**:

```bash
mkdir -p out_task3_6
python Final.py
```

Or, where a compatible MPI installation is available:

```bash
mpirun -np 4 python Final.py
```

- `Final.py` currently uses a 250 × 250 domain and mesh divisions, `dt=0.04`, `T=200`.
- It writes `phi.xdmf` and `u.xdmf` and their associated HDF5 files to `out_task3_6/`.
- Parameters are edited in the script; there is not yet a supported configuration-file or CLI interface.
- The current parameters `zet=0`, `K=0` **disable coupling** in the showcased run.
- MPI execution and any output directory need to be accessible to the solver; do not run `mpirun` with an MPI installation mismatched to PETSc/DOLFINx.

## Development-stage entry points

In a compatible environment, with the directory shown as the current working directory:

```bash
# Task 0
cd Task0
python task_0.py

# Task 1
cd ../Task1
python task_1.py
```

The source for Task 2 currently imports test helpers from the nested `Test/` directory and expects `out_task2/` to exist. From the repository root on Linux/WSL2, a possible execution setup is:

```bash
cd Task_2
mkdir -p out_task2
PYTHONPATH=Test python task_2.py
```

**Do not interpret a completed Task 2 run as verified thermal coupling until the residual check in [`VERIFICATION.md`](VERIFICATION.md) is resolved.**

The Task 3 investigation code is in `Task_3/Final.py`; examine its parameters and output folder before running. Tip post-processing scripts are located in `Task_3/Test/` and require matching saved field files.

## What to preserve when reporting a result

- Git commit SHA of the solver and analysis scripts.
- Exact package/MPI/PETSc versions and direct/iterative solver settings.
- Domain, mesh, timestep and total simulated time.
- Parameter values (`zet`, `K`, `eps_an`, `lambda_0`, `tau_0`, `D`, `u0`, etc.).
- Initial/boundary conditions and any random seed; note the MPI distribution can affect per-rank random arrays.
- Nonlinear convergence history and exact output files used in post-processing.
- Plots and raw measurements, alongside acceptance criteria for conservation, convergence and benchmark comparisons.

The repository currently includes selected result figures/CSVs, while large XDMF/HDF5 time series are not committed.
