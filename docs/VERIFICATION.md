# Verification roadmap and numerical status

This document distinguishes **tests already represented by scripts/figures** from **checks that have been demonstrated to pass under a reproducible run**. It was drafted from the repository source and artifacts; simulations have **not** been independently rerun for this documentation update.

## Verification matrix

| Feature | Available implementation/evidence | Next quantitative check | Current status |
|---|---|---|---|
| Local bulk evolution | `Task0/task_0.py`, `Task0/Testfvsphi.py` | Compare numerical trajectory with spatially uniform reference ODE solution | Development evidence present; quantitative reproduction pending |
| Diffuse interface | `Task1/task_1.py`, `Task1/Outputs/energy_track.csv`, `thickness_track.csv` | Establish energy behavior, interface width and mesh/time-step convergence | Data and figures present; new automated pass/fail test pending |
| Two-field thermal coupling | `Task_2/task_2.py`, `Task_2/Test/test2.py`, `Task_2/out_task2/diagnostics.csv` | Correct discretization discrepancy and rerun no-flux integral-balance checks | **Needs source correction and reproducible verification** |
| Boundary conditions | `Task_2/Test/test_bc.py` | Check discrete boundary flux or integrated balance to defined tolerance | Diagnostic script exists; validate implementation and tolerance |
| Four-fold anisotropy | `Final.py`, `Task_3/Final.py`, morphology figures | Compare expected symmetry/directions, isotropic limit and mesh sensitivity | Demonstration present; quantitative benchmark pending |
| Dendrite tip kinetics | `Task_3/Test/tipvelocity.py`, `Task_3/outputs/*/tip_trace.csv` | Inspect sampling/thresholding, define tip-speed convergence window, compare to cited benchmark on matching model | Traces/plots present; benchmark agreement not established |

## Important source-level consistency checks

### 1. Backward-Euler latent-heat term

The PDE written in `Final.py` is

```math
\partial_t u = D\nabla^2 u + K\partial_t\phi.
```

For test function `w_u` and an implicit Euler step, a standard no-flux weak residual *after multiplication by* `Δt` is

```math
R_u = \int_\Omega (u^{n+1}-u^n)w_u\,dx
      - K\int_\Omega (\phi^{n+1}-\phi^n)w_u\,dx
      + \Delta t D\int_\Omega \nabla u^{n+1}\cdot\nabla w_u\,dx.
```

The current `Task_2/task_2.py` code contains `- dt*0.5*(phi - phi_0)*w_u*dx`; the root `Final.py` likewise contains `- dt*k*(phi - phi_0)*w_u*dx`. This introduces an **additional `dt`** relative to the documented PDE. Because `K=0` in the root benchmark, that term is inactive there. It *is active* in Task 2 and needs to be corrected and rechecked before making a conservation claim.

**Suggested acceptance test:** With homogeneous no-flux boundary conditions, recompute the quantity

```math
E_H(t) = \left[\int_\Omega u(t)\,dx - \int_\Omega u(0)\,dx\right]
         - K\left[\int_\Omega \phi(t)\,dx - \int_\Omega \phi(0)\,dx\right],
```

and report absolute and normalized drift as a function of time, mesh size and `Δt`. Select and document tolerances *before* calling the test a pass. Do not assume the existing `diagnostics.csv` was produced by the current checked-in source without a reproducible rerun.

### 2. Output directory, imports and solver outcomes

- `Final.py` opens XDMF files in `out_task3_6/` but does not create that directory. Make it before running, or add a rank-aware `mkdir` and barrier in the solver.
- `Task_2/task_2.py` imports `test2` and `test_bc` while those scripts reside in `Task_2/Test/`. Set the module path or restructure imports before a clean run.
- Record the actual nonlinear solve convergence flag/iteration counts instead of relying on visually plausible field evolution or a generic success message.
- Check behavior for floating-point end-time stepping and MPI-dependent initialization of interface noise.

### 3. Benchmark scope

The root `Final.py` has `zet=0` and `K=0`; therefore its saved anisotropic morphology should not be described as validation of *fully coupled thermal dendritic growth*. Re-enable physically intended coefficients only after revisiting the discretization and documenting a stable benchmark setup.

## Recommended quantitative evidence to add

1. A short, automated, small-mesh smoke test that detects crashes and non-converged nonlinear steps.
2. An enthalpy-balance table with at least two time steps and two meshes, when coupling is active.
3. An isotropic limit check (`eps_an=0`) and qualitative/quantitative four-fold symmetry analysis (`eps_an>0`).
4. A reproducible tip-position definition (e.g. `φ=0` crossing), sampling intervals and a convergence study.
5. Direct comparisons with a precisely matching literature benchmark (nondimensionalization, parameters, geometry, initial conditions and measurement window all documented).

**Status convention:** "Evidence present" means data/scripts exist; "verified" should be used only after the stated acceptance criterion has been evaluated on a documented run.
