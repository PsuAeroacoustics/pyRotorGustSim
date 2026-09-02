# Unsteady BEMT Rotor-Gust Interaction

This code simulates the aerodynamic and acoustic response of a hovering rotor as it interacts with a prescribed gust. It was initially developed for assessing the effectiveness of a novel passive noise reduction strategy which incorporates narrow open-closed resonator cavities inside the rotor blade. Therefore, it also includes the ability to model and optimize the geometry of the resonators for maximum noise reduction. The repository combines:

- Blade-element/momentum-theory (BEMT) rotor aerodynamics.
- A two-term indicial-airfoil model for unsteady lift response.
- Optional resonator impedance models for filtering blade loads.
- Automatic generation of PSU-WOPWOP geometry, loading, observer, and namelist files.
- Automatic execution of PSU-WOPWOP (executable not include).
- Parameter sweeps, resonator optimization, validation scripts, and plotting utilities.

The main driver is [`src/rotor_gust_interaction.py`](src/rotor_gust_interaction.py). This is a research workflow rather than a packaged Python application: input files are JSON, generated results are written into case directories, and several scripts assume they are launched from a case directory.

## Repository Layout


Important source modules:

| Module | Purpose |
| --- | --- |
| [`src/rotor_gust_interaction.py`](src/rotor_gust_interaction.py) | Command-line simulation driver and workflow orchestration |
| [`src/geometry.py`](src/geometry.py) | Handles rotor blade geometry construction. Employs classes for each component with relevant information stored as attributes.  |
| [`src/compute_aero.py`](src/compute_aero.py) | Aerodynamic analysis (Steady trim, gust inflow, indicial response, and blade loads) |
| [`src/res_funcs.py`](src/res_funcs.py) | Includes all functions related to the resonators. |
| [`src/filter_loads.py`](src/filter_loads.py) | Filters the blade loads in accordance with the average or smeared acoustic response of the resonators. The resonator optimization is also handled here. |
| [`src/wopwop_input_configure.py`](src/wopwop_input_configure.py) | PSU-WOPWOP namelist and input-file generation |
| [`src/help_funcs.py`](src/help_funcs.py) | Case-file loading, HDF5 file I/O, WOPWOP execution |
| [`src/plot.py`](src/plot.py) | Load, filter-response, gust, and acoustic plotting utility|

## Prerequisites

## Requirements

This code has a number of dependencies that need to be setup and installed first. Some of these are provided as git submodules and are largely transparent to the user. If you are cloning a fresh repository run:
```
	git clone --recurse-submodules https://github.com/DanWeitsman/rotor_gust_interaction
```

If the repository is already cloned, but submodules were not cloned, run the following command to acquire the submodules:
```
	git submodule update --init --recursive
```


- Aerosandbox (https://github.com/peterdsharpe/AeroSandbox)
- numpy
- matplotlib
- scipy
- f90nml
- h5py
- A separately installed PSU-WOPWOP executable for running acoustic cases

All the dependencies can be installed as follows
```bash
pip install -r requirements.txt
```

### Python module path

The scripts in `src/` use direct imports such as `from geometry import *`, so make `src/` importable when running from another directory:

```bash
export PYTHONPATH="$PWD/src:$PWD/dependencies/pyWopwop:$PWD/dependencies/resonator:${PYTHONPATH}"
```

Run this from the repository root. Alternatively, use the relative `PYTHONPATH` in the example below.

### PSU-WOPWOP

Acoustic calculations require a separately installed PSU-WOPWOP executable named `wopwop3` on `PATH`. The repository does not include that executable. The bundled [`dependencies/pyWopwop`](dependencies/pyWopwop) code only writes input files and converts output; it does not replace the solver.

For observer grids large enough to trigger the driver’s parallel path, `mpirun` must also be installed and able to launch `wopwop3`.

Check the external dependency before an acoustic run:

```bash
command -v wopwop3
command -v mpirun
```

## Quick Start

The `cases/enlarged_blade` directory contains a complete example configuration. Start in that directory so generated case output stays beside its input files:

```bash
cd cases/enlarged_blade
source ../../.venv/bin/activate
export PYTHONPATH="../../src:../../dependencies/pyWopwop:../../dependencies/resonator:${PYTHONPATH}"
python ../../src/rotor_gust_interaction.py \
  --aero --acs \
  -input_geom geom.json \
  -input_param param.json \
  -observer_param observer_lgrid.json \
  -acs_param acs_param.json \
  -res_param sdof_dist_param_oblique.json
```

The example’s `run.sh` uses the same inputs and adds `--filt`:

```bash
cd cases/enlarged_blade
export PYTHONPATH="../../src:../../dependencies/pyWopwop:../../dependencies/resonator:${PYTHONPATH}"
bash run.sh
```

`--aero` creates or replaces the directory named by `case_name`. Because it removes that directory before recomputing, do not use it in a directory containing results you need to preserve.

## Driver Options

The driver accepts these switches:

| Option | Effect |
| --- | --- |
| `--aero` | Compute trimmed rotor aerodynamics and gust-induced unsteady loads. Required when no `saved_params.h5` exists. |
| `--acs` | Run PSU-WOPWOP and process acoustic output into HDF5. Requires `wopwop3`. |
| `--filt` or `-f` | Apply the resonator treatment to the blade loads. Normally used with `--aero`. |
| `--opt` or `-o` | Optimize resonator geometry/distribution while evaluating acoustic levels. |
| `--plot` or `-p` | Write load and, when applicable, filter-response plots into the case directory. |
| `--noncompact` or `-nc` | Apply the filter to an estimated chordwise pressure distribution instead of compact loads. |
| `-input_geom PATH` | Rotor/blade geometry JSON. Required by the current case loader. |
| `-input_param PATH` | Simulation and gust parameter JSON. |
| `-observer_param PATH` | Observer location or observer-grid JSON. |
| `-acs_param PATH` | PSU-WOPWOP environment/output flags JSON. |
| `-res_param PATH` | Resonator parameter JSON. |

The normal pipeline is:

```text
--aero  ->  trim and simulate loads  ->  generate WOPWOP inputs
--filt  ->  compute resonator response -> replace/filter loads
--acs   ->  run wopwop3              ->  convert pressure results to HDF5
--plot  ->  write diagnostic figures
```

If `--aero` is omitted, the driver loads `case_name/saved_params.h5` and can continue with acoustics or plotting. The case directory and all paths are resolved relative to the current working directory.

## Input Files

Each run is assembled from up to five JSON files. The example files in `cases/enlarged_blade` are the best starting point for new cases.

### Geometry: `geom.json`

```json
{
  "radius": 0.381,
  "origin": [0, 0, 0],
  "number_of_blades": 1,
  "AR": 6.2119,
  "theta_tw": 0,
  "theta_initial": 2,
  "r_c": 0.268,
  "airfoil": "naca0015",
  "airfoil_points": 200
}
```

The values are SI-based unless noted otherwise. Angles such as `theta_tw` and `theta_initial` are specified in degrees. `r_c` is the nondimensional root cutout, `AR` controls chord sizing, and `airfoil` is passed to AeroSandbox.

### Simulation and gust parameters: `param.json`

```json
{
  "case_name": "example_case",
  "computational_params": {
    "d_psi": 1,
    "spanwise_elements": 48,
    "airfoil_elements": 100,
    "number_of_revs": 2,
    "unsteady_loading": true
  },
  "flight_params": {
    "density": 1.2055,
    "kinematic_viscosity": 1.488e-5,
    "omega": 356.955,
    "sos": 341.7,
    "C_T_target": 0.0015
  },
  "gust_params": {
    "strength": 0.1170926,
    "peak_location": 0.25,
    "azimuthal_location": 90
  }
}
```

`unsteady_loading` currently requires one gust-location description:

- `azimuthal_location`: a gust azimuth in degrees.
- `gust_end_pnts`: endpoints describing the gust path in the rotor plane.
- `r_trace`: a trace parameter used to construct a curved gust path.

The gust model uses `strength`, `peak_location`, and the prescribed core-size/path parameters. `omega` is angular speed in radians per second; `d_psi` is the azimuth increment in degrees.

### Observer parameters

An observer file can describe a regular spherical grid:

```json
{
  "highPassFrequency": 1,
  "lowPassFrequency": 6250,
  "nt": 1440,
  "radius": 1.54305,
  "nbTheta": 41,
  "nbPsi": 31,
  "thetaMin": 90,
  "thetaMax": 270,
  "psiMin": -60,
  "psiMax": 60
}
```

It can also describe explicit observer coordinates using `xLoc`, `yLoc`, and `zLoc`, or a list of radii/angles. See `observer_lgrid.json`, `observer_sgrid.json`, and `observer_sgrid_2.json` for variants.

### Acoustic parameters: `acs_param.json`

This file is passed to the WOPWOP namelist generator. Common flags include:

- `loadingNoiseFlag`: loading-noise calculation.
- `thicknessNoiseFlag`: blade-surface thickness noise and blade geometry generation.
- `totalNoiseFlag`: total noise and blade geometry generation.
- `acousticPressureFlag`: acoustic pressure output.
- `ASCIIOutputFlag`, `OASPLdBFlag`, `spectrumFlag`, and `SPLdBFlag`: output products.

The complete set of supported fields is defined by `EnvironmentIn` in [`src/wopwop_input_generator.py`](src/wopwop_input_generator.py). Unknown JSON fields are tolerated by that class, but should not be relied upon without checking the generated namelist.

### Resonator parameters

Resonator files specify the treated chordwise and radial extents, resonator geometry, number of patches/elements, and optional staggered distributions. For example, `sdof_dist_param_oblique.json` contains:

```json
{
  "c_extents": [0.1, 0.3],
  "r_extents": [0.6, 1],
  "r_min": 0.0003,
  "r_max": 0.005,
  "L_min": 0.01,
  "L_max": 0.278892,
  "N_patches": 1,
  "N_res": 1,
  "OAR": 0.27,
  "staggered": true,
  "x": [0.000400737025, 0.144439679, 0.165988495,
        0.0736621405, 0.172696404, 0.189297069, 0.441024621]
}
```

The meaning and expected length of `x` depend on `N_res`, `N_patches`, and `staggered`. Use the existing files in `cases/enlarged_blade` as templates and inspect [`src/res_funcs.py`](src/res_funcs.py) before defining a new parameterization.

## Generated Output

For `case_name: "example_case"`, a run creates `example_case/` in the current working directory. Typical contents include:

```text
example_case/
├── saved_params.h5              Aerodynamic/load state and run metadata
├── cases.nam                    WOPWOP case-list namelist
├── acoustics/
│   ├── example_case.nam         Generated WOPWOP case namelist
│   ├── lifting_line_geometry.dat
│   ├── loading_blade_0.dat      One loading file per blade
│   ├── observer.ascii            Generated for explicit observers/grids
│   ├── pressure.h5              Processed WOPWOP output
│   └── pressure/, spl/          Solver output directories, depending on flags
├── Fz_tseries.png               Written by plotting workflows
├── dFz_tseries.png
└── ...
```

`saved_params.h5` is a recursive HDF5 representation of the simulation dictionary and is used to resume acoustics or plotting without repeating the aerodynamic calculation. WOPWOP result data are imported from `acoustics/pressure.h5` by `import_results_from_wopwop`.

## Parameter Sweeps

The top-level generators create per-case JSON files and a `run.sh` script from a base `param.json` and `geom.json` in the current directory:

```bash
cd cases/enlarged_blade/sweeps
python ../../../Mg_Mt_param_sweep_input_file_generator.py
bash run.sh
```

Available generators vary combinations of gust strength/core size and tip Mach number:

- `Mg_Mt_param_sweep_input_file_generator.py`: gust strength versus tip Mach number.
- `rg_Mt_param_sweep_input_file_generator.py`: gust core size versus tip Mach number while preserving the strength/core-size ratio.
- `rg_Mg_param_sweep_input_file_generator.py`: gust core size versus gust strength.

The generated script assumes `rotor_gust_interaction.py` is executable or on `PATH`. To use the source tree directly, replace that command with `python ../../../../src/rotor_gust_interaction.py` as appropriate for the sweep directory, or create a small wrapper that sets `PYTHONPATH`.

## Operating-Condition Optimization

[`src/opt_op_conditions.py`](src/opt_op_conditions.py) wraps repeated rotor-gust simulations in a bounded local optimization. It varies nondimensional gust strength, gust core size, and tip Mach number to maximize the measured acoustic difference between untreated and filtered runs.

Run it from a case directory with the same five input-file options used by the main driver:

```bash
cd cases/enlarged_blade
export PYTHONPATH="../../src:../../dependencies/pyWopwop:../../dependencies/resonator:${PYTHONPATH}"
python ../../src/opt_op_conditions.py \
  -input_geom geom.json \
  -input_param param.json \
  -observer_param observer_lgrid.json \
  -acs_param acs_param.json \
  -res_param sdof_dist_param_oblique.json
```

The optimizer rewrites the input parameter file as it evaluates candidates. Work from a copy when preserving the original case definition matters.

## Indicial-Airfoil Validation

The [`af_gust_interaction`](af_gust_interaction) directory contains smaller scripts for inspecting the isolated gust/airfoil response and comparing model variants against `cfd_data.csv`:

```bash
python af_gust_interaction/indicial_af.py
python af_gust_interaction/indicial_af_val.py
```

These scripts generate figures in the current working directory and are useful for validating the indicial-response portion independently of the rotor workflow.

## Post-Processing

Scripts in [`post`](post) import WOPWOP HDF5 data and produce figures such as OASPL carpets and parameter-sweep overlays. Many are study-specific and contain an explicit `cases_directory` and case name near the top of the file. Update those paths before running them:

```bash
python post/oaspl_carpet.py
python post/Mg_Mt_sweep_carpet.py
```

The plotting code enables LaTeX text rendering and may require a local LaTeX installation and the `Times New Roman` font. If a plotting script fails before data processing, disable or adjust those Matplotlib settings for the local environment.

## Troubleshooting

### `ModuleNotFoundError` for `geometry`, `wopwop`, or `resonator`

Run from a case directory with `src`, `dependencies/pyWopwop`, and `dependencies/resonator` on `PYTHONPATH`, as shown above.

### `FileNotFoundError` for `saved_params.h5`

Run the same case once with `--aero`. A run without `--aero` assumes that `case_name/saved_params.h5` already exists.

### `wopwop3: command not found`

Install PSU-WOPWOP separately and place its executable on `PATH`. The Python helper package in this repository cannot run acoustic propagation by itself.

### Existing results disappear after a rerun

This is expected when `--aero` is used: the driver removes the complete `case_name` directory before creating new output. Rename or copy the old case first.

### Parallel acoustic execution fails

The driver selects `mpirun wopwop3` for sufficiently large observer grids. Test serial execution first, then verify that your MPI installation can launch the same `wopwop3` binary.

### A gust run fails during setup

With `unsteady_loading: true`, ensure that `gust_params` contains exactly one supported path description: `azimuthal_location`, `gust_end_pnts`, or `r_trace`, along with the required gust strength and peak/core parameters.

## Development Notes

- Inputs and generated solver files use SI units unless a field or plotting label says otherwise.
- Acoustic OASPL calculations use a reference pressure of $20\,\mu\mathrm{Pa}$ in the analysis code.
- The repository currently has no top-level automated test suite or packaging metadata.
- Generated case directories and `dependencies/` are ignored by the repository’s `.gitignore`; keep large numerical outputs outside version control unless they are intentionally being archived.

## License and Citation

No top-level license or citation metadata is currently included. Check with the project authors before redistributing the code or using it in published work. PSU-WOPWOP and the bundled helper projects retain their own provenance and usage requirements.
