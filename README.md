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


## Requirements

This code utilizes two other repositories that need to be setup and installed first. These are provided as git submodules and are largely transparent to the user. If you are cloning a fresh repository run:
```
	git clone --recurse-submodules https://github.com/PsuAeroacoustics/pyRotorGustSim.git
```

If the repository is already cloned, but submodules were not cloned, run the following command to acquire the submodules:
```
	git submodule update --init --recursive
```

In addition the following dependencies are required. 


- Aerosandbox (https://github.com/peterdsharpe/AeroSandbox)
- numpy
- matplotlib
- scipy
- f90nml
- h5py
- A separately installed PSU-WOPWOP executable for running acoustic cases

These can be installed as follows
```bash
pip install -r requirements.txt
```

### Python module path

Before running the script, the repository containing `src/rotor_gust_interaction.py` must be added to your system PATH environment variable so that the command can be executed from any directory. This can be accomplished by adding the following lines to your `.bashrc` or `.zshrc` file:

```bash
export PATH="/path/to/src/rotor_gust_interaction.py:$PATH"
```

After updating the file, reload your shell configuration:

```bash
source ~/.bashrc
```
or 
```bash
source ~/.zshrc
```

This only needs to be done once. Thereafter, pyPostAcs.py can be executed from any directory containing the desired acoustic datasets.

### PSU-WOPWOP

Acoustic calculations require a separately installed PSU-WOPWOP executable named `wopwop3`. The path to the exicutable must also be included in you `PATH` variable. The repository does not include that executable. The bundled [`dependencies/pyWopwop`](dependencies/pyWopwop) code only writes input files and converts output; it does not replace the solver.

For observer grids large enough to trigger the driver’s parallel path, `mpirun` must also be installed and able to launch `wopwop3`.

Check the external dependency before an acoustic run:

```bash
command -v wopwop3
command -v mpirun
```

## Quick Start

The `cases/example_case` directory contains a complete example configuration. Start in that directory so generated case output stays beside its input files:

```bash
cd cases/enlarged_blade
python ../../src/rotor_gust_interaction.py \
  --aero --acs \
  -input_geom geom.json \
  -input_param param.json \
  -observer_param observer_lgrid.json \
  -acs_param acs_param.json \
  -res_param res_params.json
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
| `-input_geom PATH` | Rotor/blade geometry JSON. Required by the current case loader. |
| `-input_param PATH` | Simulation and gust parameter JSON. |
| `-observer_param PATH` | Observer location or observer-grid JSON. |
| `-acs_param PATH` | PSU-WOPWOP environment/output flags JSON. |
| `-res_param PATH` | Resonator parameter JSON. |


If `--aero` is omitted, the driver loads `case_name/saved_params.h5` and can continue with acoustics or plotting. The case directory and all paths are resolved relative to the current working directory.

## Input Files

Each run is assembled from up to five JSON files. The example files in `cases/enlarged_blade` are the best starting point for new cases.

### Geometry: `geom.json`

```json
{
  "radius": 0.381, // rotor radius [m]
  "origin": [0, 0, 0], // x,y,z coordinates of origin 
  "number_of_blades": 1, // blade count
  "AR": 6.2119, // aspect ratio 
  "theta_tw": 0, // linear twist angle (deg)
  "theta_initial": 2, // initial collective pitch angle (adjusted by trim algorithm)
  "r_c": 0.268, // nondimensional root cutout
  "airfoil": "naca0015",  //airfoil shape
  "airfoil_points": 200  // number of points defining the sectional profile of the airfoil 
}
```

### Simulation and gust parameters: `param.json`

```json
{
  "case_name": "example_case", // case name
  "computational_params": {
    "d_psi": 1, // azimuthal resolution (deg)
    "spanwise_elements": 48, // number of spanwise elements
    "airfoil_elements": 100, // number of points defining the sectional profile of the airfoil
    "number_of_revs": 2, // number of revs to simulate
    "unsteady_loading": true // set to true to include unsteady effects from gust
  },
  "flight_params": {
    "density": 1.2055, // air density (kg/m^3)
    "kinematic_viscosity": 1.488e-5, // kinematic viscosity (m^2/s)
    "omega": 356.955, // rotational rate (rad/s)
    "sos": 341.7, // sos (m/s)
    "C_T_target": 0.0015 // target thrust coefficient
  },
  "gust_params": {
    "strength": 0.1170926,  // nondimensionalized vortex circulation strength
    "peak_location": 0.25, // nondimensionalized core size (normalized by blade chord length)
    "azimuthal_location": 90 // azimuthal position of a straight gust (also supports curved gusts)
  }
}
```

`unsteady_loading` currently requires one gust-location description:

- `azimuthal_location`: a gust azimuth in degrees.
- `gust_end_pnts`: endpoints describing the gust path in the rotor plane.
- `r_trace`: a trace parameter used to construct a curved gust path.

The gust model uses `strength`, `peak_location`, and the prescribed core-size/path parameters.

### Observer parameters

The example observer file describes a spherical grid but others can be used. The variables are defined in accordance with PSU-WOPWOP as specified in the manual:

```json
{
	"highPassFrequency":1,
	"nt":1440,
	"radius": 1.54305,
	"nbTheta": 11,
	"nbPsi": 1,
	"thetaMin": 120,
	"thetaMax":240,
	"psiMin": -30,
	"psiMax": -30
}
```

### Acoustic parameters: `acs_param.json`

This file is passed to the WOPWOP namelist generator. Common flags include:

- `loadingNoiseFlag`: loading-noise calculation.
- `thicknessNoiseFlag`: blade-surface thickness noise and blade geometry generation.
- `totalNoiseFlag`: total noise and blade geometry generation.
- `acousticPressureFlag`: acoustic pressure output.
- `ASCIIOutputFlag`, `OASPLdBFlag`, `spectrumFlag`, and `SPLdBFlag`: output products.

The complete set of supported fields is defined by `EnvironmentIn` in [`src/wopwop_input_generator.py`](src/wopwop_input_generator.py). Unknown JSON fields are tolerated by that class, but should not be relied upon without checking the generated namelist.

### Resonator parameters

Resonator files specify the treated chordwise and radial extents, resonator geometry, number of patches/elements, and optional staggered distributions. For example, `res_params.json` contains:

```json
{
  "c_extents": [0.1, 0.3], // chordwise extent of treatment coverage
  "r_extents": [0.6, 1], // spanwise extent of treatment coverage
  "r_min": 0.0003, // minimum possible resonator radius (m)
  "r_max": 0.005, // maximum possible resonator radius (m)
  "L_min": 0.01, // minimum possible resonator length (m)
  "L_max": 0.278892, // maximum possible resonator length (m)
  "N_patches": 1, // number of unique impedance patches, each consisting of different resonators
  "N_res": 1, // number of unique resonators per impedance patch
  "OAR": 0.27, // open-area ratio
  "staggered": true, // set to true to stagger impedance patches along the blade span, thus leaving sections of the blade untreated
  "x": [0.000400737025, 0.144439679] // input vector to the optimizer [resonator radius, resonator radius]. The length and entries of this array differ depending on the treatment configuraiton. Inspect [`src/res_funcs.py`] for more info. 
}
```

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

## Post-Processing

Scripts in [`post`](post) import WOPWOP HDF5 data and produce figures such as OASPL carpets and parameter-sweep overlays. Many are study-specific and contain an explicit `cases_directory` and case name near the top of the file. Update those paths before running them:

```bash
python post/oaspl_carpet.py
python post/Mg_Mt_sweep_carpet.py
```
