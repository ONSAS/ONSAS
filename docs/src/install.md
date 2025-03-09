# Installation

## Running ONSAS

The user should follow these steps to install and run ONSAS:

1. Download and install [GNU-Octave](https://www.gnu.org/software/octave/) or [Matlab](https://www.mathworks.com/products/matlab.html).
1. Download the zip file of the latest ONSAS release from [this site](https://github.com/ONSAS/ONSAS/releases/latest).
1. Open GNU-Octave/Matlab and run one of the example scripts from the examples folder.

## Visualizing results

You can process the outputs using Octave, however, the open-source software [ParaView](https://www.paraview.org/) can be used to visualize the results produced by ONSAS (in vtk format).

## Generation of geometries/meshes

The user can provide the geometry of the structure using GMSH's .msh format. [GMSH](https://gmsh.info/) is an open-source tool that allows to generate high-quality meshes.

## Contributing

If you want to contribute you should create a fork and create a Pull Request in github.

### Development tools

If you want to contribute any code it is recommended that you install [miss-hit](https://github.com/florianschanda/miss_hit?tab=readme-ov-file), to check the styling of the code before pushing. It is **recommended but not mandatory** to install poetry.

Doing this in linux's bash ([see other OSs](https://docs.python.org/3/library/venv.html)): 

```
cd ONSAS/utils
python -m venv .venv
source .venv/bin/activate
```

Then install poetry, install miss-hit and check

```
pip install poetry
make install
make format_check
```
