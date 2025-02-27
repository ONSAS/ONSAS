# Installation

## Running ONSAS

The user should follow these steps to install and run ONSAS:

1. Download and install the latest version of [GNU-Octave](https://www.gnu.org/software/octave/).
1. Download the zip file of the latest ONSAS release from [these site](https://github.com/ONSAS/ONSAS/releases/latest).
1. Open GNU-Octave and run one of the example scripts from the examples folder (or create yours!).

## Visualizing results

You can process the outputs using Octave, however, the open-source software [ParaView](https://www.paraview.org/) can be used to visualize the results produced by ONSAS.

## Generation of geometries/meshes

The user can provide the geometry of the structure using two optional formats: .msh or .dxf.  [GMSH](https://gmsh.info/) is an open-source tool that allows to generate high-quality meshes. The dxf files can be used using any CAD tool.

## Development tools

If you want to contribute you should install python -> poetry -> miss-hit, to check the styling of the code before pushing. It is **recommended but not mandatory** to install poetry in a specific directory doing this in linux's bash or [whatever corresponds to you OS](https://docs.python.org/3/library/venv.html): 

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



