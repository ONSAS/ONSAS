
# Introduction

![tests](https://github.com/ONSAS/ONSAS.m/workflows/tests/badge.svg)
[![License](https://img.shields.io/badge/License-GPLv3-green.svg)](https://github.com/ONSAS/ONSAS.m/blob/master/COPYING.txt)
[![Join the chat at https://gitter.im/onsas_/community](https://badges.gitter.im/onsas_/community.svg)](https://gitter.im/onsas_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## What is ONSAS.m?

ONSAS.m is an Open Nonlinear Structural Analysis Solver for GNU-Octave/Matlab. It consists in a set of implementations of numerical methods for static/dynamic and linear/non-linear analysis of structures. The first version was developed for educational purposes and published in a Structural Analysis [handbook](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf). The current version allows to perform a diverse set of simulations and it's been used in specific research applications.

## Publications using ONSAS

Journal articles using ONSAS:

 * 2022, M. Vanzulli, J. M. Pérez Zerpa, A consistent co-rotational formulation for aerodynamic nonlinear analysis of flexible frame structures, arXiv article under-review (https://arxiv.org/abs/2204.10545)
 * 2022, M. Forets, D. Freire, J. M. Pérez Zerpa, *Combining set propagation with finite element methods for time integration in transient solid mechanics problems*, Computers & Structures [url](https://www.sciencedirect.com/science/article/abs/pii/S0045794921002212?dgcid=coauthor)

Theses using ONSAS:

 * 2021, M. Vanzulli, *Implementación de una formulación corrotacional en dinámica no lineal y aplicación al modelado de líneas de transmisión eléctrica* [url](https://www.colibri.udelar.edu.uy/jspui/handle/20.500.12008/28388)


## Some example applications

### A deployable ring

```@raw html
<img src="https://github.com/ONSAS/ONSAS.m/blob/master/docs/src/assets/deployableRing.gif?raw=true">
```

[ring](https://github.com/ONSAS/ONSAS.m/blob/master/docs/src/assets/deployableRing.gif?raw=true)


### A simple wind turbine model

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/wind.gif?raw=true" alt="wind turbine animation">
```
[wind turbine animation](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/wind.gif?raw=true)

### A tower model

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true" alt="tower animation">
```

[tower](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true)

### A uniaxial extension test

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/uniaxial.gif?raw=true" alt="uniaxial animation">
```
[uniaxial animation](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/uniaxial.gif?raw=true)

### A beam/truss pendulum

```@raw html
<img src="https://raw.githubusercontent.com/mvanzulli/Tex_CorrotationalDynamicTL_TesisMV/main/Presentacion/Videos/3.gif" alt="beam truss pendulum">
```

### A chain model

```@raw html
<img src="https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true" alt="chain animation">
```
[chain](https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true)


## Contributors and License

### License

The code is distributed under a [GNU-GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.html).

### Authors

The authorship of each version is (or tends to be) based on the criteria defined by [the JOSS journal](https://joss.readthedocs.io/en/latest/submitting.html#authorship). The co-authors have collaborated in tasks such as: design, development or extensive documentation contributions.

* [**Jorge M. Pérez Zerpa**](https://scholar.google.com.uy/citations?user=Qb476KIAAAAJ&hl=en) (**1**), leaded and managed the design and development of the code, developed the assembly functions, nonlinear truss element formulation, nonlinear static analysis function, designed and co-authored Newmark's method function, input and output functions, leaded the generation of the documentation.

* [**Mauricio Vanzulli**](https://github.com/mvanzulli) (**2**) co-developed the Newmark's method functions and scripts, developed input files for the dynamic analysis examples. Developed the nonlinear dynamic co-rotational frame element function, its validation and integration with the VIV function. 

* [**Alexandre Villié**](https://www.linkedin.com/in/alexandre-villi%C3%A9-343870187/) (**3**) developed the current Vortex-Induced-Vibrations Wake Oscillator model. Also contributed in the validation of this function in the integration with the co-rotational frame element.

* [**Joaquín Viera**](https://exportcvuy.anii.org.uy/cv/?b6b1cd2fe90a9c29279eedb0d3cc4c4d) (**1**), leaded the development of the Linear Analysis module and input files, collaborated in the design and development of the input reading and output generation modules, leaded the development of a GUI.

* [**J. Bruno Bazzano**](https://uy.linkedin.com/in/juan-bruno-bazzano-garc%C3%ADa-a045bb56) (**1**), contributed to the design/development of the buckling analysis modules, co-designed the initial version of the code, developed and implemented validation examples, validated the HHT implementation.

* [**Marcelo Forets**](https://scholar.google.fr/citations?user=XSJzDEsAAAAJ&hl=en) (**4**) made relevant contributions to application of software technology tools, in particular, the documentation generation workflow. Contributed also to the implementation of the Neo-Hookean solid model.

* [**Jean-Marc Battini**](https://scholar.google.com/citations?user=7dzVcKoAAAAJ&hl=en) (**5**), contributed functions associated with the computation of static internal forces of the nonlinear frame element.

**Affiliations**:

 1. Instituto de Estructuras y Transporte, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay
 1. Instituto de Ingeniería Mecánica y Producción Industrial, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
 1. Department of Mechanical Engineering, Polytechnique Montréal, Montréal, Canada.
 1. Centro Universitario Regional Este, Universidad de la República, Maldonado, Uruguay
 1. Department of Civil and Architectural Engineering, KTH Royal Institute of Technology, Stockholm, Sweden

### Contributions and Acknowledgments

[Santiago Correa](https://github.com/santiago-correa-89) made contributions to examples, currently available at the documentation site.

The functions in `linearStiffMatPlate3D.m` and `assemblyUniform.m` use part of the [fem_plate_example.m](https://www.fing.edu.uy/~jorgepz/files/fem_plate_example.m) code developed by Jorge Pérez Zerpa and [**Pablo Castrillo**](https://www.fing.edu.uy/~pabloc/). Professor [**Sebastian Toro**](https://scholar.google.com/citations?user=7Z3ruPAAAAAJ&hl=es), provided functions for reading dxf files, which were part of ONSAS until version 0.2.6.

Prof. Pérez Zerpa would like to thank: Prof. [Frédérick Gosselin](https://fgosselin.meca.polymtl.ca/?lang=en) for his support during the initial contributions of Alexandre Villié, Prof. [Eduardo de Souza Neto](https://scholar.google.com/citations?user=Yrk2yIMAAAAJ&hl=en) for his comments on the arc-length norm computation, and [**Pablo Blanco**](https://scholar.google.com/citations?user=X0382ScAAAAJ&hl=es),
 [**Gonzalo Ares**](https://scholar.google.com/citations?user=lCeQOH0AAAAJ&hl=en) and [**Gonzalo Maso Talou**](https://unidirectory.auckland.ac.nz/profile/g-masotalou) for so many discussions during early stages of the design of the code.

The development of ONSAS has been partially supported by funds provided by the following agencies/projects:
 - Comisión de Investigación Científica (CSIC) (project: *Definición de estrategias para la aplicación de métodos de identificación de material al diagnóstico no invasivo de Cáncer de mama*, manager, Prof. Pérez Zerpa),
 - Comisión Sectorial de Enseñanza (project: *Rediseño de prácticas de enseñanza y evaluación en Resistencia de Materiales*, manager, Prof. Pérez Zerpa),
 - Agencia Nacional de Investigación e Innovación (project VIOLETA, code `FSE_1_2016_1_131837`, manager, Prof. Gabriel Usera).


## Contact

You can send an e-mail to _jorgepz[AT]fing.edu.uy_ or join the chat in the gitter [chat room](https://gitter.im/onsas_/community).
