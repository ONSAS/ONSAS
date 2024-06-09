
# Introduction

![tests](https://github.com/ONSAS/ONSAS/workflows/octave_tests/badge.svg)
[![License](https://img.shields.io/badge/License-GPLv3-green.svg)](https://github.com/ONSAS/ONSAS/blob/master/COPYING.txt)


## What is ONSAS?

[ONSAS](https://github.com/ONSAS/ONSAS) is an Open Nonlinear Structural Analysis Solver for GNU-Octave/Matlab. It consists in a set of implementations of numerical methods for static/dynamic and linear/non-linear analysis of structures. The first version was developed for educational purposes and published in a Structural Analysis [handbook](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf) for a graduate course taught at the [School of Engineering of Universidad de la República](https://www.fing.edu.uy/). The current version allows to perform a diverse set of simulations and it has been used in specific research applications.

## Publications using ONSAS

Journal articles using ONSAS:

 * 2024, A. Villié, M. Vanzulli, J.M. Pérez Zerpa, J. Vétel, S. Etienne, F.P. Gosselin, *Modeling vortex-induced vibrations of branched structures by coupling a 3D-corotational frame finite element formulation with wake-oscillators*, **Journal of Fluids and Structures** [url](https://doi.org/10.1016/j.jfluidstructs.2024.104074)
 * 2023, M. Vanzulli, J. M. Pérez Zerpa, *A co-rotational formulation for quasi-steady aerodynamic nonlinear analysis of frame structures*, **Heliyon** [url](https://doi.org/10.1016/j.heliyon.2023.e19990)
 * 2022, M. Forets, D. Freire, J. M. Pérez Zerpa, *Combining set propagation with finite element methods for time integration in transient solid mechanics problems*, **Computers & Structures** [url](https://www.sciencedirect.com/science/article/abs/pii/S0045794921002212?dgcid=coauthor)

Theses using ONSAS:

 * 2021, M. Vanzulli, *Implementación de una formulación corrotacional en dinámica no lineal y aplicación al modelado de líneas de transmisión eléctrica* [url](https://www.colibri.udelar.edu.uy/jspui/handle/20.500.12008/28388)
 * 2021, A. Teliz, *Optimización de torres de alta tensión y su análisis frente a vientos de alta intensidad* [url](https://hdl.handle.net/20.500.12008/35985)


## Some academic examples

### A deployable ring

This nonlinear static problem is introduced in [(Goto et. al, 1992)](https://doi.org/10.1016/0020-7683(92)90024-N) and also considered in [(Battini and Pacoste, 2002)](https://doi.org/10.1016/S0045-7825(01)00352-8).

```@raw html
<img src="https://github.com/ONSAS/ONSAS/blob/master/docs/src/assets/deployableRing.gif?raw=true">
```

### Right-Angle Cantilever

This nonlinear dynamic analysis problem is introduced in [(Simo and Vu-Quoc, 1988)](https://doi.org/10.1016/0045-7825(88)90073-4) and also considered as Example 1 in [(Le et. al., 2014)](https://doi.org/10.1016/j.cma.2013.11.007).

```@raw html
<img src="https://github.com/ONSAS/ONSAS/blob/master/docs/src/assets/rightAngleCantilever.gif?raw=true" alt="right-angle animation">
```

### A simple propeller model

This problem is based on one of the examples presented in [(Vanzulli and Pérez Zerpa, 2023)](https://doi.org/10.1016/j.heliyon.2023.e19990).

```@raw html
<img src="https://github.com/ONSAS/ONSAS/blob/master/docs/src/assets/propeller.gif?raw=true" alt="propeller animation">
```

### A tower model

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true" alt="tower animation">
```

### A beam/truss pendulum

```@raw html
<img src="https://raw.githubusercontent.com/mvanzulli/Tex_CorrotationalDynamicTL_TesisMV/main/Presentacion/Videos/3.gif" alt="beam truss pendulum">
```

## Contact

You can publicly post in the [discussion section](https://github.com/ONSAS/ONSAS/discussions) or contact privately sending an e-mail to _jorgepz [AT] fing.edu.uy_ .

## Contributors and License

### License

The code is distributed under a [GNU-GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.html).

### Authors

The authorship of each version is (or tends to be) based on the criteria defined by [the JOSS journal](https://joss.readthedocs.io/en/latest/submitting.html#authorship). The co-authors have collaborated in tasks such as: design, development or extensive documentation contributions.

* [**Jorge M. Pérez Zerpa**](https://scholar.google.com.uy/citations?user=Qb476KIAAAAJ&hl=en) (**1**), leaded and managed the design and development of the code, developed the assembly functions, nonlinear truss element formulation, nonlinear static analysis function, designed and co-authored Newmark's method function, input and output functions, leaded the generation of the documentation.

* [**Mauricio Vanzulli**](https://github.com/mvanzulli) (**2**) co-developed the Newmark's method functions and scripts, developed input files for the dynamic analysis examples. Developed the nonlinear dynamic co-rotational frame element function, its validation and integration with the VIV function. 

* [**Joaquín Viera**](https://exportcvuy.anii.org.uy/cv/?b6b1cd2fe90a9c29279eedb0d3cc4c4d) (**1**), leaded the development of the Linear Analysis module and input files, collaborated in the design and development of the input reading and output generation modules, leaded the development of a GUI.

* [**Alexandre Villié**](https://www.linkedin.com/in/alexandre-villi%C3%A9-343870187/) (**3**) developed the current Vortex-Induced-Vibrations Wake Oscillator model. Also contributed in the validation of this function in the integration with the co-rotational frame element.

* [**Marcelo Forets**](https://scholar.google.fr/citations?user=XSJzDEsAAAAJ&hl=en) (**4**) made relevant contributions to application of software technology tools, in particular, the documentation generation workflow. Contributed also to the implementation of the Neo-Hookean solid model.

* [**Jean-Marc Battini**](https://scholar.google.com/citations?user=7dzVcKoAAAAJ&hl=en) (**4**), contributed functions for the computation of static internal forces of the nonlinear co-rotational frame element.

* [**J. Bruno Bazzano**](https://uy.linkedin.com/in/juan-bruno-bazzano-garc%C3%ADa-a045bb56) (**1**), contributed to the design/development of the buckling analysis modules, co-designed the initial version of the code, developed and implemented validation examples, validated the HHT implementation.

**Affiliations**:

 1. Instituto de Estructuras y Transporte, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay
 1. Instituto de Ingeniería Mecánica y Producción Industrial, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
 1. Department of Mechanical Engineering, Polytechnique Montréal, Montréal, Canada.
 1. Centro Universitario Regional Este, Universidad de la República, Maldonado, Uruguay
 1. Department of Civil and Architectural Engineering, KTH Royal Institute of Technology, Stockholm, Sweden

### Contributions and Acknowledgments

[Santiago Correa](https://github.com/santiago-correa-89) made contributions to examples, mainly beamLinearVibration, and [corrected the expression of a parameter](https://github.com/ONSAS/ONSAS/pull/699) of the $\alpha$-HHT method. Professor [**Sebastian Toro**](https://scholar.google.com/citations?user=7Z3ruPAAAAAJ&hl=es), provided functions for reading dxf files, which were part of ONSAS until version 0.2.6.

Prof. Pérez Zerpa would like to thank: Prof. [Frédérick Gosselin](https://fgosselin.meca.polymtl.ca/?lang=en) for his support during the initial contributions of Alexandre Villié, Prof. [Eduardo de Souza Neto](https://scholar.google.com/citations?user=Yrk2yIMAAAAJ&hl=en) for his comments on the arc-length norm computation, and [**Pablo Blanco**](https://scholar.google.com/citations?user=X0382ScAAAAJ&hl=es), [**Gonzalo Ares**](https://scholar.google.com/citations?user=lCeQOH0AAAAJ&hl=en) and [**Gonzalo Maso Talou**](https://unidirectory.auckland.ac.nz/profile/g-masotalou) for so many discussions during early stages of the design of the code.

The development of ONSAS has been partially supported by funds provided by the following agencies/projects:
 - Comisión de Investigación Científica (CSIC) (project: *Definición de estrategias para la aplicación de métodos de identificación de material al diagnóstico no invasivo de Cáncer de mama*, manager, Prof. Pérez Zerpa),
 - Comisión Sectorial de Enseñanza (project: *Rediseño de prácticas de enseñanza y evaluación en Resistencia de Materiales*, manager, Prof. Pérez Zerpa),
 - Agencia Nacional de Investigación e Innovación (project VIOLETA, code `FSE_1_2016_1_131837`, manager, Prof. Gabriel Usera).

