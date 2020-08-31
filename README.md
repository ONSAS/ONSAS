
# ONSAS: an Open Nonlinear Structural Analysis Solver (v. 0.1.10-RC)

[![Join the chat at https://gitter.im/onsas_/community](https://badges.gitter.im/onsas_/community.svg)](https://gitter.im/onsas_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://onsas.github.io/ONSAS_docs/dev/)


## Table of Contents
1. [About ONSAS](#aboutonsas)
1. [For users](#howtouseonsas)
1. [Contributions](#contributions)
1. [Contact](#contact)

## About ONSAS <a name="aboutonsas"></a>
------

### What is ONSAS?

ONSAS is a [GNU-Octave](https://www.gnu.org/software/octave/) code for static/dynamic and linear/non-linear analysis of structures. The first version was developed for educational purposes and was published in a [handbook](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf) of the course _Análisis no lineal de Estructuras_ taught at [Facultad de Ingeniería](https://www.fing.edu.uy/), Universidad de la República.
  
### What can ONSAS be used for?

The current version allows to perform dynamic/static nonlinear analyses of beam/truss/solid 3D structures. A reduced list of features is listed at next:

* **Elements** 2-node truss, 2-node Bernoulli frame, 4-node tetrahedron.
* **Static analysis methods** Newton-Raphson Method and Cylindrical Arc-Length Method.
* **Dynamic analysis methods** Newmark Method.
* **Loads** nodal loads, time-history user-defined loading program.

## Some examples

# A wind turbine model
![wind](https://github.com/ONSAS/ONSAS/blob/development/examples/wind.gif?raw=true)

# A chain model
![chain](https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true)

### License

The code is distributed under a GNU-GPL 3.0 license.



## How to use ONSAS <a name="howtouseonsas"></a>
------

The user should follow these steps to install and run onsas:

1. Download and install [GNU-Octave](https://www.gnu.org/software/octave/)
1. Download the ONSAS source files https://github.com/onsas/onsas/archive/v0.1.9.zip
1. Open GNU-Octave and run the _ONSAS.m_ script
1. Select one of the available input files (or create one).

We encourage the user to read the user's guide available at [https://www.fing.edu.uy/~jorgepz/onsas/mainUserGuide.html](https://www.fing.edu.uy/~jorgepz/onsas/mainUserGuide.html).


## Contributions <a name="contributions"></a>
------

### Authors
The following authors collaborated in various tasks including: design, development and testing of the code. 

* [**Jorge M. Pérez Zerpa**](https://www.fing.edu.uy/~jorgepz) <sup>1</sup>, leaded and managed the design and development of the code, developed the assembly functions, nonlinear truss element formulation, nonlinear static analysis function, designed and co-authored Newmark's method function, input and output functions, leaded the generation of the documentation.

* **J. Bruno Bazzano**<sup>1,2</sup>, leaded the development of the nonlinear/linear buckling analysis modules, co-designed the code, developed and implemented validation examples, validated the HHT implementation.

* **Joaquín Viera** <sup>1</sup>, leaded the development of the Linear Analysis module and input files, collaborated in the design and development of the input reading and output generation modules, leaded the development of GUI.

* **Mauricio Vanzulli** <sup>4</sup> co-developed the Newmark's method functions and scripts, developed input files for the dynamic analysis examples.

* [**Marcelo Forets**](https://scholar.google.fr/citations?user=XSJzDEsAAAAJ&hl=en)<sup>6</sup> developed the Neo-Hookean solid model.

The following authors contributed by :

* [**Jean-Marc Battini**](https://scholar.google.com/citations?user=7dzVcKoAAAAJ&hl=en)<sup>3</sup>, contributed functions associated with the computation of static internal forces of the nonlinear frame element.

* [**Sebastian Toro**](https://scholar.google.com/citations?user=7Z3ruPAAAAAJ&hl=es)<sup>5</sup> , provided the functions: f_LectDxf.m, f_ValGrCode.m and f_XData.m, used in the dxf import function. 

Affiliations:

1. Instituto de Estructuras y Transporte, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay
1. Bazzano & Scelza Ingenieros, Montevideo, Uruguay
1. Department of Civil and Architectural Engineering, KTH Royal Institute of Technology, Stockholm, Sweden
1. Instituto de Ingeniería Mecánica y Producción Industrial, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
1. CIMEC Santa Fe, Argentina
1. Centro Universitario Regional Este, Universidad de la República, Maldonado, Uruguay

### Contributions and Acknowledgments
The functions linearStiffMatPlate3D.m and assemblyUniform.m use part of the [fem_plate_example.m](https://gitlab.fing.edu.uy/snippets/60) code developed by Jorge Pérez Zerpa and [**Pablo Castrillo**](https://www.fing.edu.uy/~pabloc/).  J. M. Pérez Zerpa would like to thank: [**Pablo Blanco**](https://scholar.google.com/citations?user=X0382ScAAAAJ&hl=es) from the [hemolab.lncc.br](http://hemolab.lncc.br/) group at LNCC Brazil, [**Gonzalo Ares**](https://scholar.google.com/citations?user=lCeQOH0AAAAJ&hl=en) from Univ. Nacional de Mar del Plata, [**Gonzalo Maso Talou**](https://unidirectory.auckland.ac.nz/profile/g-masotalou) from the Auckland Bioengineering Institute and [**Diego Figueredo**](https://www.researchgate.net/profile/Diego_Figueredo4) for their numerous comments and suggestions. The development of this version was partially supported by funds provided by the following agencies/projects: Comisión de Investigación Científica (CSIC), Comisión Sectorial de Enseñanza ( project: _Rediseño de prácticas de enseñanza y evaluación en Resistencia de Materiales_, manager, Prof. Pérez Zerpa), Agencia Nacional de Investigación e Innovación (project VIOLETA, code FSE_1_2016_1_131837, manager, Prof. [**Usera**](https://scholar.google.com/citations?user=9U_jEd4AAAAJ&hl=en).

## Contact <a name="contact"></a>
------

You can send an e-mail to _jorgepz[AT]fing.edu.uy_ .
