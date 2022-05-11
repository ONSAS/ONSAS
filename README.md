
# ONSAS.m: an Open Nonlinear Structural Analysis Solver for GNU-Octave/Matlab

![tests](https://github.com/ONSAS/ONSAS.m/workflows/tests/badge.svg)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://onsas.github.io/ONSAS.m/dev/)
[![License](https://img.shields.io/badge/License-GPLv3-green.svg)](https://github.com/ONSAS/ONSAS.m/blob/master/COPYING)
[![Release](https://img.shields.io/github/v/release/ONSAS/ONSAS?color=yellow&include_prereleases)](https://github.com/ONSAS/ONSAS.m/releases)
[![Join the chat at https://gitter.im/onsas_/community](https://badges.gitter.im/onsas_/community.svg)](https://gitter.im/onsas_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/170120732.svg)](https://zenodo.org/badge/latestdoi/170120732)

## Table of Contents
1. [About ONSAS](#aboutonsas)
1. [For users](#howtouseonsas)
1. [Contributions](#contributions)
1. [Contact](#contact)

## About ONSAS <a name="aboutonsas"></a>

### What is ONSAS?

ONSAS is a [GNU-Octave](https://www.gnu.org/software/octave/) code for static/dynamic and linear/non-linear analysis of structures. The first version was developed for educational purposes and was published in a [handbook](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf) of the course _Análisis no lineal de Estructuras_ taught at [Facultad de Ingeniería](https://www.fing.edu.uy/), Universidad de la República.

### What can ONSAS be used for?

The current version allows to perform dynamic/static nonlinear analyses of beam/truss/solid 3D structures. A reduced list of features is listed at next:

* **Elements** 2-node truss, 2-node Bernoulli frame, 4-node tetrahedron.
* **Static analysis methods** Newton-Raphson Method and Cylindrical Arc-Length Method.
* **Dynamic analysis methods** Newmark Method and α-HHT method.
* **Loads** nodal loads, time-history user-defined loading program.

## Some examples

### A wind turbine model
![wind](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/wind.gif?raw=true)

### A truss tower model
![tower](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true)

### A chain model
![chain](https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true)

## License

The code is distributed under a GNU-GPL 3.0 license.



## How to use ONSAS <a name="howtouseonsas"></a>


The user should follow these steps to install and run onsas:

1. Download and install [GNU-Octave](https://www.gnu.org/software/octave/) and [Paraview](https://www.paraview.org/)
1. Download the ONSAS zip source files from the realeases web https://github.com/ONSAS/ONSAS/releases
1. Open GNU-Octave, move to the _examples_ directory and run one of the examples.

## Contributions <a name="contributions"></a>

### Authors of code
The following authors collaborated in various tasks including: design, development and testing of the code, with specific contributions given by the commits history: [**Jorge M. Pérez Zerpa**](https://www.fing.edu.uy/~jorgepz) <sup>1</sup>, **J. Bruno Bazzano**<sup>1,2</sup>, [**Joaquín Viera**](https://www.researchgate.net/profile/Joaquin_Viera_Sosa) <sup>1</sup>, **Mauricio Vanzulli** <sup>3</sup> and [**Marcelo Forets**](https://scholar.google.fr/citations?user=XSJzDEsAAAAJ&hl=en)<sup>4</sup>.

The following authors contributed with specific functions of files, of great relevance for the code: [**Jean-Marc Battini**](https://scholar.google.com/citations?user=7dzVcKoAAAAJ&hl=en)<sup>5</sup> (contributed functions related to the computation of static internal forces of the nonlinear frame element) and [**Sebastian Toro**](https://scholar.google.com/citations?user=7Z3ruPAAAAAJ&hl=es)<sup>6</sup> , provided the functions: f_LectDxf.m, f_ValGrCode.m and f_XData.m, used in the dxf import function.

Affiliations:

1. Instituto de Estructuras y Transporte, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay
1. Bazzano & Scelza Ingenieros, Montevideo, Uruguay
1. Instituto de Ingeniería Mecánica y Producción Industrial, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
1. Departamento de Matemática y Aplicaciones, Centro Universitario Regional del Este, Universidad de la República, Maldonado, Uruguay
1. Department of Civil and Architectural Engineering, KTH Royal Institute of Technology, Stockholm, Sweden
1. CIMEC Santa Fe, Argentina

### Contributions and Acknowledgments

The complete list of contributions and acknowledgments is available at the documentation site.

## Contact <a name="contact"></a>

If you have any question you can send an e-mail to _jorgepz [AT] fing.edu.uy_ or join the chat at   <a href="https://gitter.im/onsas_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge" target="_blank"><img alt="Gitter" src="https://img.shields.io/gitter/room/JuliaReach/Lobby?style=for-the-badge&logo=gitter&logoColor=white" /></a>
</p>
