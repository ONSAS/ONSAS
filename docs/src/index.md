
# Introduction

```@raw html
<a href="https://gitter.im/onsas_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge"><img src="https://badges.gitter.im/onsas_/community.svg" alt="Join the chat at https://gitter.im/onsas_/community">
</a>
```

## What is ONSAS?

ONSAS is an Open Nonlinear Structural Analysis Solver. It consists in a set of implementations of different numerical methods for static/dynamic and linear/non-linear analysis of structures. The first version was developed for educational purposes and published in a Structural Analysis [handbook](https://www.colibri.udelar.edu.uy/jspui/bitstream/20.500.12008/22106/1/Bazzano_P%c3%a9rezZerpa_Introducci%c3%b3n_al_An%c3%a1lisis_No_Lineal_de_Estructuras_2017.pdf).

Currently different implementations and tools aimed for diverse applications are under development. The most mature is [ONSAS.m](https://github.com/ONSAS/ONSAS.m), a GNU-Octave implementation of the solver, whose user guide is described in this documentation.

### What can ONSAS be used for?

The current version allows to perform dynamic/static nonlinear analyses of beam/truss/solid 3D structures. A reduced list of features is listed at next:

* **Elements** 2-node truss, 2-node Bernoulli frame, 4-node tetrahedron.
* **Static analysis methods** Newton-Raphson Method and Cylindrical Arc-Length Method (**to be fixed!**).
* **Dynamic analysis methods** Newmark Method and $\alpha$-HHT.
* **Loads** nodal loads, time-history user-defined loading program.

### Some examples

#### A wind turbine model

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/wind.gif?raw=true" alt="wind turbine animation">
```
[wind turbine animation](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/wind.gif?raw=true)

#### A tower model

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true" alt="tower animation">
```

[tower](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/tower.gif?raw=true)

#### A uniaxial extension test

```@raw html
<img src="https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/uniaxial.gif?raw=true" alt="uniaxial animation">
```
[uniaxial animation](https://github.com/ONSAS/ONSAS_docs/blob/master/gifs/uniaxial.gif?raw=true)

#### A beam/truss pendulum

```@raw html
<img src="https://raw.githubusercontent.com/mvanzulli/Tex_CorrotationalDynamicTL_TesisMV/main/Presentacion/Videos/3.gif" alt="beam truss pendulum">
```

#### A chain model

```@raw html
<img src="https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true" alt="chain animation">
```
[chain](https://user-images.githubusercontent.com/42485529/90902313-a6bf8d80-e3a2-11ea-8369-a9be639552f9.gif?raw=true)


#### A transmission line tower model

```@raw html
<img src="https://raw.githubusercontent.com/mvanzulli/Tex_CorrotationalDynamicTL_TesisMV/main/Presentacion/Videos/4.gif" alt="transmission line">
```

## Contributors and License

The code is distributed under a [GNU-GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.html).



### Authors

The following authors collaborated in various tasks including: design, development and testing of the code.

* [**Jorge M. Pérez Zerpa**](https://scholar.google.com.uy/citations?user=Qb476KIAAAAJ&hl=en) (**1**), leaded and managed the design and development of the code, developed the assembly functions, nonlinear truss element formulation, nonlinear static analysis function, designed and co-authored Newmark's method function, input and output functions, leaded the generation of the documentation.

* **J. Bruno Bazzano** (**1,2**), leaded the development of the nonlinear/linear buckling analysis modules, co-designed the code, developed and implemented validation examples, validated the HHT implementation.

* [**Joaquín Viera**](https://exportcvuy.anii.org.uy/cv/?b6b1cd2fe90a9c29279eedb0d3cc4c4d) (**1**), leaded the development of the Linear Analysis module and input files, collaborated in the design and development of the input reading and output generation modules, leaded the development of GUI.

* **Mauricio Vanzulli** (**3**) co-developed the Newmark's method functions and scripts, developed input files for the dynamic analysis examples.

* [**Marcelo Forets**](https://scholar.google.fr/citations?user=XSJzDEsAAAAJ&hl=en) (**4**) developed the Neo-Hookean solid model.

The following authors contributed by :

* [**Jean-Marc Battini**](https://scholar.google.com/citations?user=7dzVcKoAAAAJ&hl=en) (**5**), contributed functions associated with the computation of static internal forces of the nonlinear frame element.

* [**Sebastian Toro**](https://scholar.google.com/citations?user=7Z3ruPAAAAAJ&hl=es) (**6**), provided the functions: f_LectDxf.m, f_ValGrCode.m and f_XData.m, used in the dxf import function.

### Affiliations

 1. Instituto de Estructuras y Transporte, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay
 1. Bazzano & Scelza Ingenieros, Montevideo, Uruguay
 1. Instituto de Ingeniería Mecánica y Producción Industrial, Facultad de Ingeniería, Universidad de la República, Montevideo, Uruguay.
 1. Centro Universitario Regional Este, Universidad de la República, Maldonado, Uruguay
 1. Department of Civil and Architectural Engineering, KTH Royal Institute of Technology, Stockholm, Sweden
 1. CIMEC Santa Fe, Argentina

### Contributions and Acknowledgments

The functions in `linearStiffMatPlate3D.m` and `assemblyUniform.m` use part of the
[fem_plate_example.m](https://www.fing.edu.uy/~jorgepz/files/fem_plate_example.m) code
developed by Jorge Pérez Zerpa and [**Pablo Castrillo**](https://www.fing.edu.uy/~pabloc/).

 J. M. Pérez Zerpa would like to thank: [**Pablo Blanco**](https://scholar.google.com/citations?user=X0382ScAAAAJ&hl=es)
 from the [hemolab.lncc.br](http://hemolab.lncc.br/) group at LNCC Brazil,
 [**Gonzalo Ares**](https://scholar.google.com/citations?user=lCeQOH0AAAAJ&hl=en) from Univ. Nacional de Mar del Plata, [**Gonzalo Maso Talou**](https://unidirectory.auckland.ac.nz/profile/g-masotalou) from
 the Auckland Bioengineering Institute and [**Diego Figueredo**](https://www.researchgate.net/profile/Diego_Figueredo4)
 for their numerous comments and suggestions.

 The development of ONSAS has been partially supported by funds provided by the following agencies/projects:
 Comisión de Investigación Científica (CSIC), Comisión Sectorial de Enseñanza (project: *Rediseño de prácticas de enseñanza y evaluación en Resistencia de Materiales*, manager, Prof. Pérez Zerpa), Agencia Nacional de Investigación e Innovación
 (project VIOLETA, code `FSE_1_2016_1_131837`, manager, Prof. [**Usera**](https://scholar.google.com/citations?user=9U_jEd4AAAAJ&hl=en).

### Contact

You can send an e-mail to _jorgepz[AT]fing.edu.uy_ .
