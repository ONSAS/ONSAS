### Theory and References

#### Element Formulations

ONSAS provides a diverse set of finite elements formulations for structural analysis. We strongly suggest the user to explore the references that we used to implement each one of the formulations.

For 1D structures, truss and frame elements are available, supporting both linear and large-displacement analysis via a co-rotational formulation. The nonlinear static truss co-rotational implementation is based in (Crisfield, 1997). The nonlinear static frame formulation is based on (Battini and Pacoste, 2002), while the dynamic implementation follows the consistent co-rotational approach by Le et al. (2014). For planar frame problems softening hinges are also implemented based on (Jukic, 2013).

For continuum problems, 2D (triangle) and 3D (tetrahedron) elements are implemented, allowing for Linear Elastic, Saint-Venant-Kirchhoff, and Neo-Hookean material models. For plane triangle elements, elasto-plastic analysis is available considering a formulation based on (de Souza Neto, et.al., 2008). 

#### Key References

 * (Bathe, 2014) Klaus-Jurgen Bathe.  Finite Element Procedures . 2014.
 * (Bazzano and Pérez Zerpa, 2017) J. B. Bazzano and J. Perez Zerpa.  Introducción al Análisis No Lineal de Estructuras. 2017.
 * (Battini and Pacoste, 2002) Co-rotational beam elements with warping effects in instability problems, Computer Methods in Applied Mechanics and Engineering, 191 (17-18). 2020.
 * (Crisfield, 1997) M. A. Crisfield, Non-Linear Finite Element Analysis Solids and Structure, Volume 2, Advanced Topics, , Wiley, 1997.
 * (Holzapfel, 2000) G. Holzapfel, Nonlinear Solid Mechanics, A continuum approach for Engineering, 2000, Wiley.
 * (Jukic, 2013) M. Jukić and B. Brank and A. Ibrahimbegović, Embedded discontinuity finite element formulation for failure analysis of planar reinforced concrete beams and frames, Engineering Structures, 50, 2013.
 * (Le et.al., 2014) Thanh-Nam Le and Jean-Marc Battini and Mohammed Hjiaj. A consistent 3D corotational beam element for nonlinear dynamic analysis of flexible structures, Computer Methods in Applied Mechanics and Engineering, 269, 2014.
 * (de Souza Neto, et.al., 2008) E. A. de Souza Neto and D. Perić and D. R. J. Owen, Computational Methods for Plasticity: Theory and Applications, 2008, Wiley.


