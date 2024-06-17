# Creating structural models

The data and properties of each structural model are defined through a set of definitions in a .m script. These properties are stored in [struct](https://octave.org/doc/v5.2.0/Structures.html#Structures) data structures. The following structs must be defined and provided as input to the ONSAS function in this order:

 1. `materials`
 1. `elements`
 1. `boundaryConds`
 1. `initialConds`
 1. `mesh`
 1. `numericalMethod`
 1. `otherParams`

Each struct has its own _fields_ with specific names, used to store each corresponding property or information. Each field is obtained or assigned using _structName.fieldName_. A description of each struct and its fields follows at next. It is highly recommended to read the current sectiong following one of the examples presented in the documentation. 

## The `materials` struct

The materials struct contains the information of the material behavior considered for each element.

### `material.modelName`

This is field contains a string of the material model used to compute the internal forces of the structure. The models implemented in ONSAS are:

 * `'elastic-linear'`: used for linear behavior with small strains and small displacements. The scalar parameters of this model are $p_1=E$ the Young modulus and $p_2=\nu$ the Poisson's ratio.

 * `'elastic-SVK'`: used for a Saint-Venant-Kirchhoff material where the parameters $p_1$ and $p_2$ are the Lamé parameters with the strain-energy density function given by the following equation (where $\textbf{E}$ is the Green-Lagrange strain tensor)
```math
\Psi( \textbf{E} ) = \frac{p_1}{2} tr(\textbf{E})^2 + p_2 tr(\textbf{E}^2)
\quad
p_1 = \frac{ E \nu }{ (1+\nu) (1-2\nu) }
\quad
p_2 = \frac{ E }{ 2 (1+\nu) }
```

 * `'elastic-NHC'`: used for a Neo-Hookean compressible material. The model implemented is given by
```math
\Psi( \textbf{C} ) = \frac{p_1}{2} ( tr(\textbf{C})-3 -2 L( \sqrt{det(\textbf{C})} ) ) + \frac{p_2}{2} \left( \sqrt{det(\textbf{C})}-1 \right)^2
 \quad
 p_1 = \frac{ E }{ 2 (1+\nu) }
 \quad
 p_2 = \frac{ E }{ 3 (1-2 \nu) }
```

 * `'elastic-rotEngStr'`: used for 1D elements (truss or frame) under large displacements.

 * `'elastic-rotLogStr'`: used for 1D elements (truss) under large displacements.

 * `'plastic'`: an ElastoPlastic material with isotropic hardening given by the von mises flow rule for the plane strain element. The parameters are introduced as: REVISAR!! $p_1=E$ , $p_2 = K$ and $p_3=\sigma_{Y,0}$.

 * `'plastic-rotEngStr'`: an ElastoPlastic material .

 * `'plastic-rotLogStr'`: an ElastoPlastic material .

### `materials.modelParams`

A cell structure with vectors with the material properties of each material used in the model. The $i$-th entry of the cell, contains a vector like this:
```math
[ p_1 \dots p_{n_P} ]
```
where $n_P$ is the number of parameters of the constitutive model and $\mathbf{p}$ is the vector of constitutive parameters.

### `material.density`

This is a cell with the scalar values of the densities of the materials used in the model.

### `material.nodalMass`

This fields sets a vector of nodal masses components $[m_x, m_y, m_z]$ that is assigned to nodes.

## The `elements` struct

The elements struct contains the information about the type of finite elements used and their corresponding parameters.

### `elements.elemType`

A cell structure with the string-names of the elements used: `node`, `truss`, `frame`, `triangle` or `tetrahedron`. Other auxiliar types such as `edge` are also available

### `elements.elemTypeParams`
A cell structure with auxiliar params information, required for some element types:

 * `triangle` vector with parameters, the first parameter is an integer indicating if plane stress (1) or plane strain (2) case is considered.

### `elements.massMatType`
The `massMatType` field sets, for frame or truss elements, whether consistent or lumped mass matrix is used for the inertial term in dynamic analyses. The `massMatType` field should be set as a string variable: `'consistent'` or `'lumped'`,  and if it is not declared then by default the `'lumped'` mass matrix is set.
 
### `elements.elemCrossSecParams`
This is a cell structure with the information of the geometry of the element.

#### 1D elements

For `truss` or `frame` elements, this cell has two entries, first a string with a name of the type of cross section, and in the second entry a vector of real parameters setting the shape of that section:
```math
\{ crossSectionTypeString, \,\, [ crossSectionParam_{1}, \,\,\dots,\,\, crossSectionParam_{n} ] \}
```
with $n$ being the number of parameters of the cross section type, and `crossSectionTypeString` the type of cross section. The possible cross section strings and their corresponding properties are:

 - `generic`  :general sections, where areas and inertias are provided as parameters according to the vector: $[A \,\, J \,\, I_{yy} \,\, I_{zz} \,\, I_{\rho}(1,1) \,\, I_{\rho}(2,2) \,\, I_{\rho}(3,3) ] $ where $A$ is the area, $I_{ii}$ is the second moment of inertia of the cross-section respect to $i$ direction, $J$ is the polar moment of inertia and $I_{\rho}$ is the inertia tensor.
 - `rectangle`: rectangular sections where thicknesses ``t_y`` and ``t_z`` are provided as the vector $[t_y, t_z]$
 - `circle` : circular sections where diameter is provided.
 - `pipe` : circular hollow section where external and internal diameters are provided as first and second entries of the vector of elementCrossSecParams.

For `edge` elements the thickness is expected (for 2D load computations).

See the `crossSectionProps.m` function for more details.



#### 2D elements

For 2D elements such as `triangle` in this field a float number representing the thickness of the element is set.   

### `elements.aeroCoefFunctions`
If a frame aerodynamic analysis is desired, the drag, lift and pitch moment functions should be defined using this field. This field should contain a cell with either the strings of the functions or the definition of anonymous functions for draf lif and pitch moment in that order. Each function must receive as first input the incidence angle and as second the Reynolds number. For some `elemCrossSecParams` like `'circle'` internal built-in functions are set as default thus there is no need to set this field.


### `elements.chordVector`
A vector with the three coordinates of the aerodynamic chord vector (the system of coordinates considered for this is the local reference system at the undeformed configuration)


## The `boundaryConds` struct

### `boundaryConds.loadsCoordSys`
cell containing the coordinates system for the loads applied in each BC, each entry should be a `'global'` string or a `'local'`, or an empty array if no load is applied in that BC setting `[]`.

### `boundaryConds.loadsTimeFact`
cell with the inline function definitions of load factors of the loads applied of an empty array.

### `boundaryConds.loadsBaseVals`
cell with the (row) vector of the components of the load case
```math
[ f_x,  \, m_x, \, f_y, \, m_y, \, f_z, \, m_z ]
```
where $f_i$ are the components of forces and $m_i$ are the moments. Both forces or moments are considered per unit of length in the case of `truss`/`frame`/`edge` elements, or per unit of area in the case of `triangle`.

### `boundaryConds.userLoadsFileName`
string with the filename of the `.m` function file provided by the user that can be used to apply forces not given by time-varying loadFactors. This function file should be placed in the example folder and it must receive two arguments:  t (the time) and UsCell (a cell with: {the current displacement, velocity and acceleration} ). The function should one forces vector with the size of all the degrees of freedom of the problem (in global coordinates).

### `boundaryConds.imposDispDofs`
cell with vectors of the local degrees of freedom imposed (integers from 1 to 6)

### `boundaryConds.imposDispVals`
cell with vectors of the values of displacements imposed.

### `boundaryConds.springDofs`
vector with the local degrees of freedom of the node with springs (integers from 1 to 6)

### `boundaryConds.springVals`
vector with the values of the springs stiffnesses.

## The `mesh` struct

The mesh struct contains the finite element mesh information.

### `mesh.nodesCoords`
matrix with the coordinates of all the nodes of the mesh. The $i$-th row contains the three coordinates of the node $i$: $[x_i , \, y_i ,\, z_i]$,

### `mesh.conecCell`
[cell array](https://octave.org/doc/v5.2.0/Cell-Arrays.html) with the elements and node-connectivity information. The $\{i,1\}$ entry contains the vector with the MEB (Material, Element, boundaryConds) indexes and the nodes of the $i$-th element. The structure of the vector at each entry of the cell is:
```math
 [ materialInd, \, elementInd, \, boundaryCondInd, \, node_1 \dots node_{n} ]
```
where the first three indexes are natural numbers and $n$ is the number of nodes required by the type of element. If no property is assigned the $0$ index can be used, for instance, nodes used to introduced loads should be defined with `materialIndex = 0`.


## The `initialConds` struct

If initial conditions are homogeneous, then an empty struct should be defined using `initialConds = struct() ;`. Otherwise the fields that can be set are:

 - `initialConds.U`: a vector of the displacements at time 0.
 - `initialConds.Udot`: a vector of the velocities  at time 0.
 - `initialConds.Udotdot`: a vector of the accelerations at time 0.


## The `analysisSettings` struct

This struct contains the parameters required to apply the numerical method for the resolution of the nonlinear equations:

 * `methodName`: string with the name of the method used: `'newtonRaphson'`,`'arcLength'`, `'newmark'`,`'alphaHHT'`.
 * `stopTolDeltau`: float with tolerance for convergence in relative norm of displacements increment
 * `stopTolForces`: float with tolerance for convergence in relative norm of residual loads
 * `stopTolIts`: integer with maximum number of iterations per time step
 * `deltaT`: time step
 * `finalTime`: final time of simulation
 * `incremArcLen`: radius of the cylinder for the arcLength method. if scalar is provided then this is fixed during all times. if a vector is provided then for each time $t_i$ the entry $i$ of the vector will be used as radius. 
 * `ALdominantDOF`: if this is not set or set as `[]` (zero) then the cylindical ArcLength method based on DeSouzaNeto's Computational Methods for plasticity is used, if a non-empty vector is provided, then the dominant dof arc length variant based on Jirásek & Bazant, Inelastic Analysis of Structures, 2002, Chapter 22, is used, the first entry is the dof and the second the scaling factor
 * `deltaNM`: delta parameter of newmark method. If this parameter is not declared then the classic Trapezoidal Newmark delta = $1/2$ is set.
 * `alphaNM`: alpha parameter of newmark method. If this parameter is not declared then the classic  Trapezoidal Newmark alpha = $1/4$ is set.
 * `alphaHHT`: alpha parameter of alpha-HHT method. If this parameter is not declared then alpha=$-0.05$ is set.
 * `posVariableLoadBC`: (parameter used by the arcLength method) this parameter is an integer with the entry of the _boundaryConds_ cell corresponding with the loads vector affected by the load factor
 * `iniDeltaLamb`: (parameter used by the arcLength method) this parameter sets the initial increment for the load factor $\lambda$.

another additional optional parameters are:

 * `booleanSelfWeight`: a boolean indicating if self weight loads are considered or not. The loads are computed using the density of the material and in the $-z$ global direction.
 * `iniMatUs`: a matrix with initial solutions for each time step.
 * `addedMassBoolean`: if this parameter is set `'true'` the fluid density is considered in the intertial forces term for frame elements.

the aerodynamic-frame element parameters set are
* `fluidProps`: is a row cell with the density $\rho_f$, viscosity $\nu_f$ and the function with the fluid velocity  

```math
\{ \rho_f; \,\, \nu_f; \,\, 'fluidVelocity'\}
```
 * `geometricNonLinearAero`: a boolean to take into account geometric nonlinearities or (reconfiguration) for each element.
 * `numGaussPointsAeroForce`:  number of Gauss integration points per element for the aerodynamic forces vector. Default is 4. 
 * `computeAeroStiffnessMatrix`: a boolean to compute the aerodynamic forces stiffness matrix using a central difference algorithm. Default is `'false'`, since can affect performance.  

## The `otherParams` struct

  * `problemName`: string with the name of the problem, to be used in outputs.
  * `plots_format`: string indicating the format of the output. Use `'vtk'` for vtk output. default: no output.
  * `plots_deltaTs_separation`: integer number __N__ such that the time between vtk plots is __N x deltaT__.
  * `controlDofs`: matrix with information of the degrees of freedom to compute and control. Each row should contain this form: `[ node localdof ]`.
  * `storeBoolean`: boolean to store the results of the current iteration such as the displacements, tangent matrices, normal forces and stresses. [default: 1]
  * `nodalDispDamping`: scalar value of a linear viscous damping factor applied for all the displacement degrees of freedom [default: 0]
