% Copyright 2025, ONSAS Authors (see documentation)
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
%
% md Function for computation of the nodal contact force and tangent stiffness matrix of a
% md frictionless node-to-rigid-plane contact element, for plane problems in the x-y plane.
% md
% md The contact constraint is regularized using the penalty method (see Wriggers,
% md Computational Contact Mechanics, 2nd ed., section 6.3). Being $\textbf{n}$ the unit
% md normal of the rigid plane, pointing towards the body, and $\textbf{x}_0$ a point of
% md that plane, the gap function of a node placed at $\textbf{x} = \textbf{X} + \textbf{u}$ is
% md```math
% md g(\textbf{u}) = \textbf{n} \cdot ( \textbf{x} - \textbf{x}_0 )
% md```
% md thus negative values of $g$ represent penetration. The penalized contact contribution to
% md the residual force vector and to the tangent stiffness matrix are given by
% md```math
% md \textbf{f} = \epsilon_N g \textbf{n} , \qquad
% md \textbf{K} = \epsilon_N \textbf{n} \otimes \textbf{n}
% md```
% md for $g \leq 0$, and zero otherwise. Since the force is returned as an internal force
% md contribution, the resulting tangent matrix is positive semi-definite.
% md
% md Inputs:
% md  - elemNodesxyzRefCoords: reference coordinates of the node of the element
% md  - elemDisps: x-y displacements of the node of the element
% md  - modelName: name of the contact model, only 'rigid-plane-penalty' is available
% md  - modelParams: [ penaltyParam  x0  y0  nx  ny ]
% md
% md Outputs:
% md  - fs: cell with the internal, viscous and inertial force vectors
% md  - ks: cell with the tangent, damping and mass matrices
% md  - gap: value of the gap function, negative if the node penetrates the rigid plane
%
function [fs, ks, gap] = contactNodeRigidPlane(elemNodesxyzRefCoords, elemDisps, modelName, modelParams)

  assert(strcmp(modelName, 'rigid-plane-penalty'), ...
         'material model not implemented for the contact-node element.');

  penaltyParam = modelParams(1);
  pointOnPlane = modelParams(2:3)';
  normalVec    = modelParams(4:5)';

  assert(penaltyParam > 0, 'the penalty parameter of the contact element must be positive.');
  assert(norm(normalVec) > 0, 'the normal vector of the rigid plane must not be null.');

  normalVec = normalVec / norm(normalVec);

  % md current position of the node in the x-y plane
  posCurr = elemNodesxyzRefCoords(1:2)' + elemDisps;

  % md gap function: negative values represent penetration of the rigid plane
  gap = normalVec' * (posCurr - pointOnPlane);

  % md the contact contribution is added only for closed gaps. The non-strict inequality
  % md keeps the stiffness of nodes exactly placed on the rigid plane, avoiding a singular
  % md system for bodies that are supported only by the contact elements.
  if gap <= 0
    Finte = penaltyParam * gap * normalVec;
    Ke    = penaltyParam * (normalVec * normalVec');
  else
    Finte = zeros(2, 1);
    Ke    = zeros(2, 2);
  end

  % md no mass nor damping is associated with the contact element
  fs = { Finte, [], zeros(2, 1) };
  ks = { Ke, zeros(2, 2), zeros(2, 2) };
