%md This functions computes manually the aerodynamic loads submitted to the hole beam
function f = myLinearAeroLA(t)
global numElements
%geometric parameers
l = 20 ; d = .5 ; 
%fluid parameters
rhoA = 1.225 ; 
% inititialize external force
f = zeros( (numElements + 1)*6, 1) ;
% read the wnd velocity
windVel = feval('windVelLA', 0, t) ;
% the angle of incidence is 
betaRel = acos(dot([0 1 0] , [0 0 1] ));
% the drag, lift and moment coefficients are:
c_d = feval('dragCoefFunctionLA'  , betaRel );
c_l = feval('liftCoefFunctionLA'  , betaRel );
c_m = feval('momentCoefFunctionLA',-betaRel );
%md Then the dynamic pressures $q_0$ defined above are expressed such that: 
q = 1/2 * rhoA * (windVel(3)^2 + windVel(2)^2) ;
%md next the loads per unit of length are  
qz = q * c_d * d ; 
qy = q * c_l * d ;
qm = q * c_m * d ; 

% add forces into f vector
for elem = 1:numElements
  dofsElem = (elem - 1)*6 + 1 : (elem - 1)*6 + 12 ;
  lelem = l / numElements ;
  % torsional laods
  Tx = qm * lelem / 2       ;
  dofsMx = dofsElem(2:6:end);
  %y Loads
  Fy = -qy * lelem / 2      ;
  dofsFY = dofsElem(3:6:end);
  Mz = qy * lelem^2 / 12    ;
  dofsMz = dofsElem(6:6:end);
  %z Loads
  Fz = qz * lelem / 2       ;
  dofsFz = dofsElem(5:6:end);
  My = qz * lelem^2 / 12    ;
  dofsMy = dofsElem(4:6:end);

  % add into f vectos
  f(dofsMx) = f(dofsMx) + Tx          ; 
  f(dofsMy) = f(dofsMy) + [-My ; My ] ; 
  f(dofsMz) = f(dofsMz) + [-Mz ; Mz ] ; 
  f(dofsFY) = f(dofsFY) + Fy          ; 
  f(dofsFz) = f(dofsFz) + Fz          ; 
end
