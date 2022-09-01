% This file contains the parameters fo the hydro_frame_force example
% For Example 5
accDir = pwd ;
addpath( genpath( [ accDir '/../../../ONSAS.m/src'] ) );
addpath( genpath( [ accDir '/../source'] ) );
%
% Time parameters
%
%load('ZsolV4Nelem=20_FT=100_E=5e+10_dt=0.0028_vmax=0.2_eps=0.3_A=12.mat');
%--------
% Numerical parameters
%
dt = 0.0028 ; finalTime = 100*dt; numElements = 3;
%--------
% Fluid properties
%
nuFluid = 1e-6; va = 0.035;; rhoFluid = 1000;  
%
nameFuncVel = 'windVelUniform'; 
nameLiftFunc = 'liftCoef'; % 0.3
nameDragFunc = 'dragCoef'; % 1.2
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
%--------   
%
% Material and geometric properties
% material
E = 5e10;
rho = 2*rhoFluid; 
nu = .3; 
cu = 0; % structure damping
%
l = 1 ; d = 0.001; I = pi * d^4 / 64 ;
%
St = 0.2; 
fw0 = St*va/d ; %natural shedding frequency
Tw0 = 1/fw0;
Re = d*va/nuFluid;
%UrLeclercq = (St*(l^2)*va/d)*sqrt((ms+ma)/(E*I));
%--------