% Numerical parameters
%
dt = 0.0028 ; finalTime = 10*dt; numElements = 3;
%--------
% Fluid properties
%
nuFluid = 1e-6; rhoFluid = 1000;  
%
nameFuncVel = 'windVelUniform'; 
nameLiftFunc = 'liftCoefVIV'; 
nameDragFunc = 'dragCoef'; 
cL0 = feval(nameLiftFunc, 0, 1); cD0 = feval(nameDragFunc, 0, 1); 
%--------   
%
% material
E = 5e10;
rho = 2*rhoFluid; 
nu = .3; 
cu = 0; 
%
l = 1 ; d = 0.001; I = pi * d^4 / 64 ;
%
St = 0.2; 
%--------
% With WOMV3
    % uzsol = 1.0e-07 *[  0    0         0         0         0         0         0            0         0         0         0;
    %                     0    0.0012    0.0058    0.0144    0.0264    0.0409    0.0569       0.0738    0.0907    0.1073    0.1232;
    %                     0    0.0011    0.0052    0.0126    0.0227    0.0347    0.0485       0.0642    0.0817    0.1010    0.1216;
    %                     0    0.0010    0.0050    0.0127    0.0229    0.0349    0.0491       0.0659    0.0852    0.1060    0.1275];
% With WOMV4
uzsol = 1.0e-07 *[ 0         0         0         0         0         0         0        0         0         0         0;
                   0   -0.0002   -0.0008   -0.0012   -0.0008    0.0004    0.0020        0.0035    0.0044    0.0042    0.0026;
                   0    0.0005    0.0023    0.0050    0.0082    0.0110    0.0127        0.0131    0.0131    0.0138    0.0153;
                   0    0.0011    0.0047    0.0111    0.0214    0.0356    0.0524        0.0719    0.0970    0.1284    0.1637];
