function [W, Wint] = BEMdynamicInflow(Wold, Wqs, Wqsold, Wintold, V0, r, R, a, deltaT)

% Dynamic Wake Model (Snel and Schepers, 1995). Oye Model
k = 0.6; 

Tau1 = ( 1.1 / ( 1 - 1.3*a ) ) * R / norm(V0);
Tau2 = ( 0.39 - ( 0.26*(r/R)^2 ) )*Tau1;

% Calculate right hand of first order differential equation using backward
% difference
H = Wqs + k*Tau1*(Wqs - Wqsold) / deltaT;
              
% Solve intermedial Induced Velocity solving first differential equation
Wint = H + ( Wintold - H ) * exp(-deltaT/Tau1);
              
% Solve analitically second first order differential equation        
W = Wint + ( Wold - Wint ) * exp(-deltaT/Tau2);
          
end