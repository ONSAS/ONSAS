function [W, Wint] = uBEMdynamicInflow(Wold, Wqs, Wqsold, Wintold, VrelNorm, r, R, a, deltaT)

% Dynamic Wake Model (Snel and Schepers, 1995). Oye Model
Woldy    = dot(Wold   , [0 , 1 , 0])     ;  Woldz    = dot(Wold   , [0 , 0 , 1])   ;
Wqsy     = dot(Wqs    , [0 , 1 , 0])     ;  Wqsz     = dot(Wqs    , [0 , 0 , 1])   ;
Wqsyold  = dot(Wqsold , [0 , 1 , 0])     ;  Wqszold  = dot(Wqsold , [0 , 0 , 1])   ;
Wintyold = dot(Wintold, [0 , 1 , 0])     ;  Wintzold = dot(Wintold, [0 , 0 , 1])   ;

k = 0.6; 

Tau1 = ( 1.1 / ( 1 - 1.3*a ) ) * R / sqrt(VrelNorm);
Tau2 = ( 0.39 - ( 0.26*(r/R)^2 ) )*Tau1;

% Calculate right hand of first order differential equation using backward
% difference
Hy = Wqsy + k*Tau1*(Wqsy - Wqsyold) / deltaT;
Hz = Wqsz + k*Tau1*(Wqsz - Wqszold) / deltaT;
              
% Solve intermedial Induced Velocity solving first differential equation
Winty = Hy + ( Wintyold - Hy ) * exp(-deltaT/Tau1);
Wintz = Hz + ( Wintzold - Hz ) * exp(-deltaT/Tau1);
              
% Solve analitically second first order differential equation        
Wy = Winty + ( Woldy - Winty ) * exp(-deltaT/Tau2);
Wz = Wintz + ( Woldz - Wintz ) * exp(-deltaT/Tau2);

W    = [0, Wy, Wz]';
Wint = [0, Winty, Wintz]';
          
end