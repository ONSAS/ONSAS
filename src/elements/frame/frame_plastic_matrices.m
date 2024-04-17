

function [ Kfd, Kfalfa, Khd, Khalfa] = frame_plastic_matrices(E, A, )


    
Kfd    = zeros(6, 6) ;
Kfalfa = zeros(6, 1) ;
Khd    = zeros(1, 6) ;
Khalfa = 0 ;

Bu = [-1/l 1/l] ;

npi = length(kp) ;
Ms  = zeros(npi, 1) ;
tM  = 0 ;

for ip = 1:npi

  N = bendingInterFuns (xpi(ip), l, 2) ;

  Bv = [N(1) N(3)] ;  Btheta = [N(2) N(4)] ;
  
  Bd = [ Bu   0 0 0 0    ; ...
         0  0 Bv  Btheta ] ;

  Kfd     = Bd'*[E*A 0; 0 Cep_np1(ip)]*Bd ;
  
  Kfalfa  = Bd'*[E*A 0; 0 Cep_np1(ip)]*[0 Ghats(ip)]' ;
  
  Khd     = [0 Ghats(ip)]*[E*A 0; 0 Cep_np1(ip)]*Bd ;
  
  Khalfa  = Ghats(ip)*Cep*Ghats(ip) ;
  
  epsilon = Bu*uvector ;
  
  Fi      = Bd' * [E*A*epsilon; M1xpi] ;

  % stiffness matrices / integration (Gauss-Lobatto)
  Kfd    = Kfd    + Kfdj    * wpi(ii) ;
  Kfalfa = Kfalfa + Kfalfaj * wpi(ii) ;
  Khd    = Khd    + Khdj    * wpi(ii) ;
  Khalfa = Khalfa + Khalfaj * wpi(ii) ;
    
  % internal forces / integration (Gauss-Lobatto)
  Fint = Fint + Fi*wpi(ii) ;

end

Khalfa = Khalfa + Ks ;
