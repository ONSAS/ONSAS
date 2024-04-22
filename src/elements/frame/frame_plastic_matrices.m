function [ Kfd, Kfalfa, Khd, Khalfa, Fint] = frame_plastic_matrices(E, Ks, A, l, uvector, xpi, wpi, Mnp1, kp_np1, Cep_np1, Ghats)
    
Kfd    = zeros(6, 6) ;
Kfalfa = zeros(6, 1) ;
Khd    = zeros(1, 6) ;
Khalfa = 0 ;

Bu = [-1/l 1/l] ;

npi = length(kp_np1) ;

for ip = 1:npi

  N = bendingInterFuns (xpi(ip), l, 2) ;

  Bv = [N(1) N(3)] ;  Btheta = [N(2) N(4)] ;
  
  Bd = [ Bu   0 0 0 0    ; ...
         0  0 Bv  Btheta ] ;

  Kfd     = Bd'*[E*A 0; 0 Cep_np1(ip)]*Bd ;
  
  Kfalfa  = Bd'*[E*A 0; 0 Cep_np1(ip)]*[0 Ghats(ip)]' ;
  
  Khd     = [0 Ghats(ip)]*[E*A 0; 0 Cep_np1(ip)]*Bd ;
  
  Khalfa  = Ghats(ip)*Cep_np1*Ghats(ip) ;
  
  epsilon = Bu*uvector ;
  
  Fi      = Bd' * [E*A*epsilon; Mnp1(ip)] ;

  % stiffness matrices / integration (Gauss-Lobatto)
  Kfd    = Kfd    + Kfd    * wpi(ip) ;
  Kfalfa = Kfalfa + Kfalfa * wpi(ip) ;
  Khd    = Khd    + Khd    * wpi(ip) ;
  Khalfa = Khalfa + Khalfa * wpi(ip) ;
    
  % internal forces / integration (Gauss-Lobatto)
  Fint = Fi + Fi*wpi(ip) ;

end

Khalfa = Khalfa + Ks ;

end