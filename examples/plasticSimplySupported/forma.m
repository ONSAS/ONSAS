% u = [ v1, theta2, v2, theta2]
function [xs, deformada] = forma( u , alpha, xd, L)

xs = linspace(0,L, 15)';

Ns = bendingInterFuns(xs , L, 0 ) ;

deformada = Ns*u + Mhat(xs, xd,L)*alpha; 

function  M = Mhat(x,xd,L)

M = (L-x).^2 .*  (L*x - (L+2*x) * xd) / (L^3)  + (x-xd) .*(-1+(x>xd));