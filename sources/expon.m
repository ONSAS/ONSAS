function R = expon(t);
  I = eye(3,3);  al = norm(t) ;
  if al == 0
    R = I ;
  else
    Rsk = skew(t);
    R = I+sin(al)/al*Rsk+2*(sin(al/2)/al)^2*Rsk^2;
  end
end
