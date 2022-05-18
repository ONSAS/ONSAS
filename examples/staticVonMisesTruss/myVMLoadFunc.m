% von mises user load example

function f = myVMLoadFunc(t);

f = zeros(6*3,1);
f(11) = -1.5e8*t ;
