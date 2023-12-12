function [ aoa, cl, cd, cm, fstat, clinv, clfullysep ] = readFoilData(filename)
    
strucdata=importdata(filename);

aoa(:,1)=strucdata(:,1);
cl(:,1)=strucdata(:,2);
cd(:,1)=strucdata(:,3);
cm(:,1)=strucdata(:,4);
fstat(:,1)=strucdata(:,5);
clinv(:,1)=strucdata(:,6);
clfullysep(:,1)=strucdata(:,7);

end