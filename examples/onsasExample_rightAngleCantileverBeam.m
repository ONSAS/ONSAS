% ------------------------------------------------------------------------------
% example Right-angle cantilever
%
% ------------------------------------------------------------------------------

clear all, close all
dirOnsas = [ pwd '/..' ] ;   problemName = 'rightAngleCantileverBeam' ;

% Material and geometrical properties-----------------------------------
E   =  1e6  ;
nu  = -0.5  ;
A   =  1    ;
I   =  1e-3 ;
L   =  10   ;
rho =  1    ;

% inconsistent
J   = I ;

materialsParams = {[ rho 1 E nu ]} ;

crossSecsParams = {[ 2 A I I J 20 10 10 ]} ;
							%Inertia dyadic tensor

%Node Mat
nElemsPerBeam = 10 ;
auxCoords = linspace(0,L,nElemsPerBeam+1)' ;

Nodes = [ zeros(nElemsPerBeam+1,1)       auxCoords  zeros(nElemsPerBeam+1,1) ; ...
          -auxCoords(2:end) ones(nElemsPerBeam,1)*L       zeros(nElemsPerBeam,1) ] ;

aux = (1:(2*nElemsPerBeam+1))' ;
          
%Conec Mat-------

%Nodes types 2 conec
					%M % E %L %C %S %NumNod	
fixedNod  = 		[0   1  0  0  1		1 ] ; % fixed node
loadedNod = 		[0   1  1  0  0 nElemsPerBeam+1 ] ; % loaded node


%Build Connec Element
Mat_Elem_MELCS = zeros(2*nElemsPerBeam,5);
for i=1:size(ones(2*nElemsPerBeam,5),1)
Mat_Elem_MELCS (i,:)  = [ 1 2 0 1 0] ;
end

Conec_elem = [ Mat_Elem_MELCS aux(1:(end-1)), aux(2:end) ] ;

Mat_NodeConec = cell(2+size(Conec_elem,1),1) ;

Conec{1, 1} = fixedNod ; % fixed node
Conec{2, 1} = loadedNod ; % loaded node
for i=1:size(Conec_elem,1)
  Conec{2+i, 1} =  Conec_elem(i,:) ;
end


elementsParams  = { 1; 3} ; %Nodo=1 Truss=1 Beam=3



% Method and tolerances	------------------------------------------------
% tolerances
stopTolDeltau = 0           ; 
stopTolForces = 1e-7        ;
stopTolIts    = 30          ;

timeIncr   =  0.25    ;
%~ finalTime  = .5    ;
%~ finalTime  = 3 ;    
finalTime  = 20 ;    


DeltaNW    =  0.5               ;
AlphaNW    =  0.25              ;
%~ numericalMethodParams = [ 3 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts DeltaNW AlphaNW] ;

%~ alphaHHT = 0 ;
alphaHHT = -0.05 ;
numericalMethodParams = [ 4 timeIncr finalTime stopTolDeltau stopTolForces stopTolIts alphaHHT ] ;



% ----------------------------------------------------------------------
% Loads and springs 	------------------------------------------------
%Loads
nLoadSteps = finalTime/timeIncr ;
                %Global %Variable
loadsParams   = {[ 1        1   0 0 0 0 0 1 0 ]} ;
loadFactorsFunc = @(t) 50*t*(t<1) + (100-50*t)*(t>=1)*(t<2) + 0 ;

%Springs
springsParams    = {[ inf  inf  inf  inf  inf  inf ]} ;






%POST PROCESS-----------------------------------------------------------
%Booleans store and plot
printFlag = 0 ;
reportBoolean = 0 ;
storeBoolean = 0 ;

%Plot parameters
controlDofs = [ nElemsPerBeam+1 3  1 ] ;
%~ plotParamsVector = [0 ];
plotParamsVector = [ 3 ]; sectPar = [ 12 .25 .25 ] ;


acdir = pwd ; cd(dirOnsas); ONSAS, cd(acdir) ;

lw  = 2   ; ms  = 5.5 ;
lw2 = 3.2 ; ms2 = 23 ;
plotfontsize = 22 ;

figure
%~ plot(controlDisps,'b--o','linewidth',lw,'markersize',ms);
plot(timesVec, controlDisps,'b--o','linewidth',lw,'markersize',ms);
grid on
labx = xlabel('Time (s)');   laby = ylabel(sprintf('Displacement node: %2i dof %1i', controlDofs(1), controlDofs(2) ) ) ;
set(gca, 'linewidth', 1.2, 'fontsize', plotfontsize )
set(labx, 'FontSize', plotfontsize); set(laby, 'FontSize', plotfontsize) ;

cd(dirOnsas); cd(outputDir);
%print('rightAngle','-dpdflatex','-tight')
%~ print('rightAngle','-dpdflatex')
%~ print('rightAngle.png','-dpng')
cd(acdir);


figure
grid on
plot(timesVec, loadFactors,'r','linewidth',lw,'markersize',ms);
