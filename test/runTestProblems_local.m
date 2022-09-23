
close all, clear all;

if isunix, dirSep = '/'; else dirSep = '\'; end
addpath( [ pwd  dirSep '..' dirSep  'src' dirSep ] ); octaveBoolean = isThisOctave ;


keyfiles = { 'static_von_mises_truss/static_von_mises_truss.m'   ...
           ; 'uniformCurvatureCantilever/uniformCurvatureCantilever.m' ...
           ; 'linearCylinderPlaneStrain/linearCylinderPlaneStrain.m'    ...
           ; 'uniaxialExtension/uniaxialExtension.m'                    ...
           ; 'nonlinearPendulum/nonlinearPendulum.m'                    ...
           ; 'springMass/springMass.m'                                  ...
           ; 'cantileverSelfWeight/onsasExample_cantileverSelfWeight.m' ...
           ; 'simpleWindTurbine/simpleWindTurbine.m'                    ...
           ; 'frameLinearAnalysis/frameLinearAnalysis.m'                ...
           ; 'linearAerodynamics/linearAerodynamics.m'                  ...    
           ; 'reconfigurationBeam/circularReconfiguration.m'            ...    
           ; 'beamLinearVibration/beamLinearVibration.m'                ...
           ; 'VIVtest/ONSAS_VIVtest.m'                                  ...
           ; 'cantileverLinearHardening/cantileverLinearHardening.m' 		...   
           }

current  = 1 ;   verifBoolean = 1 ;  testDir = pwd ;

while current <= length(keyfiles) && verifBoolean == 1

  % run current example
  fprintf([' === running script: ' keyfiles{current} '\n' ]);

  % save key files data to avoid clear all commands
  save( '-mat', 'exData.mat', 'current', 'keyfiles', 'dirSep', 'testDir' );

  run( [ pwd dirSep '..' dirSep 'examples' dirSep keyfiles{current} ] ) ;

  if verifBoolean
    status = 'PASSED';
  else
    status = 'FAILED';
  end

  % reload key files data and increment current
  load('exData.mat') ;

  fprintf([' === test problem %2i:  %s === \n\n'], current, status );

  current = current + 1 ;
  delete('exData.mat');
  cd ( testDir )
end

if verifBoolean ==1
  fprintf('all test examples PASSED!\n')
else
  error('test examples not passed.')
end
