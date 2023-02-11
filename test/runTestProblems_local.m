
close all, clear all;
addpath( [ pwd filesep '..' filesep  'src' filesep ] ); octaveBoolean = isThisOctave ;

keyfiles = { ...
             'beamLinearVibration/beamLinearVibration.m'               ...
           ; 'cantileverSelfWeight/onsasExample_cantileverSelfWeight.m'...
           ; 'cantileverModalAnalysis/cantileverModalAnalysis.m'   ...
           ; 'cantileverLinearHardening/cantileverLinearHardening.m' 	 ...
           ; 'eulerColumn/eulerColumn.m'                               ...
           ; 'frameLinearAnalysis/frameLinearAnalysis.m'               ...
           ; 'linearAerodynamics/linearAerodynamics.m'                 ...
           ; 'linearCylinderPlaneStrain/linearCylinderPlaneStrain.m'   ...
           ; 'uniformCurvatureCantilever/uniformCurvatureCantilever.m' ...
           ; 'uniaxialExtension/uniaxialExtension.m'                   ...
           ; 'nonlinearPendulum/nonlinearPendulum.m'                   ...
           ; 'reconfigurationBeam/circularReconfiguration.m'           ...
           ; 'springMass/springMass.m'                                 ...
           ; 'simplePropeller/simplePropeller.m'                   ...
           ; 'static_von_mises_truss/static_von_mises_truss.m'   ...
           ; 'uniaxialCompression/uniaxialCompression.m'               ...
           ; 'VIVtest/ONSAS_VIVtest.m'                                 ...
           }

current  = 1 ;   verifBoolean = 1 ;  testDir = pwd ;

num_tests = length(keyfiles) ;
while (current <= num_tests) && (verifBoolean == 1)

  % run current example
  fprintf([' === running script: ' keyfiles{current} '\n' ]);

  aux_time = cputime();

  % save key files data to avoid clear all commands
  save( '-mat', 'exData.mat', 'current', 'keyfiles', 'testDir', 'aux_time' );

  run( [ pwd filesep '..' filesep 'examples' filesep keyfiles{current} ] ) ;

  if verifBoolean
    status = 'PASSED';
  else
    status = 'FAILED';
  end

  % reload key files data and increment current
  load('exData.mat') ; num_tests = length(keyfiles) ;

  aux_time = cputime() - aux_time ; keyfiles{current,2} = aux_time ;

  fprintf([' === test problem %2i:  %s in %8.1e s === \n\n'], current, status, aux_time );

  current = current + 1 ;
  delete('exData.mat');
  cd ( testDir )
end

if verifBoolean ==1
  fprintf('all test examples PASSED!\n')
else
  error('test examples not passed.')
end
