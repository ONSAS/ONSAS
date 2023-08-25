% Copyright 2023 Jorge M. Perez Zerpa, Joaquin Viera, Mauricio Vanzulli.
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

close all, clear all;
addpath( genpath( [ pwd filesep '..' filesep  'src' filesep ]) ); octaveBoolean = isThisOctave ;

setenv('TESTS_RUN','yes')

keyfiles = { ...
             'beamLinearVibration/beamLinearVibration.m'               ...
           ; 'beamTrussJoint/beamTrussJoint.m'                         ...
           ; 'cantileverModalAnalysis/cantileverModalAnalysis.m'       ...
           ; 'cantileverPlate/cantileverPlate.m'                       ...
           ; 'cantileverSelfWeight/cantileverSelfWeight.m'						 ...
           ; 'dragBeamReconfiguration/dragBeamReconfiguration.m'       ...
           ; 'eulerColumn/eulerColumn.m'                               ...
           ; 'frameLinearAnalysis/frameLinearAnalysis.m'               ...
           ; 'linearAerodynamics/linearAerodynamics.m'                 ...
           ; 'nonlinearPendulum/nonlinearPendulum.m'                   ...
           ; 'platePatchTest/platePatchTest.m'                         ...
           ; 'ringPlaneStrain/ringPlaneStrain.m'                       ...
           ; 'simplePropeller/simplePropeller.m'                       ...
           ; 'springMass/springMass.m'                                 ...
           ; 'staticVonMisesTruss/staticVonMisesTruss.m'               ...
           ; 'staticPlasticVonMisesTruss/staticPlasticVonMisesTruss.m' ...
           ; 'uniaxialCompression/uniaxialCompression.m'               ...
           ; 'uniaxialExtension/uniaxialExtension.m'                   ...
           ; 'uniformCurvatureCantilever/uniformCurvatureCantilever.m' ...
           ; 'VIVCantilever/VIVCantilever.m'                           ...
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
