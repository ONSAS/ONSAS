% function for testing ONSAS using moxunit
% ----------------------------------------
function test_suite=runTestProblems_moxunit_disp
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  initTestSuite;

% function test_1
  % cd( [ 'examples' filesep ...
  %       'beamLinearVibration' ] );
  % beamLinearVibration
  % assertEqual( verifBoolean, true );

% function test_2
%  	cantileverModalAnalysis
%   assertEqual( verifBoolean, true );

% function test_3
%   addpath(genpath([ 'examples' filesep 'cantileverSelfWeight' ] ));
%   cantileverSelfWeight
%   assertEqual( verifBoolean, true );

% function test_4
%   addpath(genpath([ 'examples' filesep 'dragBeamReconfiguration' ] ));
%   dragBeamReconfiguration
%   assertEqual( verifBoolean, true );

% function test_5
% 	eulerColumn
%   assertEqual( verifBoolean, true );

% function test_6
%   frameLinearAnalysis
%   assertEqual( verifBoolean, true );

% function test_7
%   linearAerodynamics
%   assertEqual( verifBoolean, true );

% function test_8
%   addpath(genpath([ 'examples' filesep 'linearCylinderPlaneStrain' ] ) ) ;
%   linearCylinderPlaneStrain
%   assertEqual( verifBoolean, true );

% function test_9
%   nonlinearPendulum
%   assertEqual( verifBoolean, true );

% function test_10
%   springMass
%   assertEqual( verifBoolean, true );

% function test_11
%   simplePropeller
%   assertEqual( verifBoolean, true );

function test_12
  staticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_13
  uniaxialCompression
  assertEqual( verifBoolean, true );

function test_14
  addpath(genpath([ 'examples' filesep 'uniaxialExtension' ] ) ) ;
  uniaxialExtension
  assertEqual( verifBoolean, true);

% function test_15
%   cd( [ '..' filesep ...
%         'uniformCurvatureCantilever' ] );
%   uniformCurvatureCantilever
%   assertEqual( verifBoolean, true);

% function test_16
%   cd( [ '..' filesep ...
%         'VIVCantilever' ] );
%   VIVCantilever
%   assertEqual( verifBoolean, true );

function test_17
  assertEqual( gaussIntegrationTest, true);
