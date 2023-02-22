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

function test_2
%   cd( [ '..' filesep ...
%         'cantileverModalAnalysis' ] );
 	cantileverModalAnalysis
  assertEqual( verifBoolean, true );

function test_3
  addpath(genpath([ 'examples' filesep 'cantileverSelfWeight' ] ));
  cantileverSelfWeight
  assertEqual( verifBoolean, true );

% function test_4
%   cd( [ '..' filesep ...
%         'dragBeamReconfiguration' ] );
%   dragBeamReconfiguration
%   assertEqual( verifBoolean, true );

% function test_5
%   cd( [ '..' filesep ...
%         'eulerColumn' ] );
% 	eulerColumn
%   assertEqual( verifBoolean, true );

% function test_6
%   cd( [ '..' filesep ...
%         'frameLinearAnalysis' ] );
%   frameLinearAnalysis
%   assertEqual( verifBoolean, true );

% function test_7
%   cd( [ '..' filesep ...
%         'linearAerodynamics' ] );
%   linearAerodynamics
%   assertEqual( verifBoolean, true );

% function test_8
%   cd( [ '..' filesep ...
%         'linearCylinderPlaneStrain' ] );
%   linearCylinderPlaneStrain
%   assertEqual( verifBoolean, true );

% function test_9
%   cd( [ '..' filesep ...
%         'nonlinearPendulum' ] );
%   nonlinearPendulum
%   assertEqual( verifBoolean, true );

% function test_10
%   cd( [ '..' filesep ...
%         'springMass' ] );
%   springMass
%   assertEqual( verifBoolean, true );

% function test_11
%   cd( [ '..' filesep ...
%         'simplePropeller' ] );
%   simplePropeller
%   assertEqual( verifBoolean, true );

% function test_12
%   cd( [ '..' filesep ...
%         'staticVonMisesTruss' ] );
%   staticVonMisesTruss
%   assertEqual( verifBoolean, true );

function test_13
%   cd( [ '..' filesep ...
%         'uniaxialCompression' ] );
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
  % cd( [ '..' filesep '..' filesep 'test'] );
  assertEqual( gaussIntegrationTest, true);
  % cd( [ '..' ] );
