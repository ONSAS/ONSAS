% function for testing ONSAS using moxunit
% ----------------------------------------
function test_suite=runTestProblems_moxunit_disp
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  setenv('TESTS_RUN', 'yes')
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

%function test_10
  % setenv('TESTS_RUN', 'yes')
  % springMass
 % assertEqual( verifBoolean, true );

%function test_11
  % addpath(genpath([ 'examples' filesep 'simplePropeller' ] ) ) ;
  % setenv('TESTS_RUN', 'yes')
  % simplePropeller
 % assertEqual( verifBoolean, true );

% function test_12
%   addpath(genpath([ 'examples' filesep 'staticVonMisesTruss' ] ) ) ;
%   staticVonMisesTruss
%   assertEqual( verifBoolean, true );

function test_13
  % setenv('TESTS_RUN', 'yes')
  uniaxialCompression
  assertEqual( verifBoolean, true );

% function test_14
%   addpath(genpath([ 'examples' filesep 'uniaxialExtension' ] ) ) ;
%   uniaxialExtension
%   assertEqual( verifBoolean, true);

function test_15  
  uniformCurvatureCantilever
  assertEqual( verifBoolean, true);

% function test_16
%   cd( [ '..' filesep ...
%         'VIVCantilever' ] );
%   VIVCantilever
%   assertEqual( verifBoolean, true );

function test_17
  assertEqual( gaussIntegrationTest, true);
