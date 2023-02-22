% function for testing ONSAS using moxunit
% ----------------------------------------
function test_suite=runTestProblems_moxunit_disp
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end
  initTestSuite;

function test_1
  cd( [ 'examples' filesep ...
        'beamLinearVibration' ] );
  beamLinearVibration
  assertEqual( verifBoolean, true );

function test_2
pwd
  cd( [ '..' filesep ...
        'cantileverModalAnalysis' ] );
	cantileverModalAnalysis
  assertEqual( verifBoolean, true );

function test_3
  cd( [ 'examples' filesep ...
        'cantileverSelfWeight' ] );
  cantileverSelfWeight
  assertEqual( verifBoolean, true );

function test_4
  cd( [ 'examples' filesep ...
        'dragBeamReconfiguration' ] );
  dragBeamReconfiguration
  assertEqual( verifBoolean, true );

function test_5
  cd( [ 'examples' filesep ...
        'eulerColumn' ] );
	eulerColumn
  assertEqual( verifBoolean, true );

function test_6
  cd( [ 'examples' filesep ...
        'frameLinearAnalysis' ] );
  frameLinearAnalysis
  assertEqual( verifBoolean, true );

function test_7
  cd( [ 'examples' filesep ...
        'linearAerodynamics' ] );
  linearAerodynamics
  assertEqual( verifBoolean, true );

function test_8
  cd( [ 'examples' filesep ...
        'linearCylinderPlaneStrain' ] );
  linearCylinderPlaneStrain
  assertEqual( verifBoolean, true );

function test_9
  cd( [ 'examples' filesep ...
        'nonlinearPendulum' ] );
  nonlinearPendulum
  assertEqual( verifBoolean, true );

function test_10
  cd( [ 'examples' filesep ...
        'springMass' ] );
  springMass
  assertEqual( verifBoolean, true );

function test_11
  cd( [ 'examples' filesep ...
        'simplePropeller' ] );
  simplePropeller
  assertEqual( verifBoolean, true );

function test_12
  cd( [ 'examples' filesep ...
        'staticVonMisesTruss' ] );
  staticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_13
  cd( [ 'examples' filesep ...
        'uniaxialCompression' ] );
  uniaxialCompression
  assertEqual( verifBoolean, true );

function test_14
  cd( [ 'examples' filesep ...
        'uniaxialExtension' ] );
  uniaxialExtension
  assertEqual( verifBoolean, true);

function test_15
  cd( [ 'examples' filesep ...
        'uniformCurvatureCantilever' ] );
  uniformCurvatureCantilever
  assertEqual( verifBoolean, true);

function test_16
  cd( [ 'examples' filesep ...
        'VIVCantilever' ] );
  VIVCantilever
  assertEqual( verifBoolean, true );

function test_17
  assertEqual( gaussIntegrationTest, true);
