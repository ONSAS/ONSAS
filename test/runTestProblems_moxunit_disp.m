% function for testing ONSAS using moxunit
% ----------------------------------------
function test_suite=runTestProblems_moxunit_disp
  % initialize tests
  try
    test_functions=localfunctions()
  catch
  end

  % set auxiliar environment variable
  setenv('TESTS_RUN', 'yes') 

  % initialize the MOxUnit test suite
  initTestSuite;

function test_1
  beamLinearVibration
  assertEqual( verifBoolean, true );

function test_2
 	cantileverModalAnalysis
  assertEqual( verifBoolean, true );

function test_3
  cantileverSelfWeight
  assertEqual( verifBoolean, true );

function test_4
  dragBeamReconfiguration
  assertEqual( verifBoolean, true );

function test_5
	eulerColumn
  assertEqual( verifBoolean, true );

function test_6
  frameLinearAnalysis
  assertEqual( verifBoolean, true );

function test_7
  linearAerodynamics
  assertEqual( verifBoolean, true );

function test_8
  ringPlaneStrain
  assertEqual( verifBoolean, true );

function test_9
  nonlinearPendulum
  assertEqual( verifBoolean, true );

function test_10
  springMass
  assertEqual( verifBoolean, true );

function test_11
  simplePropeller
  assertEqual( verifBoolean, true );

function test_12
  staticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_13
  uniaxialCompression
  assertEqual( verifBoolean, true );

function test_14
  uniaxialExtension
  assertEqual( verifBoolean, true);

function test_15  
  uniformCurvatureCantilever
  assertEqual( verifBoolean, true);

function test_16
  VIVCantilever
  assertEqual( verifBoolean, true );

function test_17
  beamTrussJoint
  assertEqual( verifBoolean, true );

function test_18
  staticPlasticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_19
  platePatchTest
  assertEqual( verifBoolean, true );

function test_20
  cantileverPlate
  assertEqual( verifBoolean, true );

function test_21
  addedMassPendulum
  assertEqual( verifBoolean, true );
  
function test_22
  assertEqual( gaussIntegrationTest, true);
