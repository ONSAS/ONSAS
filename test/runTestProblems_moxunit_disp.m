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
  beamLinearVibration
  assertEqual( verifBoolean, true );

function test_2
	cantilever_ModalAnalysis
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
  linearCylinderPlaneStrain
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

function test_18
  assertEqual( gaussIntegrationTest, true);
