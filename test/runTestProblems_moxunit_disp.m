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
  onsasExample_staticVonMisesTruss
  assertEqual( verifBoolean, true );

function test_2
  onsasExample_uniformCurvatureCantilever
  assertEqual( verifBoolean, true );

function test_3
  linearCylinderPlaneStrain
  assertEqual( verifBoolean, true );
  
function test_4
  uniaxialExtension
  assertEqual( verifBoolean, true );

function test_5
  nonlinearPendulum
  assertEqual( verifBoolean, true );

function test_6
  springMass
  assertEqual( verifBoolean, true );

function test_7
  onsasExample_cantileverSelfWeight
  assertEqual( verifBoolean, true );

function test_8
  simpleWindTurbine
  assertEqual( verifBoolean, true );

function test_9
  frameLinearAnalysis
  assertEqual( verifBoolean, true );

function test_10
  linearAerodynamics
  assertEqual( verifBoolean, true );

function test_11
  beamLinearVibration
  assertEqual( verifBoolean, true );

function test_12
  ONSAS_VIVtest
  assertEqual( verifBoolean, true );

function test_13
  assertEqual( gaussIntegrationTest, true);
  
function test_14
	cantileverLinearHardening
  assertEqual( verifBoolean, true);
