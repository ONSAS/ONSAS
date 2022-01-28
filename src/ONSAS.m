
function [ matUs, loadFactorsMat ] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams )

%mdFirst the input structs are converted to structs with the model information
[ modelCurrSol, modelProperties, BCsData ] = ONSAS_init( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

%mdAfter that the structs are used to perform the numerical time analysis
[ matUs, loadFactorsMat ] = ONSAS_solve( modelCurrSol, modelProperties, BCsData ) ;

%md Finally the report is generated
outputReport( modelProperties.outputDir, modelProperties.problemName )
