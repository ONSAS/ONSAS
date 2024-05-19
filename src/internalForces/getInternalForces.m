
function matFint = getInternalForces( forcesStruct, elems, fieldNames)

nElems  = length( elems      ) ;
nFields = length( fieldNames ) ;

matFint = zeros( nElems, nFields );

for i =1:nElems
    for j = 1:nFields
        getfield( forcesStruct(elems(i)), fieldNames{j} );
    end
end