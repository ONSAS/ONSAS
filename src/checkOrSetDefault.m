function [outVar ] = checkOrSetDefault ( structName, fieldName, default )

if isfield( structName, fieldName )
  outVar = getfield( structName, fieldName ) ; 
else
  outVar = default ;
end
