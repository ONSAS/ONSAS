function [ aoa, cl, cd, cm ] = uBEMreadFoilData( polars )

% For uniform test
if iscell(polars) && numel(polars) >= 3
    format     = table2array(readtable( polars{3} ) );
    [ row, ~ ] = size(format)          ;
    for i = 1:length(polars)
        if i == 1 || i == 2   % First two columns are cylinder shapes
            aoa(:,i)   =    format(:,1) ;  
            for j = 1:row
                strucdata    = table2array( readtable( polars{i} ) )  ;
                cl(j,i)      = strucdata(1,2);
                cd(j,i)      = strucdata(1,3);
                cm(j,i)      = strucdata(1,4);
            end
        else
            strucdata       = table2array( readtable(  polars{i} ) );
            if length(strucdata(:,1)) ~= row
                [uniqueValues, ~, index] = unique(strucdata(:, 1));
                rowsToKeep               = accumarray(index, 1) == 1;
                filteredData             = strucdata(rowsToKeep, :) ;
                aoa(:,i)                 = format(:, 1) ;
                for j = 1:row
                    cl(j,i)    =  interp1( filteredData(:,1), filteredData(:,2),  format(j, 1), 'linear', 'extrap' );
                    cd(j,i)    =  interp1( filteredData(:,1), filteredData(:,3),  format(j, 1), 'linear', 'extrap' );
                    cm(j,i)    =  interp1( filteredData(:,1), filteredData(:,4),  format(j, 1), 'linear', 'extrap' );
                end
            else
                aoa(:,i)        = strucdata(:,1) ;
                cl(:,i)         = strucdata(:,2) ;
                cd(:,i)         = strucdata(:,3) ;
                cm(:,i)         = strucdata(:,4) ;
            end
        end
    end
else
     strucdata  = table2array( readtable( polars{1} ) )  ;
     aoa        = strucdata(:,1) ;
     cl         = strucdata(:,2) ;
     cd         = strucdata(:,3) ;
     cm         = strucdata(:,4) ;
end
end
