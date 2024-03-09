function [ aoa, cl, cd, cm ] = BEMreadFoilData( path, polars )

AoAlist       = (-180:5:180 )';

for i = 1:length(polars)
    if i == 1 || i == 2   % First two columns are cylinder shapes
        aoa(:,i)   =    AoAlist(:,1) ;  
        for j = 1:length(AoAlist)
            strucdata    = table2array( readtable( fullfile( path, polars{i} ) ) )  ;
            cl(j,i)      = strucdata(1,2);
            cd(j,i)      = strucdata(1,3);
            cm(j,i)      = strucdata(1,4);
        end
    else
        strucdata       = table2array( readtable(  polars{i} ) );
        if length(strucdata(:,1)) ~= length(AoAlist)
            [uniqueValues, ~, index] = unique(strucdata(:, 1));
            rowsToKeep               = accumarray(index, 1) == 1;
            filteredData             = strucdata(rowsToKeep, :) ;
            aoa(:,i)                 = AoAlist(:, 1) ;
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
end
