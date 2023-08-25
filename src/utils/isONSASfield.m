function isONSASfield(ONSASstruct, structName, chars)

    fields = fieldnames(ONSASstruct) ;
    len = length(chars) ; % as cell
    for i = 1:len
        chars{i}
        fields{i}
        if strcmp(chars{i},fields{i}) ~= 1
            error("The struct %s does not contain the struct %s, to see the fields please check: https://onsas.github.io/ONSAS.m/dev/howtouse/creatingModels/", structName, fields{i}) 
        end     
    end

end