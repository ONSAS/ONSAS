function fileName = auxReadDXF(nomArch)

%~ cd '../input'
cd input
fid = fopen(nomArch);
cd '../sources'
fileName = 'auxTxt.txt' ;

auxFile = fopen(fileName, 'w+');

line = fgetl(fid);
s={};
while ischar(line)
	s=[s;line];
	if length(line) == 0
		line = 'lala';
	else
		line = fgetl(fid);
	end
end
fclose(fid);

rows = size(s,1);
for i=1:rows
	fprintf(auxFile, '%s\n', s{i});
end

cd ..
fclose(auxFile);
