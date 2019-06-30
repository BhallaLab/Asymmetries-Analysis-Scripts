function output = sectionMapRead(f,point)

fseek(f,point,'bof');

a = fread(f,1,'int32');
b = fread(f,1,'int32');
c = fread(f,1,'long');

output = [a,b,c];
end