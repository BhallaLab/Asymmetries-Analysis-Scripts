function [fileFormat, numSweeps,channelCount,prots,data] = abfLoader(recFile)

fid = fopen(recFile1,'r');
if fid<0
    error('Something is wrong with the file')
end
fseek(fid,0,'bof');
fileFormat = fread(fid,4,'char=>char')';
    if fileFormat == "ABF2" || fileFormat == "ABF1"
        disp('')
    else
        error('Not a valid ABF file')
    end

fseek(fid,12,'bof');
numSweeps = fread(fid,1,'uint32');

%% Header (later)

%% Protocol Section
protSectionCounts = readStruct(fid,76);
fseek(fid,protSectionCounts(1)*512,'bof');
prots = fread(fid,512,'uint8');

%% Data Section
dataSectionCounts = readStruct(fid,236);
channelCount = readStruct(fid,92);channelCount = channelCount(3);
fseek(fid,dataSectionCounts(1)*512,'bof');
data = fread(fid,dataSectionCounts(3),'int16');

fclose(fid);
end

function output = readStruct(f,point)

fseek(f,point,'bof');

a = fread(f,1,'int32');
b = fread(f,1,'int32');
c = fread(f,1,'long');

output = [a,b,c];
end
