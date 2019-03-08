clearvars
clc

%% To set the default directories

importDir = uigetdir(); %Opens up a folder selection box
exportDir = importDir;

%% read from the import directory and process
cd(importDir) %Make it the current working directory
allFiles = dir('Expt02*.atf'); %list all the files 


for i=1:length(allFiles)
    %Fetch the file
    TraceFile = allFiles(i).name;
    ExptID = strsplit(TraceFile,'.'); %Separate file name from extension
    ExptID = ExptID{1}; % This variable is going to help create files later
    
    %The data in Tracefiles is tab deliminted and starts at row#12 and column#1
    data = dlmread(TraceFile,'\t',11,0);
    
    %Load the file
    %To check the size of the the variable 'data' inside each file
    %For consistency every file has a single variable called 'data'
    RecSize = size(data);
    
    % Generate a single trace for Ch#1 and Ch#2
    PatchTrace = [];
    PolygonTrace = [];

    %PatchTrace is every even column starting from 2
    %PolygonTrace is every odd column starting from 2
    for k=2:RecSize(2)
        if mod(k,2)==0
            PatchTrace = [PatchTrace;data(:,k)];
        else
            PolygonTrace = [PolygonTrace;data(:,k)];
        end
    end
        
    %Generate the file containing the two required variables
    exportFileName = strcat(exportDir,'\',ExptID);
    save(exportFileName, 'P*Trace')
    
    clear exportFileName
    clear k
    clear PathName
    clear TraceFile
    clear RecSize
    clear data
end
%---------------------------------------------------------------------