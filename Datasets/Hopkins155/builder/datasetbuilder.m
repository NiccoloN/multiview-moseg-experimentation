clear all;

cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('../'));

% List all files in the current directory
files = dir("../old");

% Initialize a cell array to store file names
SeqList = cell(1, numel(files));

% Loop through each file and store its name
for i = 1:numel(files)
    if ~files(i).isdir
        SeqList{i} = files(i).name;
    end
end

% Remove empty entries from the cell array
SeqList = SeqList(~cellfun('isempty', SeqList));

% Specify the file name to remove
fileToRemove = '.DS_Store'; 

% Find the index of the file to remove
idxToRemove = find(strcmp(SeqList, fileToRemove));

% Remove the file from the list if found
if ~isempty(idxToRemove)
    SeqList(idxToRemove) = [];
    disp(['File "', fileToRemove, '" removed from the list.']);
else
    disp(['File "', fileToRemove, '" not found in the list.']);
end

% Specify the file name to remove
fileToRemove = 'SeqList.mat'; 

% Find the index of the file to remove
idxToRemove = find(strcmp(SeqList, fileToRemove));

% Remove the file from the list if found
if ~isempty(idxToRemove)
    SeqList(idxToRemove) = [];
    disp(['File "', fileToRemove, '" removed from the list.']);
else
    disp(['File "', fileToRemove, '" not found in the list.']);
end

% Specify the file name to remove
fileToRemove = 'datasetbuilder.m'; 

% Find the index of the file to remove
idxToRemove = find(strcmp(SeqList, fileToRemove));

% Remove the file from the list if found
if ~isempty(idxToRemove)
    SeqList(idxToRemove) = [];
    disp(['File "', fileToRemove, '" removed from the list.']);
else
    disp(['File "', fileToRemove, '" not found in the list.']);
end

save("../SeqList", "SeqList");

% Iterate through each file in the list
for i = 1:numel(SeqList)
    % Load the data from the file
    fileName = SeqList{i};
    temp = load(fullfile('old/', fileName));
    
    Data.nFrames = temp.frames;
    Data.nSparsePoints = temp.points;
    Data.ySparse = temp.y;
    Data.GtLabel = temp.s;
    Data.visibleSparse = ones(temp.points, temp.frames);
    
    save(fullfile('../', fileName), 'Data');
end