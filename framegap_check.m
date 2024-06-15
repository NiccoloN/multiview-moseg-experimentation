%% Preliminary operations
clear all
close all
clc
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('./'));

%% Set parameters
sel_dataset = 'KT3DMoSeg'; % KT3DMoSeg or Hopkins155
run_range = 1; % It can be either set as an array of the desired runs or to 'all'
                   % If it is not set to 'all', the error matrix will be
                   % computed in order to get arrays of best alpha, lambda
                   % and gamma values
frame_gaps = 1;  %Framegaps to be analyzed

%% Load remaining parameters
FG_Pars = []; % An array of structs containing parameters of the wanted frame gaps
for fg_ind = 1:length(frame_gaps)
    % Select runs
    if ischar(run_range) && strcmp(run_range, 'all')
        % Initialize an empty array to store the read numbers
        run_numbers = [];    
        % List all directories in the Runs folder
        run_folders = dir(fullfile('Runs', sel_dataset, frame_gaps(fg_ind), 'Run*'));    
        % Loop through each directory
        for i = 1:numel(run_folders)
            % Extract the number from the folder name (assuming it's in the format "RunX")
            run_number = str2double(extractAfter(run_folders(i).name, 'Run'));        
            % Append the number to the array
            run_numbers = [run_numbers, run_number];
        end
    elseif isnumeric(run_range)
        run_numbers = run_range;
    else
        % Handle invalid input
        error('Invalid input for run_range. It must be either a string or an array of numbers.');
    end
    
    % Load or evaluate parameters
    if ischar(run_range) && strcmp(run_range, 'all')
        % Load common parameters
        load(fullfile('Runs', sel_dataset, sprintf('Framegap%d', frame_gaps(fg_ind)), 'FG_Parameters.mat'));
    else
        % Read parameters of the first run
        load(fullfile('Runs', sel_dataset, sprintf('Framegap%d/Run%d', frame_gaps(fg_ind),...
            run_range(1)), 'Parameters.mat'), 'Par');
        prev_Par = Par;
        
        %Load parameters from the other runs and ensure they are equal
        for run_ind = run_range
            folder_name = sprintf('Framegap%d/Run%d', frame_gaps(fg_ind), run_ind);
            load(fullfile('Runs', sel_dataset, folder_name, 'Parameters.mat'), 'Par');
            if ~isequal(Par, prev_Par)
                error('Parameters are not consistent across runs for the selected frame gap. Problem encountered in Run%d', run_ind);
            end
            prev_Par = Par;
        end    
        seq_range = Par.seq_range;
        max_hypos = Par.max_hypos;
        alpha_range = Par.alpha_range;
        lambda_range = Par.lambda_range;
        gamma_range = Par.gamma_range;
    
        % Calculate Error_matrix
        Error_matrix = error_analisys(sel_dataset, run_range, seq_range, max_hypos, alpha_range, lambda_range, gamma_range, frame_gaps(fg_ind));
        set(0, 'DefaultFigureVisible', 'off'); %Deactivate figure visibility
        sel_lambda = lambda_graph(Error_matrix, alpha_range, lambda_range);
        sel_gamma = gamma_graph(Error_matrix, alpha_range, gamma_range);
        sel_alphas = alpha_graph(Error_matrix, alpha_range, lambda_range, gamma_range, sel_lambda, sel_gamma);
        set(0, 'DefaultFigureVisible', 'on'); %Reactivate figure visibility
    
        %Populate FG_Par struct
        FG_Par = Par;
        FG_Par.run_range = run_range;
        FG_Par.Error_matrix = Error_matrix;
        FG_Par.best_lambda = sel_lambda;
        FG_Par.best_gamma = sel_gamma;
        FG_Par.best_alphas = sel_alphas;
    end
    FG_Pars = [FG_Pars FG_Par];
end

%% Sequence by sequence graph for each frame gap
for fg_ind = 1:length(frame_gaps)
    seq_graph(FG_Pars(fg_ind).Error_matrix, FG_Pars(fg_ind).alpha_range,...
        FG_Pars(fg_ind).lambda_range, FG_Pars(fg_ind).gamma_range,...
        FG_Pars(fg_ind).best_alphas, FG_Pars(fg_ind).best_lambda,...
        FG_Pars(fg_ind).best_gamma, FG_Pars(fg_ind).seq_range)
end

%% Sequence by sequence frame gap comparison
best_seq_graph(FG_Pars);
