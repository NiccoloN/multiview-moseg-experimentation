%% Preliminary operations
clear
close all
clc
warning('off', 'all')
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('./'));

calc_matrix = false;

%% Set parameters
[sel_dataset, ork_method, clustering_method, run_range, frame_gap_range, ork_philosophy, sel_seqs] = select_parameters();
display_graphs = true; %Turn it to false if you do not want to display figures for alpha, lambda and gamma graphs

run_parameters = struct();
run_parameters.sel_dataset = sel_dataset;
run_parameters.ork_method = ork_method;
run_parameters.clustering_method = clustering_method;
run_parameters.frame_gap_range = frame_gap_range;
run_parameters.ork_philosophy = ork_philosophy;

%% General checks on selected parameters
% Detect if there are many simulations differing only for the selected sequences and act accordingly
subfolders = dir(fullfile('Runs', sel_dataset));
subfolders = subfolders([subfolders.isdir]);
subfolder_names = {subfolders.name};
similar_subfolders = subfolder_names(startsWith(subfolder_names, sprintf('ORK%s-%s_CLUSTERING%s_FRAMEGAP%d_SEQS', ork_philosophy, ork_method, clustering_method, frame_gap_range)));
switch numel(similar_subfolders) 
    case 0
        error('There is no simulation performed with the selected parameters. Please run one and retry')
    case 1
        main_folder_name = fullfile('Runs', sel_dataset, similar_subfolders{1});
        path_name = fullfile(main_folder_name, 'Perf_Parameters.mat');
    otherwise
        error('There is more than one simulation with the same parameters, but evaluated on a different sequence range. Please remove the unwanted one and retry')
end

% Check if the selected run_range is compatible with the number of present runs
run_folders = dir(fullfile(main_folder_name, 'Run*'));
run_numbers = [];
for i = 1:numel(run_folders)
    folder_name = run_folders(i).name;
    run_number = sscanf(folder_name, 'Run%d');
    run_numbers = [run_numbers, run_number];
end
if isequal(run_range, run_numbers)
    disp('All the performed runs are being considered for the performance check.');
elseif ismember(run_range, run_numbers)
    not_used_runs = setdiff(run_numbers, run_range);
    disp('Not all runs are being used for the performance check.');
    disp('The following runs are not being used:');
    disp(not_used_runs);
else
    not_performed_runs = setdiff(run_range, run_numbers);
    error('Not all the selected runs have been performed.\nThe following runs have not been performed:\n%s\nThe following runs can be analyzed:\n%s', num2str(not_performed_runs), num2str(run_numbers));
end

if display_graphs
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% Error matrix and graphs generation
% Check if data already exists
if exist(path_name, 'file') == 2
    % Load the existing struct from the file
    temp = load(path_name, 'Perf_Parameters');
    Perf_Parameters = temp.Perf_Parameters;
    clear temp
    calc_matrix = false;
else
    % There is no existing data
    disp('No existing data has been found. Evaluating and saving new data');
    calc_matrix = true;
end

if calc_matrix
    disp('Evaluating new Error_matrix...')

    % Read parameters from the runs and ensure they are equal
    load(fullfile(main_folder_name, sprintf('Run%d', run_range(1)), 'Parameters.mat'), 'Par');
    prev_Par = Par;
    for run_ind = run_range
        folder_name = fullfile(main_folder_name, sprintf('Run%d', run_ind));
        load(fullfile(folder_name, 'Parameters.mat'), 'Par');
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
    clear prev_Par run_ind
    
    % Calculate Error_matrix
    [Error_matrix, Clusters_matrix] = error_analisys(main_folder_name, run_range, seq_range, max_hypos, alpha_range, lambda_range, gamma_range);
    disp('Done');

    % Lambda-alpha graph
    sel_lambda = lambda_range(1);
    if length(lambda_range) > 1
        sel_lambda = lambda_graph(Error_matrix, alpha_range, lambda_range);
    end
    
    % Gamma-alpha graph
    sel_gamma = gamma_range(1);
    if length(gamma_range) > 1
        sel_gamma = gamma_graph(Error_matrix, alpha_range, gamma_range);
    end
    
    % Alpha graph
    sel_alphas = alpha_graph(Error_matrix, alpha_range, lambda_range, gamma_range, sel_lambda, sel_gamma, sel_seqs, run_parameters);
    
    % Save data
    Perf_Parameters = Par;
    Perf_Parameters.run_range = run_range;
    Perf_Parameters.Error_matrix = Error_matrix;
    Perf_Parameters.Clusters_matrix = Clusters_matrix;
    Perf_Parameters.best_lambda = sel_lambda;
    Perf_Parameters.best_gamma = sel_gamma;
    Perf_Parameters.best_alphas = sel_alphas;

    save(path_name, 'Perf_Parameters');
    disp('New parameters saved successfully.');

else
    % Load saved data
    % Get the field names of Perf_Parameters
    field_names = fieldnames(Perf_Parameters);    
    % Loop through each field and assign it to a variable with the same name in the workspace
    for i = 1:length(field_names)
        assignin('base', field_names{i}, Perf_Parameters.(field_names{i}));
    end
    % Display graphs
    sel_lambda = lambda_range(1);
    sel_gamma = gamma_range(1);
    if length(lambda_range) > 1
        sel_lambda = lambda_graph(Error_matrix, alpha_range, lambda_range);
    end
    if length(gamma_range) > 1
        sel_gamma = gamma_graph(Error_matrix, alpha_range, gamma_range);
    end
    sel_alphas = alpha_graph(Error_matrix, alpha_range, lambda_range, gamma_range, sel_lambda, sel_gamma, sel_seqs, run_parameters);
end

set(0, 'DefaultFigureVisible', 'on'); %Reactivate figure visibility

%% Sequence by sequence graphs
seq_graph(Error_matrix, Clusters_matrix, alpha_range, lambda_range, gamma_range, sel_alphas, sel_lambda, sel_gamma, sel_seqs, run_parameters)
seq_bestalpha_graph(Error_matrix, Clusters_matrix, alpha_range, lambda_range, gamma_range, sel_lambda, sel_gamma, sel_seqs, run_parameters)
% lambda_graph(Error_matrix, alpha_range, lambda_range, sel_seqs)
% gamma_graph(Error_matrix, alpha_range, gamma_range, sel_seqs)

%% Find best alphas for view_points
[Best_alphas_matrix, Best_alpha_ind_matrix] = find_best_alphas_matrix(Error_matrix, alpha_range, lambda_range, gamma_range, sel_seqs, sel_lambda, sel_gamma);

%% Helper functions

 function [sel_dataset, ork_method, clustering_method, run_range, frame_gap_range, ork_philosophy, sel_seqs] = select_parameters()

    % Define options for each parameter
    dataset_options = {'KT3DMoSeg', 'Hopkins155'};
    ork_options = {'defaultORK', 'altORK1', 'altORK2', 'altORK3'};
    clustering_options = {'kmeans', 'dbscan', 'linkage'};
    ork_philosophy_options = {'simplesum', 'multiframe','multiframe_sum','multiframe_score','single_frame_gap'};
    
    % Get screen size
    screen_size = get(0, 'ScreenSize');
    screen_width = screen_size(3);
    screen_height = screen_size(4);
    
    % Calculate position of the figure to center it
    figure_width = 600;
    figure_height = 620;
    figure_position = [(screen_width - figure_width) / 2, (screen_height - figure_height) / 2, figure_width, figure_height];
    
    % Create a figure
    f = figure('Name', 'Select Parameters', 'Position', figure_position, 'NumberTitle', 'off');
    
    % Calculate width for each listbox
    listbox_width = (figure_width - 150) / 3; % Subtracting 50 for padding and dividing by 3 for three listboxes
    
    % Create listboxes for each parameter
    lb_dataset = uicontrol('Style', 'listbox', 'String', dataset_options, 'Position', [50, 480, listbox_width, 100]);
    lb_ork = uicontrol('Style', 'listbox', 'String', ork_options, 'Position', [50 + listbox_width + 25, 480, listbox_width, 100]);
    lb_clustering = uicontrol('Style', 'listbox', 'String', clustering_options, 'Position', [50 + 2*(listbox_width + 25), 480, listbox_width, 100]);
    
    % Add titles to each listbox
    title_dataset = uicontrol('Style', 'text', 'String', 'Select Dataset:', 'Position', [50, 580, listbox_width, 20], 'HorizontalAlignment', 'center');
    title_ork = uicontrol('Style', 'text', 'String', 'ORK Method:', 'Position', [50 + listbox_width + 25, 580, listbox_width, 20], 'HorizontalAlignment', 'center');
    title_clustering = uicontrol('Style', 'text', 'String', 'Clustering Method:', 'Position', [50 + 2*(listbox_width + 25), 580, listbox_width, 20], 'HorizontalAlignment', 'center');
    
    % Create a listbox for selecting ork_philosophy
    lb_ork_philosophy = uicontrol('Style', 'listbox', 'String', ork_philosophy_options, 'Position', [50, 400, 500, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for ork_philosophy
    title_ork_philosophy = uicontrol('Style', 'text', 'String', 'Select ORK Philosophy:', 'Position', [50, 440, 500, 20], 'HorizontalAlignment', 'left');
    
    % Create an edit box for run range
    edit_run_range = uicontrol('Style', 'edit', 'String', '1:3', 'Position', [50, 320, 500, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for run range
    title_run_range = uicontrol('Style', 'text', 'String', 'Enter a run range:', 'Position', [50, 360, 500, 20], 'HorizontalAlignment', 'left');
    
    % Create an edit box for frame gap range
    edit_frame_gap_range = uicontrol('Style', 'edit', 'Position', [50, 240, 500, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for frame gap range
    title_frame_gap_range = uicontrol('Style', 'text', 'String', 'Select the frame gap range (just one integer number):', 'Position', [50, 280, 500, 20], 'HorizontalAlignment', 'left');
    
    % Create an edit box for selecting sequences
    edit_sel_seqs = uicontrol('Style', 'edit', 'String', 'all', 'Position', [50, 160, 500, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for selecting sequences
    title_sel_seqs = uicontrol('Style', 'text', 'String', 'Select sequences (e.g., "1:3" or "all"):', 'Position', [50, 200, 500, 20], 'HorizontalAlignment', 'left');
    
    % Create an OK button
    btn_ok = uicontrol('Style', 'pushbutton', 'String', 'OK', 'Position', [220, 40, 120, 30], ...
                       'Callback', @okCallback);
    
    % Initialize variables to store selections
    sel_dataset = '';
    ork_method = '';
    clustering_method = '';
    run_range = 1:3;
    frame_gap_range = '';
    ork_philosophy = '';
    sel_seqs = '';
    
    % Callback function for OK button
    function okCallback(~, ~)
        sel_dataset = dataset_options{get(lb_dataset, 'Value')};
        ork_method = ork_options{get(lb_ork, 'Value')};
        clustering_method = clustering_options{get(lb_clustering, 'Value')};
        
        % Parse the string entered for run range
        run_range_str = get(edit_run_range, 'String');
        if contains(run_range_str, ':') % Check if a range is entered
            run_range_parts = strsplit(run_range_str, ':');
            run_range = str2double(run_range_parts{1}):str2double(run_range_parts{2});
        else
            run_range = str2double(run_range_str); % Convert single number directly
        end
        
        % Parse the string entered for frame gap range
        frame_gap_range = str2double(get(edit_frame_gap_range, 'String'));
        
        % Get the selected ORK philosophy
        ork_philosophy = ork_philosophy_options{get(lb_ork_philosophy, 'Value')};
        
        % Parse the string entered for selected sequences
        sel_seqs_str = get(edit_sel_seqs, 'String');
        if strcmpi(sel_seqs_str, 'all')
            if strcmp(sel_dataset, 'KT3DMoSeg')
                sel_seqs = 1:22;
            elseif strcmp(sel_dataset, 'Hopkins155')
                sel_seqs = 1:156;
            else
                error('Invalid dataset selection. sel_seqs can only be "all" for KT3DMoSeg or Hopkins155 datasets.');
            end
        elseif contains(sel_seqs_str, ':')
            sel_seqs_parts = strsplit(sel_seqs_str, ':');
            sel_seqs = str2double(sel_seqs_parts{1}):str2double(sel_seqs_parts{2});
        elseif ~isempty(regexp(sel_seqs_str, '^\d+$', 'once'))
            sel_seqs = str2double(sel_seqs_str);
        else
            error('Invalid sequence values selection. sel_seqs can only be "all", an array declared as start:end, or an integer.');
        end
        
        % Check if sel_seqs is within the valid range for selected dataset
        if strcmp(sel_dataset, 'KT3DMoSeg') && (min(sel_seqs) < 1 || max(sel_seqs) > 22)
            error('Invalid sequence selection. sel_seqs must be within the range 1:22 for KT3DMoSeg dataset.');
        elseif strcmp(sel_dataset, 'Hopkins155') && (min(sel_seqs) < 1 || max(sel_seqs) > 156)
            error('Invalid sequence selection. sel_seqs must be within the range 1:156 for Hopkins155 dataset.');
        end

        % Close the figure
        close(f);
    end

    % Wait for the figure to close
    uiwait(f);
end