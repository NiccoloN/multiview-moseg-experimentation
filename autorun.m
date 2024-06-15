%% Autorun script
clc
clear all
% close all
warning('off', 'all')

gcp; % creates parallel pool if it doesn't exist
tstart = tic;

%% Set path
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('./'));

%% Set parameters
% Set general parameters
[sel_dataset, ork_method, clustering_method, reuse_kernels, run_range, frame_gaps, ork_philosophy] = select_parameters();

% Internal parameters:
use_default_parameters = true; % Set to true to force parameters to match ones used in the other runs with the same framegap
                               % (suggested). This way, this run results can be aggregated with results of the other ones
                               % for a cumulative analisys

% Manually set internal parameters NOT SUGGESTED
seq_range = 1;       % 1:22 for K3DMoSeg, 1:156 for Hopkins155
max_hypos = 500;
alpha_range = 5:15;
% lambda_range = [0 0.002 0.004 0.006 0.008 0.01 0.02 0.03 0.04];
% gamma_range = [0 0.002 0.004 0.006 0.008 0.01 0.02 0.03 0.04];
lambda_range = 0.01;
gamma_range = 0.01;

% Or load default parameters SUGGESTED
if use_default_parameters
    [seq_range, max_hypos, alpha_range, lambda_range, gamma_range] = load_parameters(sel_dataset);
end

% Generate parameters struct
Par = struct();
Par.seq_range = seq_range;
Par.max_hypos = max_hypos;
Par.alpha_range = alpha_range;
Par.lambda_range = lambda_range;
Par.gamma_range = gamma_range;

%% Load data
temp = load(fullfile('Datasets/', sel_dataset, '/SeqList.mat'));
seq_list = temp.SeqList;

data = cell(1, length(seq_list));
for s_i = 1:length(seq_list)
    seq_name = seq_list{s_i};
    filepath = fullfile('Datasets/', sel_dataset,'/', seq_name);
    temp = load(filepath);
    data{s_i} = temp.Data;
end

%% Motion segmentation with random sampling
for r_i = run_range
    for frame_gap = frame_gaps
        clc
        fprintf('Run %d, with frame gap %d in progress...\n\n', r_i, frame_gap);
    
        %Clear old results
        if isfolder('Results')
            rmdir('Results','s')
        end
    
        %Make run folder
        folder_name = sprintf('ORK%s-%s_CLUSTERING%s_FRAMEGAP%d_SEQS%d-%d/Run%d', ork_philosophy, ork_method, clustering_method, frame_gap, seq_range(1), seq_range(end), r_i);
        target_folder = fullfile('Runs/', sel_dataset, folder_name);
    
        if not(isfolder(target_folder))
            mkdir(target_folder);
        end
        
        %Save parameters struct
        save(fullfile(target_folder, 'Parameters.mat'), 'Par');
    
        %Save ground truth for ML approach
        for s_i = seq_range
    
            seq_name = seq_list{s_i};
            target_folder = fullfile('Results/ML',seq_name);
            
            if not(isfolder(target_folder))
                mkdir(target_folder);
            end
    
            gt = data{s_i}.GtLabel;
            save(fullfile(target_folder, 'GroundTruth'), 'gt');
        end
    
        if reuse_kernels ~= 1
            switch ork_philosophy
                case 'multiframe_score'
                    frame_gap_range = 1:frame_gap;

                    %Hypothesis generation
                    fprintf('Affine Hypotheses:\n');
                    hypos_A = hypo_conspts('affine', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Affine Hypothesis: Done\n\n');
                
                    fprintf('Homography Hypotheses:\n');
                    hypos_H = hypo_conspts('homography', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Homography Hypothesis: Done\n\n');

                    fprintf('Fundamental Hypotheses:\n');
                    hypos_F = hypo_conspts('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Fundamental Hypothesis: Done\n\n');
                
                    %Applying ORK        
                    fprintf('Affine ORK:\n');
                    kernels_A = ork_multiframe_score('affine', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_A, ork_method);
                    fprintf('Affine ORK: Done\n\n');
                
                    fprintf('Homography ORK:\n');
                    kernels_H = ork_multiframe_score('homography', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_H, ork_method);
                    fprintf('Homography ORK: Done\n\n');
                
                    fprintf('Fundamental ORK:\n');
                    kernels_F = ork_multiframe_score('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_F, ork_method);
                    fprintf('Fundamental ORK: Done\n\n');

                case 'multiframe_sum'
                    frame_gap_range = 1:frame_gap;

                    %Hypothesis generation
                    fprintf('Affine Hypotheses:\n');
                    hypos_A = hypo_conspts('affine', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Affine Hypothesis: Done\n\n');
                
                    fprintf('Homography Hypotheses:\n');
                    hypos_H = hypo_conspts('homography', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Homography Hypothesis: Done\n\n');
            
                    fprintf('Fundamental Hypotheses:\n');
                    hypos_F = hypo_conspts('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Fundamental Hypothesis: Done\n\n');
                
                    %Applying ORK        
                    fprintf('Affine ORK:\n');
                    kernels_A = ork_simplesum('affine', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_A, ork_method);
                    fprintf('Affine ORK: Done\n\n');
                
                    fprintf('Homography ORK:\n');
                    kernels_H = ork_simplesum('homography', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_H, ork_method);
                    fprintf('Homography ORK: Done\n\n');
                
                    fprintf('Fundamental ORK:\n');
                    kernels_F = ork_simplesum('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_F, ork_method);
                    fprintf('Fundamental ORK: Done\n\n');
                case 'multiframe'
                    frame_gap_range = 1:frame_gap;

                    %Hypothesis generation
                    fprintf('Affine Hypotheses:\n');
                    hypos_A = hypo_conspts('affine', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Affine Hypothesis: Done\n\n');
                
                    fprintf('Homography Hypotheses:\n');
                    hypos_H = hypo_conspts('homography', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Homography Hypothesis: Done\n\n');
            
                    fprintf('Fundamental Hypotheses:\n');
                    hypos_F = hypo_conspts('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Fundamental Hypothesis: Done\n\n');
                
                    %Applying ORK        
                    fprintf('Affine ORK:\n');
                    kernels_A = ork_multiframe('affine', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_A, ork_method);
                    fprintf('Affine ORK: Done\n\n');
                
                    fprintf('Homography ORK:\n');
                    kernels_H = ork_multiframe('homography', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_H, ork_method);
                    fprintf('Homography ORK: Done\n\n');
                
                    fprintf('Fundamental ORK:\n');
                    kernels_F = ork_multiframe('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_F, ork_method);
                    fprintf('Fundamental ORK: Done\n\n');
                
                case 'simplesum'
                    frame_gap_range = 1:frame_gap;

                    %Hypothesis generation
                    fprintf('Affine Hypothesis:\n');
                    hypos_A = hypo_anypts('affine', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Affine Hypothesis: Done\n\n');
                
                    fprintf('Homography Hypothesis:\n');
                    hypos_H = hypo_anypts('homography', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Homography Hypothesis: Done\n\n');

                    fprintf('Fundamental Hypothesis:\n');
                    hypos_F = hypo_anypts('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Fundamental Hypothesis: Done\n\n');
                
                    %Applying ORK            
                    fprintf('Affine ORK:\n');
                    kernels_A = ork_simplesum('affine', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_A, ork_method);
                    fprintf('Affine ORK: Done\n\n');
                
                    fprintf('Homography ORK:\n');
                    kernels_H = ork_simplesum('homography', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_H, ork_method);
                    fprintf('Homography ORK: Done\n\n');
                
                    fprintf('Fundamental ORK:\n');
                    kernels_F = ork_simplesum('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_F, ork_method);
                    fprintf('Fundamental ORK: Done\n\n');

                case 'single_frame_gap'
                    frame_gap_range = frame_gap;

                    %Hypothesis generation
                    fprintf('Affine Hypothesis:\n');
                    hypos_A = hypo_anypts('affine', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Affine Hypothesis: Done\n\n');
                
                    fprintf('Homography Hypothesis:\n');
                    hypos_H = hypo_anypts('homography', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Homography Hypothesis: Done\n\n');
            
                    fprintf('Fundamental Hypothesis:\n');
                    hypos_F = hypo_anypts('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range);
                    fprintf('Fundamental Hypothesis: Done\n\n');
                
                    %Applying ORK            
                    fprintf('Affine ORK:\n');
                    kernels_A = ork_simplesum('affine', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_A, ork_method);
                    fprintf('Affine ORK: Done\n\n');
                
                    fprintf('Homography ORK:\n');
                    kernels_H = ork_simplesum('homography', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_H, ork_method);
                    fprintf('Homography ORK: Done\n\n');
                
                    fprintf('Fundamental ORK:\n');
                    kernels_F = ork_simplesum('fundamental', seq_list, seq_range, data, max_hypos, frame_gap_range, hypos_F, ork_method);
                    fprintf('Fundamental ORK: Done\n\n');

                otherwise
                    error('The chosen ORK philosophy is not valid');
            end
    
        else
    
            fprintf('Loading existing kernels...\n\n');
            target_folder = fullfile('Runs/', sel_dataset, folder_name, sprintf('Run%d', r_i));
    
            kernels_A = read_kernels('affine', target_folder, seq_list);
            kernels_H = read_kernels('homography', target_folder, seq_list);
            kernels_F = read_kernels('fundamental', target_folder, seq_list);
        end
    
        %Motion segmentation
    
        best_alphas = zeros(1, 3);
    
        fprintf('Affine Motion Segmentation:\n');
        best_alphas(1) = moseg_single('Affine', seq_list, seq_range, data, max_hypos, ...
            frame_gap_range, alpha_range, kernels_A, clustering_method);
        fprintf('Affine Motion Segmentation: Done\n\n');
    
        fprintf('Homography Motion Segmentation:\n');
        best_alphas(2) = moseg_single('Homography', seq_list, seq_range, data, max_hypos, ...
            frame_gap_range, alpha_range, kernels_H, clustering_method);
        fprintf('Homography Motion Segmentation: Done\n\n');
    
        fprintf('Fundamental Motion Segmentation:\n');
        best_alphas(3) = moseg_single('Fundamental', seq_list, seq_range, data, max_hypos, ...
            frame_gap_range, alpha_range, kernels_F, clustering_method);
        fprintf('Fundamental Motion Segmentation: Done\n\n');
    
        fprintf('Kernel Addition Motion Segmentation:\n');
        moseg_keradd(seq_list, seq_range, data, max_hypos, frame_gap_range, ...
            alpha_range, best_alphas, kernels_A, kernels_H, kernels_F, clustering_method);
        fprintf('Kernel Addition Motion Segmentation: Done\n\n');
    
        fprintf('Co-Regularization Motion Segmentation:\n');
        moseg_coreg(seq_list, seq_range, data, max_hypos, frame_gap_range, ...
            alpha_range, best_alphas, lambda_range, kernels_A, kernels_H, kernels_F, clustering_method);
        fprintf('Co-Regularization Motion Segmentation: Done\n\n');
    
        fprintf('Subset Constraint Motion Segmentation:\n');
        moseg_subset(seq_list, seq_range, data, max_hypos, frame_gap_range, ...
            alpha_range, best_alphas, gamma_range, kernels_A, kernels_H, kernels_F, clustering_method);
        fprintf('Subset Constraint Motion Segmentation: Done\n\n');
    
        %Save results
        % movefile('Results/Kernels', fullfile('Runs/', sel_dataset, folder_name));
        % movefile('Results/ML', fullfile('Runs/', sel_dataset, folder_name));
        movefile('Results/MoSeg', fullfile('Runs/', sel_dataset, folder_name));
        fprintf('End of Simulation %d with frame gap %d\n\n', r_i, frame_gap);
    end
end  

executionMinutes = toc(tstart) / 60

%% Helper functions
function kernels = read_kernels(method, target_folder, seq_list)

    kernels = cell(1, length(seq_list));
    files = dir(fullfile( target_folder, 'Kernels', method));
    ignore = 0;
    for i = 1:numel(files)
        if startsWith(files(i).name, 'ORK')
            filename = files(i).name;
            temp = load(fullfile(target_folder, 'Kernels', method, filename));
            kernels{i-ignore} = temp.K;
        else
            ignore = ignore + 1;
        end
    end
end

function [seq_range, max_hypos, alpha_range, lambda_range, gamma_range] = load_parameters(sel_dataset)
        load(fullfile('Runs', sel_dataset, 'default_parameters.mat'), 'Parameters')
        seq_range = Parameters.seq_range;
        max_hypos = Parameters.max_hypos;
        alpha_range = Parameters.alpha_range;
        lambda_range = Parameters.lambda_range;
        gamma_range = Parameters.gamma_range;
end

function [sel_dataset, sel_ork, sel_clustering_method, reuse_kernels, run_range, frame_gaps, ork_philosophy] = select_parameters()

    % Define options for each parameter
    dataset_options = {'KT3DMoSeg', 'Hopkins155'};
    ork_options = {'defaultORK', 'altORK1', 'altORK2', 'altORK3'};
    clustering_options = {'kmeans', 'dbscan', 'linkage'};
    reuse_options = {'No', 'Yes'};
    ork_philosophy_options = {'simplesum', 'multiframe', 'multiframe_sum', 'multiframe_score', 'single_frame_gap'};
    
    % Get screen size
    screen_size = get(0, 'ScreenSize');
    screen_width = screen_size(3);
    screen_height = screen_size(4);
    
    % Calculate position of the figure to center it
    figure_width = 600;
    figure_height = 550;  % Increased to accommodate the new dialog box
    figure_position = [(screen_width - figure_width) / 2, (screen_height - figure_height) / 2, figure_width, figure_height];
    
    % Create a figure
    f = figure('Name', 'Select Parameters', 'Position', figure_position, 'NumberTitle', 'off');
    
    % Calculate width for each listbox
    listbox_width = (figure_width - 150) / 4;
    
    % Create listboxes for each parameter
    lb_dataset = uicontrol('Style', 'listbox', 'String', dataset_options, 'Position', [50, 400, listbox_width, 100]);
    lb_ork = uicontrol('Style', 'listbox', 'String', ork_options, 'Position', [50 + listbox_width + 25, 400, listbox_width, 100]);
    lb_clustering = uicontrol('Style', 'listbox', 'String', clustering_options, 'Position', [50 + 2*(listbox_width + 25), 400, listbox_width, 100]);
    lb_reuse = uicontrol('Style', 'listbox', 'String', reuse_options, 'Position', [50 + 3*(listbox_width + 25), 400, listbox_width - 25, 100]);
    
    % Add titles to each listbox
    title_dataset = uicontrol('Style', 'text', 'String', 'Select Dataset:', 'Position', [50, 510, listbox_width, 20], 'HorizontalAlignment', 'center');
    title_ork = uicontrol('Style', 'text', 'String', 'ORK Method:', 'Position', [50 + listbox_width + 25, 510, listbox_width, 20], 'HorizontalAlignment', 'center');
    title_clustering = uicontrol('Style', 'text', 'String', 'Clustering Method:', 'Position', [50 + 2*(listbox_width + 25), 510, listbox_width, 20], 'HorizontalAlignment', 'center');
    title_reuse = uicontrol('Style', 'text', 'String', 'Reuse Kernels?', 'Position', [50 + 3*(listbox_width + 25), 510, listbox_width - 25, 20], 'HorizontalAlignment', 'center');
    
    % Create a listbox for selecting ork_philosophy
    lb_ork_philosophy = uicontrol('Style', 'listbox', 'String', ork_philosophy_options, 'Position', [100, 300, 400, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for ork_philosophy
    title_ork_philosophy = uicontrol('Style', 'text', 'String', 'Select ORK Philosophy:', 'Position', [100, 340, 400, 20], 'HorizontalAlignment', 'left');
    
    % Create an edit box for run range
    edit_run_range = uicontrol('Style', 'edit', 'String', '1:3', 'Position', [100, 200, 400, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for run range
    title_run_range = uicontrol('Style', 'text', 'String', 'Enter a run range (as a single value or an array start:end):', 'Position', [100, 240, 400, 20], 'HorizontalAlignment', 'left');
    
    % Create an edit box for frame gap range
    edit_frame_gap_range = uicontrol('Style', 'edit', 'Position', [100, 100, 400, 40], 'HorizontalAlignment', 'left');
    
    % Add a title for frame gap range
    title_frame_gap_range = uicontrol('Style', 'text', 'String', 'Select the frame gap range (as a single value or an array start:end):', 'Position', [100, 140, 400, 20], 'HorizontalAlignment', 'left');
    
    % Create an OK button
    btn_ok = uicontrol('Style', 'pushbutton', 'String', 'OK', 'Position', [250, 20, 100, 40], ...
                       'Callback', @okCallback);
    
    % Initialize variables to store selections
    sel_dataset = '';
    sel_ork = '';
    sel_clustering_method = '';
    reuse_kernels = '';
    run_range = 1:3;
    frame_gaps = '';
    ork_philosophy = '';
    
    % Callback function for OK button
    function okCallback(~, ~)
        sel_dataset = dataset_options{get(lb_dataset, 'Value')};
        sel_ork = ork_options{get(lb_ork, 'Value')};
        sel_clustering_method = clustering_options{get(lb_clustering, 'Value')};
        reuse_kernels = strcmp(reuse_options{get(lb_reuse, 'Value')}, 'Yes');
        
        % Parse the string entered for run range
        run_range_str = get(edit_run_range, 'String');
        if contains(run_range_str, ':') % Check if a range is entered
            run_range_parts = strsplit(run_range_str, ':');
            run_range = str2double(run_range_parts{1}):str2double(run_range_parts{2});
        else
            run_range = str2double(run_range_str); % Convert single number directly
        end
        
        % Parse the string entered for frame gap range
        frame_gap_range_str = get(edit_frame_gap_range, 'String');
        if contains(frame_gap_range_str, ':') % Check if a range is entered
            frame_gap_range_parts = strsplit(frame_gap_range_str, ':');
            frame_gaps = str2double(frame_gap_range_parts{1}):str2double(frame_gap_range_parts{2});
        else
            frame_gaps = str2double(frame_gap_range_str); % Convert single number directly
        end
        
        % Get the selected ORK philosophy
        ork_philosophy = ork_philosophy_options{get(lb_ork_philosophy, 'Value')};
        
        % Close the figure
        close(f);
    end

    % Wait for the figure to close
    uiwait(f);
end