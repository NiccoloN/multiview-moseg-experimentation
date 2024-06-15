%% Preliminary operations and values to be set
clear all
clc

% Parameters
sel_seq = 4;
frame_range = [1 10 20];
method = 'KerAdd';
label_colors = ["red", "blue", "green", "yellow", "magenta"];
sel_run = 1;
sel_alpha = 8;
[sel_dataset, ork_method, clustering_method, run_range, frame_gap_range, ork_philosophy, sel_seqs] = select_parameters();

%% Set Path
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath('..'));

%% Load data
dataset_path = '../Datasets/KT3DMoSeg';

temp = load(fullfile(dataset_path, 'SeqList.mat'));
seq_list = temp.SeqList;

seq_name = seq_list{sel_seq};
seq_path = fullfile(dataset_path, 'OriginalSequences/', seq_name);

temp = load(fullfile(dataset_path, seq_name));
data = temp.Data;

% Detect if there are many simulations differing only for the selected sequences and act accordingly
subfolders = dir(fullfile('../Runs', sel_dataset));
subfolders = subfolders([subfolders.isdir]);
subfolder_names = {subfolders.name};
similar_subfolders = subfolder_names(startsWith(subfolder_names, sprintf('ORK%s-%s_CLUSTERING%s_FRAMEGAP%d_SEQS', ork_philosophy, ork_method, clustering_method, frame_gap_range)));
switch numel(similar_subfolders) 
    case 0
        error('There is no simulation performed with the selected parameters. Please run one and retry')
    case 1
        folder_name = fullfile('Runs', sel_dataset, similar_subfolders{1});
    otherwise
        error('There is more than one simulation with the same parameters, but evaluated on a different sequence range. Please remove the unwanted one and retry')
end

% temp = load(fullfile(folder_name, sprintf('Run%d', sel_run), 'ML', seq_name, method, sprintf('ML_nhpf-500_alpha-%d.mat', sel_alpha)));
temp = load(fullfile(folder_name, sprintf('Run%d', sel_run), 'ML', seq_name, method, sprintf('ML_nhpf-500_alphas-%d,%d,%d.mat', sel_alpha, sel_alpha, sel_alpha)));
clusters = temp.clusters;

%% View points

figure
rows = length(frame_range);

for f_i = 1:rows

    im = imread(fullfile(seq_path, sprintf('%06d.png', frame_range(f_i))));

    subplot(rows, 2, f_i * 2 - 1)
    imshow(im);
    hold on

    points = data.ySparse(:,:,frame_range(f_i));
    points = points.';
    label = data.GtLabel;
    
    for i = 1:size(points, 1)
        plot(points(i,1), points(i,2), 'o', Color=label_colors(label(i)));
    end

    subplot(rows, 2, f_i * 2)
    imshow(im);
    hold on
    
    for i = 1:size(points, 1)
        plot(points(i,1), points(i,2), 'o', Color=label_colors(clusters(i)));
    end
end

sgtitle(sprintf('Dataset: %s\nORK %s evaluated using %s method\nClustering: %s\nFrame gap: up to %d frames\nSel alpha: %d', sel_dataset, ork_philosophy, ork_method, clustering_method, frame_gap_range, sel_alpha))

%% Helper functions
function [sel_dataset, ork_method, clustering_method, run_range, frame_gap_range, ork_philosophy, sel_seqs] = select_parameters()

    % Define options for each parameter
    dataset_options = {'KT3DMoSeg', 'Hopkins155'};
    ork_options = {'defaultORK', 'altORK1', 'altORK2', 'altORK3'};
    clustering_options = {'kmeans', 'dbscan', 'linkage'};
    ork_philosophy_options = {'simplesum', 'multiframe','single_frame_gap'};
    
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