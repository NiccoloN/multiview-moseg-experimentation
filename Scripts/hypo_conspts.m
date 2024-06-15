% -------------------------------------------------------------------------
%                           Hypo consecutive pts
% -------------------------------------------------------------------------
% This function generates hypotheses using points selected from frames
% distanced by values between 1 and frame_gap_range.
% It is MANDATORY for an hypothesis to be generated that the points
% remain visible for the frames in between.

%% Sample affine hypotheses

function [Hypos] = hypo_conspts(method, seq_list, seq_range, Data, max_hypos, frame_gap_range)
    
    %% Sample hypotheses from all sequences
    model_type = lower(method);

    Hypos = cell(1, length(seq_list));

    parfor s_i = seq_range
        
        Seq_data = Data{s_i};
        num_frames = Seq_data.nFrames;
        
        %%% Initialize hypotheses
        Seq_hypos = struct;
        Seq_hypos.H = [];
        Seq_hypos.r = [];
        Seq_hypos.v = [];
        Seq_hypos.supp = [];
        Seq_hypos.cond_numb = [];
       
        for f_i = 1:num_frames-min(frame_gap_range)
        
            %% Prepare candidate data
            r = f_i;
            final_frame = min(r+max(frame_gap_range), num_frames);
            max_frame_gap = final_frame - r;
            
            %%% Select points visible on both frames
            visible_pts_ind = all(Seq_data.visibleSparse(:,r:final_frame), 2);

            for frame_gap = frame_gap_range(frame_gap_range <= max_frame_gap)
                v = r+frame_gap;

                y1 = Seq_data.ySparse(:,visible_pts_ind,r);
                y2 = Seq_data.ySparse(:,visible_pts_ind,v);
        
                %% Normalise raw correspondences.
                dat_img_1 = normalise2dpts(y1);
                dat_img_2 = normalise2dpts(y2);
                normalized_data = [dat_img_1; dat_img_2];
                
                % Maximum CPU seconds allowed
                lim = 20;
            
                % Storage.
                par = cell(2,1);
                res = cell(2,1);
                inx = cell(2,1);
                tim = cell(2,1);
                hit = cell(2,4);
                met = char('Random','Multi-GS');

                if frame_gap == min(frame_gap_range)
                    % Random sampling.
                    [par{1}, res{1}, inx{1}, tim{1}] = randomSampling(lim,normalized_data,max_hypos,model_type);
                    % Guided-sampling using the Multi-GS method. (alternative sampling strategy)
                    % [par{1}, res{1}, inx{1}, tim{1}] = multigsSampling(lim,normalized_data,max_hypos,10,model_type);
                    supp = inx{1};
                else
                    % Reuse support points
                    [par{1}, res{1}] = fit_hypos(normalized_data,supp,max_hypos,model_type);
                    inx{1} = supp;
                end
                
                cond_numbers = zeros(1, size(par{1}, 2));
                for hypo_ind = 1:size(par{1}, 2)
                    matrix_coeffs = par{1}(:, hypo_ind);
                    matrix = [matrix_coeffs(1:3), matrix_coeffs(4:6), matrix_coeffs(7:9)];
                    cond_numbers(hypo_ind) = cond(matrix); % Horizontal array of condition numbers for each hypothesis matrix
                end
    
                % Guided-sampling using the Multi-GS method. (alternative sampling strategy)
                % [ par{2} res{2} inx{2} tim{2} ] = multigsSampling(lim,normalized_data,max_NumHypoPerFrame,10,model_type);
                
                %% Accumulate hypotheses
                Seq_hypos.H = [Seq_hypos.H  par{1}];                            %I parametri del modello
                Seq_hypos.r = [Seq_hypos.r ; r*ones(size(par{1},2),1)];         %Un vettore pieno di numeri r (vedi su per cosa è r)
                Seq_hypos.v = [Seq_hypos.v ; v*ones(size(par{1},2),1)];         %Un vettore pieno di numeri v (vedi su per cosa è v)
                Seq_hypos.supp = [Seq_hypos.supp  inx{1}];                      %Il vettore degli indici dei punti che usa per fittare l'ipotesi
                Seq_hypos.cond_numb = [Seq_hypos.cond_numb; cond_numbers];      %Una matrice (n°frame-frame_gap)*(n°hypo) di condition numbers
            end
        end
        
        Hypos{s_i} = Seq_hypos;
        fprintf('Finish %d-th seq\n',s_i);
    end
end
