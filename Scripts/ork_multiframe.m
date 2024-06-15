% -------------------------------------------------------------------------
%                             ORK multiframe
% -------------------------------------------------------------------------
% This function generates a subaffinity matrix for each initial frame,
% built considering hypotheses performances along the next frame_gap_range
% number of frames.
% Then, every subaffinity matrix is summed to obtain the final affiniy
% matrix K
%
% WARNING: this function works ONLY with hypotheses generated with the
% hypo_conspts method.

function[Kernels] = ork(method, seq_list, seq_range, data, max_hypos, frame_gap_range, hypos, ork_method)
    
    model_type = lower(method);  
    [fitfn, resfn, degenfn, psize, numpar] = getModelParam(model_type);

    Kernels = cell(1, length(seq_list));
    normalized_bad_hypos = zeros(length(seq_range), 1);
    
    for s_i = seq_range
        
        seq_name = seq_list{s_i};
        seq_data = data{s_i};
        seq_hypos = hypos{s_i};

        hypo_occurencies = zeros(size(seq_hypos.cond_numb));
        
        %% Kernel result save path
        save_path = fullfile('Results/Kernels/',model_type);
        
        if not(isfolder(save_path))
            mkdir(save_path);
        end
        
        kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
            seq_name,max_hypos));
        
        %% Compute kernel by accumulating all frame pairs
        K = zeros(seq_data.nSparsePoints);
        max_frame_gap = max(frame_gap_range);

        for f_i = 1:seq_data.nFrames-max_frame_gap
            Res_framegap = zeros(seq_data.nSparsePoints, max_hypos, max_frame_gap);

            for frame_gap_ind = frame_gap_range
            mdl_idx = find(seq_hypos.r == f_i & seq_hypos.v == f_i+frame_gap_ind)'; 

                for h_i = mdl_idx
                    r = seq_hypos.r(h_i);
                    v = seq_hypos.v(h_i);
                    
                    %% Select points visible on both frames
                    visible_pts_ind = seq_data.visibleSparse(:,r) & seq_data.visibleSparse(:,v);
                    
                    y1 = seq_data.ySparse(:,visible_pts_ind,r);
                    y2 = seq_data.ySparse(:,visible_pts_ind,v);
                    
                    %% Normalise raw correspondences
                    dat_img_1 = normalise2dpts(y1);
                    dat_img_2 = normalise2dpts(y2);
                    normalized_data = [dat_img_1 ; dat_img_2];
                    
                    %% Calculate residual
                    Res_allpts = zeros(seq_data.nSparsePoints, 1);
                    Res_allpts(visible_pts_ind) = feval(resfn,seq_hypos.H(:,h_i),normalized_data);
                    Res_framegap(:, mod(h_i-1, max_hypos)+1, frame_gap_ind) = Res_allpts;           
                end
            end
            
            Res = sum(Res_framegap,3);
            visibility_counter = zeros(seq_data.nSparsePoints, 1);
            for framegap_ind = frame_gap_range
                v = f_i+framegap_ind;
                visible_pts_index = seq_data.visibleSparse(:,f_i) & seq_data.visibleSparse(:,v);
                visibility_counter = visibility_counter + visible_pts_index;
            end
            valid_pts_ind = find(visibility_counter~=0);
            Res = Res(valid_pts_ind,:)./visibility_counter(valid_pts_ind);
            
            %% Compute ORK kernel
            [Res_sorted, resinx] = sort(Res,2);
            
            h = round(0.1 * size(Res,2));
            selected_resinx = resinx(:, 1:h);

            for h_i = 1:size(seq_hypos.cond_numb, 2)
                hypo_occurencies(f_i, h_i) = sum(selected_resinx == h_i, "all");
            end

            max_res_diffs = ones(1, size(Res, 2));
            for hypo_i = 1:size(Res, 2)
                [inliers, ~] = find(selected_resinx == hypo_i);
                hypo_residuals = unique(Res(inliers, hypo_i));
                if (length(hypo_residuals) > 1)
                    max_res_diffs(hypo_i) = max(hypo_residuals) - min(hypo_residuals);
                end
            end
            
            K_temp = zeros(seq_data.nSparsePoints);

            i = 1; % questo ci vuole perché MATLAB è stupido in culo
            K_ORK = zeros(size(Res, 1));
            switch ork_method
                case 'defaultORK'
                    K_ORK = computeIntersection(resinx', resinx', h);
                case 'altORK1'
                    parfor i = 1:size(K_ORK, 1)
                        for j = 1:i
                            intersection = intersect(selected_resinx(i, :), selected_resinx(j, :));
                            res_similarities = ones(1, length(intersection)) - abs(Res(i, intersection) - Res(j, intersection)) ./ max_res_diffs(intersection);
                            K_ORK(i, j) = sum(res_similarities) / h;
                        end
                    end
                case 'altORK2'
                    h_altORK2 = 2*h;
                    parfor i = 1:size(K_ORK, 1)
                        for j = 1:i
                            selected_resinx_i = selected_resinx(i, :);
                            selected_resinx_j = selected_resinx(j, :);
                            intersection = intersect(selected_resinx_i, selected_resinx_j);
                            [~, importances_i] = sort(selected_resinx_i(ismember(selected_resinx_i, intersection)));
                            [~, importances_j] = sort(selected_resinx_j(ismember(selected_resinx_j, intersection)));
                            hypo_importance = (importances_i + importances_j) / (h_altORK2 * 2);
                            hypo_importance = importance_function(hypo_importance);
                            res_similarities = ones(1, length(intersection)) - abs(Res(i, intersection) - Res(j, intersection)) ./ max_res_diffs(intersection);
                            K_ORK(i, j) = sum(res_similarities .* hypo_importance);    
                        end
                    end
                case 'altORK3'
                    h_altORK3 = 2*h;
                    parfor i = 1:size(K_ORK, 1)
                        for j = 1:i
                            selected_resinx_i = selected_resinx(i, :);
                            selected_resinx_j = selected_resinx(j, :);
                            intersection = intersect(selected_resinx_i, selected_resinx_j);
                            [~, importances_i] = sort(selected_resinx_i(ismember(selected_resinx_i, intersection)));
                            [~, importances_j] = sort(selected_resinx_j(ismember(selected_resinx_j, intersection)));
                            hypo_importance = (importances_i + importances_j) / (h_altORK3 * 2);
                            hypo_importance = importance_function(hypo_importance);
                            res_similarities = ones(1, length(intersection));
                            K_ORK(i, j) = sum(res_similarities .* hypo_importance);    
                        end
                    end
                otherwise
                    error('The chosen ORK method is not valid');
            end                    
            for i = 1:size(K_ORK, 1)
                for j = 1:i
                    K_ORK(j, i) = K_ORK(i, j);
                    assert(~isnan(K_ORK(i, j)));
                end
            end

            K_temp(valid_pts_ind,valid_pts_ind) = K_ORK;            
            K = K + K_temp;
        end

        bad_cond_numb = seq_hypos.cond_numb;
        bad_cond_numb(bad_cond_numb < median(bad_cond_numb, "all")) = 0;
        bad_cond_numb(bad_cond_numb ~= 0) = 1;

        bad_hypo_occurencies = bad_cond_numb.*hypo_occurencies;
        total_bad_hypo = sum(bad_hypo_occurencies, "all");
        normalized_bad_hypos(s_i) = total_bad_hypo / sum(hypo_occurencies, "all");
        
        %% Save Results
        Kernels{s_i} = K;
        save(kernel_filepath,'K');
        
        fprintf('Finish %d-th seq\n',s_i);
    end

    save(fullfile(save_path, 'normalized_bad_hypos'), 'normalized_bad_hypos');
end

function y = importance_function(x)
    y = 1 - ones(size(x)) ./ (ones(size(x)) + exp(-10 * (x - 0.5)));
    % y = 1 - x;
end
