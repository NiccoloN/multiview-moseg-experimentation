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

function[Kernels] = ork_multiframe_score(method, seq_list, seq_range, data, max_hypos, frame_gap_range, hypos, ork_method)
    
    model_type = lower(method);  
    [fitfn, resfn, degenfn, psize, numpar] = getModelParam(model_type);

    Kernels = cell(1, length(seq_list));
    normalized_bad_hypos = zeros(length(seq_range), 1);

    h = round(0.1 * max_hypos);
    
    for s_i = seq_range
        
        seq_name = seq_list{s_i};
        seq_data = data{s_i};
        seq_hypos = hypos{s_i};

        cock_normalizer = zeros(seq_data.nSparsePoints);

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

        for f_i = 1:seq_data.nFrames-min(frame_gap_range) % use min to consider also the last frames
            final_frame = min(f_i + max(frame_gap_range), seq_data.nFrames);
            max_frame_gap = final_frame - f_i;
            reduced_frame_gap_range = frame_gap_range(frame_gap_range <= max_frame_gap);

            Res_framegap = zeros(seq_data.nSparsePoints, max_hypos, max_frame_gap);
            hypo_inds = -ones(seq_data.nSparsePoints, max_hypos, max_frame_gap);
            visible_pts_ind = zeros(seq_data.nSparsePoints, length(reduced_frame_gap_range),'logical');

            for frame_gap = reduced_frame_gap_range
                r = f_i;
                v = f_i+frame_gap;
                mdl_idx = find(seq_hypos.r == r & seq_hypos.v == v)';
                
                %% Select points visible on both frames
                visible_pts_ind(:, frame_gap) = seq_data.visibleSparse(:,r) & seq_data.visibleSparse(:,v);
                
                y1 = seq_data.ySparse(:,visible_pts_ind(:, frame_gap),r);
                y2 = seq_data.ySparse(:,visible_pts_ind(:, frame_gap),v);
                
                %% Normalise raw correspondences
                dat_img_1 = normalise2dpts(y1);
                dat_img_2 = normalise2dpts(y2);
                normalized_data = [dat_img_1 ; dat_img_2];

                for h_i = mdl_idx
                    %% Calculate residual
                    Res_allpts = zeros(seq_data.nSparsePoints, 1);
                    Res_allpts(visible_pts_ind(:, frame_gap)) = feval(resfn,seq_hypos.H(:,h_i),normalized_data);
                    
                    Res_framegap(:, mod(h_i-1, max_hypos)+1, frame_gap) = Res_allpts;
                end

                [Res_sorted, resinx] = sort(Res_framegap(:,:,frame_gap), 2);
                hypo_inds(visible_pts_ind(:, frame_gap), :, frame_gap) = resinx(visible_pts_ind(:, frame_gap),:);
            end

            Scores = zeros(seq_data.nSparsePoints, max_hypos, max_frame_gap);
            for framegap_ind = reduced_frame_gap_range
                best_hypo_inds = hypo_inds(:,1:h,framegap_ind);
                point_inds = 1:seq_data.nSparsePoints;
                for point_ind = point_inds(visible_pts_ind(:, frame_gap))
                    for other_framegap_ind = reduced_frame_gap_range(reduced_frame_gap_range~=framegap_ind)
                        [~,~,other_scores] = intersect(best_hypo_inds(point_ind,:), hypo_inds(point_ind,:,other_framegap_ind), 'stable');
                        Scores(point_ind, best_hypo_inds(point_ind,:), other_framegap_ind) = max_hypos+1-other_scores;
                    end
                    own_scores = 1:h;
                    Scores(point_ind, best_hypo_inds(point_ind,:), framegap_ind) = max_hypos+1-own_scores;
                end
            end
            
            Scores = sum(Scores, 3);
            visible_pts_ind = all(visible_pts_ind, 2);

            %% Compute ORK kernel
            [~, resinx] = sort(Scores,2,'descend');
            resinx = resinx(visible_pts_ind,:);
            
            % h = round(0.1 * size(Res,2));
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
            K_ORK = zeros(nnz(visible_pts_ind));
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

            K_temp(visible_pts_ind,visible_pts_ind) = K_ORK;            
            K = K + K_temp;

            cock_normalizer(visible_pts_ind,visible_pts_ind) = cock_normalizer(visible_pts_ind,visible_pts_ind) + 1;
        end

        % K = K ./ (cock_normalizer + eps);

        Pts_occ = double(seq_data.visibleSparse);    % point occurrence across all frames
        cooc_normalizer = Pts_occ*Pts_occ'+ 0.1;     % points cooccurrence

        K = K./cooc_normalizer;

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

