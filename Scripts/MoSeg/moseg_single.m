%% script to conduct motion segmentation

function [best_alpha] = moseg_single(method, seq_list, seq_range, data, max_hypos, frame_gap_range, alpha_range, kernels, clustering_method)
    
    best_alpha = -1;
    best_mean_error = inf;
    
    %% Evaluate all power scaling parameters
    for alpha = alpha_range
        
        %%% motion segmentation result save path
        result_path = fullfile('Results/MoSeg/',method);
        
        if not(isfolder(result_path))
            mkdir(result_path);
        end
        
        error_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
            max_hypos,alpha));
        
        error = [];
        cluster_idx = [];
        
        %% Evaluate all sequences
        for s_i = seq_range
            
            seq_name = seq_list{s_i};
            Seq_data = data{s_i};
            K = kernels{s_i};

            ml_path = fullfile('Results/ML',seq_name,method);
        
            if not(isfolder(ml_path))
                mkdir(ml_path);
            end

            U_filepath = fullfile(ml_path,sprintf('ML_nhpf-%d_alpha-%g.mat',...
                max_hypos,alpha));
            
            %% Normalize affinity matrix by the Cooccurrence of points       
            % Pts_occ = double(Seq_data.visibleSparse);    % point occurrence across all frames
            % cooc_normalizer = Pts_occ*Pts_occ'+ 0.1;     % points cooccurrence
            % 
            % K = K./cooc_normalizer;
            
            %% Sparsify affinity matrix
            K = func_Adapt_eNN(K,alpha)+eps;
            
            %% Do spectral clustering
            nMotion = max(Seq_data.GtLabel);
            [cluster_idx{s_i}, eig, LapKer, red_lap_trace] = SpectralClustering_svd(K,nMotion,'normalized', clustering_method);
            
            if isrow(cluster_idx{s_i})
                cluster_idx{s_i} = cluster_idx{s_i}';
            end
            
            %% Eval Classification Error Rate
            [error(s_i), true_cluster_labels] = Misclassification(cluster_idx{s_i},Seq_data.GtLabel);
            
            cluster_idx_copy1 = cluster_idx{s_i};
            cluster_idx_copy2 = cluster_idx_copy1;

            for clust_ind = 1:length(true_cluster_labels)
                indices = (cluster_idx_copy2 == true_cluster_labels(clust_ind));
                cluster_idx_copy1(indices) = clust_ind;
            end

            cluster_idx{s_i} = cluster_idx_copy1;

            % fprintf('seq-%d error=%.2f%%\n',s_i,100*error(s_i));

            clusters = cluster_idx{s_i};
            save(U_filepath,'clusters', 'red_lap_trace');

            % figure
            % hold on
            % grid on
            % label_colors = ["red", "blue", "green"];
            % for i = 1:size(clusters,1)
            %     p = LapKer(i, :);
            %     marker = '.';
            %     if clusters(i) ~= Seq_data.GtLabel(i)
            %         marker = 'x';
            %     end
            %     if length(p) == 2
            %         plot(p(1), p(2), marker, Color=label_colors(clusters(i)));
            %     else
            %         plot3(p(1), p(2), p(3), marker, Color=label_colors(clusters(i)));
            %     end
            % end
        end
        
        mean_error = mean(error(seq_range));
        if (mean_error < best_mean_error)
            best_mean_error = mean_error;
            best_alpha = alpha;
        end
        
        fprintf('alpha = %d mean error = %.2f%%\n',alpha,100*mean_error);
        
        %% Save Results
        
        save(error_filepath,'error','cluster_idx');
    end
end
