function [] = moseg_coreg(seq_list, seq_range, data, max_hypos, frame_gap_range, alpha_range, best_alphas, lambda_range, ...
    kernels_A, kernels_H, kernels_F, clustering_method)
    
    for alpha = alpha_range

        alphas = [alpha, alpha, alpha];
        moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, alphas, lambda_range, kernels_A, kernels_H, kernels_F, clustering_method)
    end

    moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, best_alphas, lambda_range, kernels_A, kernels_H, kernels_F, clustering_method)
end

function moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, alphas, lambda_range, ...
    kernels_A, kernels_H, kernels_F, clustering_method)

    model_type = 'CoReg';
    
    for lambda = lambda_range

        %% motion segmentation result save path
        result_path = fullfile('Results/MoSeg/',model_type);
        
        if not(isfolder(result_path))
            mkdir(result_path);
        end
        
        error_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alphas-%d,%d,%d_lambda-%g.mat',...
            max_hypos,alphas(1),alphas(2),alphas(3),lambda));
        
        %% motion segmentation on all sequences
        error = [];
        cluster_idx = cell(1, length(seq_range));
        
        parfor s_i = seq_range
            
            seq_name = seq_list{s_i};
            Seq_data = data{s_i};

            ml_path = fullfile('Results/ML',seq_name,model_type);
    
            if not(isfolder(ml_path))
                mkdir(ml_path);
            end

            U_filepath = fullfile(ml_path,sprintf('ML_nhpf-%d_alphas-%d,%d,%d.mat',...
                max_hypos,alphas(1),alphas(2),alphas(3)));
            
            %% Normalize affinity matrix by the Cooccurrence of points            
            % Pts_occ = double(Seq_data.visibleSparse);    % point occurrence across all frames
            % cooc_normalizer = Pts_occ*Pts_occ'+ 0.1;     % points cooccurrence

            K_A = kernels_A{s_i};
            % K_A = K_A./(cooc_normalizer+0.1);
            [K_A,~] = func_Adapt_eNN(K_A,alphas(1));
            
            K_H = kernels_H{s_i};
            % K_H = K_H./(cooc_normalizer+0.1);
            [K_H,~] = func_Adapt_eNN(K_H,alphas(2));
            
            K_F = kernels_F{s_i};
            % K_F = K_F./(cooc_normalizer+0.1);
            [K_F,~] = func_Adapt_eNN(K_F,alphas(3));
            
            %% CoRegularization motion segmentation
                                  
            K = [];
            K(:,:,1) = K_A+eps;
            K(:,:,2) = K_H+eps;
            K(:,:,3) = K_F+eps;
            
            nMotion = max(Seq_data.GtLabel);
            epsilon = 1e-6;
            MaxItr = 15;
            
            %% CoRegularization
            [U_CoReg, ~ , ~, ~, red_lap_trace] = func_CoRegularize_eig(K,nMotion,lambda,epsilon,MaxItr);
            
            %% Normalize
            U_All = [];
            for k_i = 1:size(U_CoReg,3)
                U_All = [U_All func_L2Normalize(U_CoReg(:,:,k_i),2)];
            end
            
            %%% Normalize
            U = func_L2Normalize(U_All,2);
            
            switch clustering_method
                case 'kmeans'
                    cluster_idx{s_i} = kmeans(U, nMotion, 'replicates',500, 'start', 'cluster', ...
                        'EmptyAction', 'singleton');
                case 'dbscan'
                    for t = 0.01:0.01:0.5
                        Grps = dbscan(U,t,1);
                        labels = unique(Grps);
                        labels = labels(labels ~= -1);
                        if length(labels) == nMotion
                            % fprintf('best threshold: %g\n', t);
                            break;
                        end
                    end
                    cluster_idx{s_i} = Grps;
                case 'linkage'
                    cluster_idx{s_i} = cluster(linkage(U), 'maxclust', nMotion);
                otherwise
                    fprintf('Error: The chosen clustering method is not valid (ma matlab puzza e non posso mettere un errore in un parfor)');
            end

            if isrow(cluster_idx{s_i})
                cluster_idx{s_i} = cluster_idx{s_i}';
            end
            
            %%% Evaluate Classification Error
            [error(s_i), true_cluster_labels] = Misclassification(cluster_idx{s_i},Seq_data.GtLabel);
        
            cluster_idx_copy1 = cluster_idx{s_i};
            cluster_idx_copy2 = cluster_idx_copy1;

            for clust_ind = 1:length(true_cluster_labels)
                indices = (cluster_idx_copy2 == true_cluster_labels(clust_ind));
                cluster_idx_copy1(indices) = clust_ind;
            end

            cluster_idx{s_i} = cluster_idx_copy1;
                        
            % fprintf('Sequence %s Error = %.2f%% \n',seq_name,100*error(s_i));

            clusters = cluster_idx{s_i};
            parsave1(U_filepath, clusters, red_lap_trace)
        end
        
        fprintf('alphas = [%d, %d, %d], lambda = %.3f, mean error = %.2f%%\n',alphas,lambda,100*mean(error(seq_range)));
        
        %% Save Results
        save(error_filepath,'cluster_idx','error');
    end
end

function parsave1(U_filepath, clusters, red_lap_trace)
    save(U_filepath,'clusters', 'red_lap_trace');
end
