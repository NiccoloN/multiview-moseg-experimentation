function [] = moseg_keradd(seq_list, seq_range, data, max_hypos, frame_gap_range, alpha_range, best_alphas, ...
    kernels_A, kernels_H, kernels_F, clustering_method)

    for alpha = alpha_range

        alphas = [alpha, alpha, alpha];
        moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, alphas, kernels_A, kernels_H, kernels_F, clustering_method)
    end

    moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, best_alphas, kernels_A, kernels_H, kernels_F, clustering_method)
end

function moseg(seq_list, seq_range, data, max_hypos, frame_gap_range, alphas, ...
    kernels_A, kernels_H, kernels_F, clustering_method)

    model_type = 'KerAdd';

    %% motion segmentation result save path
    result_path = fullfile('Results/MoSeg/',model_type);
    
    if not(isfolder(result_path))
        mkdir(result_path);
    end
    
    error_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alphas-%d,%d,%d.mat',...
        max_hypos,alphas(1),alphas(2),alphas(3)));
    
    %% motion segmentation on all sequences
    error = [];
    cluster_idx = [];
    
    for s_i = seq_range
        
        seq_name = seq_list{s_i};
        Seq_data = data{s_i};

        ml_path = fullfile('Results/ML',seq_name,model_type);
    
        if not(isfolder(ml_path))
            mkdir(ml_path);
        end

        U_filepath = fullfile(ml_path,sprintf('ML_nhpf-%d_alphas-%d,%d,%d.mat',...
            max_hypos,alphas(1),alphas(2),alphas(3)));
        
        %% Normalize affinity matrix by the Cooccurrence of points
        % Pts_occ = double(Seq_data.visibleSparse);
        % cooc_normalizer = Pts_occ*Pts_occ'+ 0.1;

        K_A = kernels_A{s_i};
        % K_A = K_A./(cooc_normalizer+0.1);
        [K_A,Mask] = func_Adapt_eNN(K_A,alphas(1));
        
        K_H = kernels_H{s_i};
        % K_H = K_H./(cooc_normalizer+0.1);
        [K_H,Mask] = func_Adapt_eNN(K_H,alphas(2));
        
        K_F = kernels_F{s_i};
        % K_F = K_F./(cooc_normalizer+0.1);
        [K_F,Mask] = func_Adapt_eNN(K_F,alphas(3));
        
        %% Sum up all kernels
        AffinityMat = K_A + K_H + K_F + eps;
        % AffinityMat = func_Adapt_eNN(AffinityMat,alphas(1)) + eps;
        
        %% Do spectral clustering
        nMotion = max(Seq_data.GtLabel);
        [cluster_idx{s_i}, eig, U, red_lap_trace] = SpectralClustering_svd(AffinityMat,nMotion,'normalized', clustering_method);
        
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
    end
    
    fprintf('alphas = [%d, %d, %d] mean error = %.2f%%\n',alphas(1),alphas(2),alphas(3),100*mean(error(seq_range)));

    %% Save Results
    save(error_filepath,'error','cluster_idx');
end
