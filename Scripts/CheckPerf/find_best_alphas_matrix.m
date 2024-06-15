
function [Best_alphas_matrix, Best_alpha_ind_matrix] = find_best_alphas_matrix(Error_matrix, alpha_range, lambda_range, gamma_range, sel_seqs, sel_lambda, sel_gamma)

    methods = {'Affine', 'Homography', 'Fundamental', 'KerAdd', 'CoReg', 'Subset'};
    % % Get the dimensions of Error_matrix
    % [num_runs, num_methods, num_seqs, num_alphas, num_lambdas, num_gammas] = size(Error_matrix);
    
    Best_alphas_matrix = zeros(length(methods), length(sel_seqs));
    Best_alpha_ind_matrix = zeros(length(methods), length(sel_seqs));
    Averaged_matrix = mean(Error_matrix, 1);
    sel_lambda_ind = lambda_range == sel_lambda;
    sel_gamma_ind = gamma_range == sel_gamma;

    for method_ind = 1:length(methods)
        for seq_ind = sel_seqs
            Sequences_methods_errors = squeeze(Averaged_matrix(:, method_ind, seq_ind, :, sel_lambda_ind, sel_gamma_ind));            
            [min_error, min_alpha_ind] = min(Sequences_methods_errors);
            best_alpha = alpha_range(min_alpha_ind);
            Best_alphas_matrix(method_ind, seq_ind) = best_alpha;
            Best_alpha_ind_matrix(method_ind, seq_ind) = min_alpha_ind;
        end
    end
end