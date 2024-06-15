% -------------------------------------------------------------------------
% Error matrix generation
% -------------------------------------------------------------------------
% This function reads files contained in main_folder_name and organizes
% them in a 6-dimensional error matrix.
% 
% The matrix structure follows:
% dim.1: index of run (referring to run_range array)
% dim.2: used method (Ordered as: Affine, Homography, Fundamental, KerAdd,
%                       CoReg, Subset)
% dim.3: sequence number
% dim.4: index of alpha (referring to alpha_range array)
% dim.5: index of lambda (referring to lambda_range array)
% dim.6: index of gamma (referring to gamma_range array)
% 
% The required parameters are:
% main_folder_name: directory where folders called 'RunX' are located,
%                   where X is the run number
% seq_range: array of considered sequences
% max_hypos: number of maximum hypotesis per frame
% alpha_range: array of values of alpha
% lambda_range: array of values of lambda
% gamma_range: array of values of gamma

function [Error_matrix, Clusters_matrix] = error_analisys(main_folder_name, run_range, seq_range, max_hypos, alpha_range, lambda_range, gamma_range)


methods = {'Affine', 'Homography', 'Fundamental', 'KerAdd', 'CoReg', 'Subset'};
file_name = sprintf('Error_RandSamp_nhpf-%d_alpha', max_hypos);

Error_matrix = zeros(length(run_range), length(methods), length(seq_range), length(alpha_range), length(lambda_range), ...
    length(gamma_range));

Clusters_matrix = cell(length(run_range), length(methods), length(seq_range), length(alpha_range), length(lambda_range), ...
    length(gamma_range));

num_seqs = length(seq_range);
num_alphas = length(alpha_range);
num_lambdas = length(lambda_range);
num_gammas = length(gamma_range);

for ind = 1:length(run_range)
    run_ind = run_range(ind);
    folder_name = fullfile(main_folder_name, sprintf('/Run%d', run_ind));

    for method_ind = 1:length(methods)
        sel_method = methods{method_ind};
        target_folder = fullfile(folder_name, 'MoSeg', sel_method);

        parfor alpha_ind = 1:num_alphas
            for lambda_ind = 1:num_lambdas
                for gamma_ind = 1:num_gammas

                    % old runs
                    % file_extensions = {
                    %         sprintf('-%d', alpha_range(alpha_ind))
                    %         sprintf('-%d', alpha_range(alpha_ind)) 
                    %         sprintf('-%d', alpha_range(alpha_ind))
                    %         sprintf('-%d', alpha_range(alpha_ind))
                    %         sprintf('-%d_lambda-%g', alpha_range(alpha_ind), lambda_range(lambda_ind))
                    %         sprintf('-%d_gamma-%g', alpha_range(alpha_ind), gamma_range(gamma_ind))
                    %     };
                    % new runs
                    file_extensions = {
                            sprintf('-%d', alpha_range(alpha_ind))
                            sprintf('-%d', alpha_range(alpha_ind)) 
                            sprintf('-%d', alpha_range(alpha_ind))
                            sprintf('s-%d,%d,%d', alpha_range(alpha_ind), alpha_range(alpha_ind), alpha_range(alpha_ind))
                            sprintf('s-%d,%d,%d_lambda-%g', alpha_range(alpha_ind), alpha_range(alpha_ind), alpha_range(alpha_ind), lambda_range(lambda_ind))
                            sprintf('s-%d,%d,%d_gamma-%g', alpha_range(alpha_ind), alpha_range(alpha_ind), alpha_range(alpha_ind), gamma_range(gamma_ind))
                        };
                    target_file_name = fullfile(target_folder, ...
                        strcat(file_name, file_extensions{method_ind}, '.mat'));
                    target_file = load(target_file_name);

                    for seq_ind = 1:num_seqs
                        Error_matrix(ind, method_ind, seq_ind, alpha_ind, lambda_ind, gamma_ind) = target_file.error(seq_range(seq_ind));
                        Clusters_matrix(ind, method_ind, seq_ind, alpha_ind, lambda_ind, gamma_ind) = target_file.cluster_idx(seq_range(seq_ind));
                    end
                end
            end
        end
    end
end

end