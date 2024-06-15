% Error matrix generation obtained with the best results comparing
% different frame gaps
% -----------------------------------------------------------------
% Parameters
% sel_dataset: dataset used
% run_range: number of runs under study
% seq_range: sequences under study
% max_hypos: numer of hypothesis
% alpha_range: range of alphas used
% lambda_range: range of lambdas used
% gamma_range: range of gammas used
% sel_framegap1: first frame_gap to compare
% sel_framegap2: second frame_gap to compare


function [Error_matrix_best] = comp_framegap(sel_dataset, run_range, seq_range, max_hypos, alpha_range, lambda_range, gamma_range, sel_framegap1, sel_framegap2)


methods = {'Affine', 'Fundamental', 'Homography', 'KerAdd', 'CoReg', 'Subset'};
file_name = sprintf('Error_RandSamp_nhpf-%d_alpha-', max_hypos);
file_extensions = {'', '', '', '', '', ''};

Error_matrix_best = zeros(length(run_range), length(methods), length(seq_range), length(alpha_range), length(lambda_range), ...
    length(gamma_range));


for run_ind = run_range
    folder_name1 = sprintf('Run%d_framegap%d', run_ind, sel_framegap1);
    folder_name2 = sprintf('Run%d_framegap%d', run_ind, sel_framegap2);

    for method_ind = 1:length(methods)
        sel_method = methods{method_ind};
        target_folder1 = fullfile('Runs', sel_dataset, folder_name1, 'MoSeg', sel_method);
        target_folder2 = fullfile('Runs', sel_dataset, folder_name2, 'MoSeg', sel_method);

        for alpha_ind = 1:length(alpha_range)
            
            for lambda_ind = 1:length(lambda_range)
                file_extensions{5} = strcat('_lambda-', num2str(lambda_range(lambda_ind)));

                for gamma_ind = 1:length(gamma_range)
                    file_extensions{6} = strcat('_gamma-', num2str(gamma_range(gamma_ind)));
                    target_file_name1 = fullfile(target_folder1, strcat(file_name, num2str(alpha_range(alpha_ind)), ...
                        file_extensions{method_ind}, '.mat'));
                    target_file_name2 = fullfile(target_folder2, strcat(file_name, num2str(alpha_range(alpha_ind)), ...
                        file_extensions{method_ind}, '.mat'));
                    target_file1 = load(target_file_name1);
                    target_file2 = load(target_file_name2);

                    for seq_ind = 1:length(seq_range)
                        if target_file2.error(seq_ind) <= target_file1.error(seq_ind)
                            Error_matrix_best(run_ind, method_ind, seq_ind, alpha_ind, lambda_ind, gamma_ind) = target_file2.error(seq_ind);
                        else
                            Error_matrix_best(run_ind, method_ind, seq_ind, alpha_ind, lambda_ind, gamma_ind) = target_file1.error(seq_ind);
                        end
                    end
                end
            end
        end
    end
end

end

