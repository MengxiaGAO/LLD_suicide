function [pval_pos, pval_neg, pval_both, prediction_r] = CPM_LOOCV_permutation_p_value_diagonal_regress (no_iterations, thresh_pos, thresh_neg, all_mats, all_behav, all_cov, corr_type)
%% This function is to calculate the permutation p value of n-fold CPM
% Input:
%       no_iterations:  the number of permutations (N, a number)
%       thresh_pos:     the optimal p threshold of positive network (e.g., thresh_pos(1,1))
%       thresh_neg:     the optimal p threshold of negative network (e.g., thresh_neg(1,1))
%       all_mats:       an X*Y*Z 3D matrix
%       all_behav:      a Z*1 vector
%       all_cov:        a Z*M vector contains different covariables
%       corr_type:      correlation type: 'Pearson', 'Kendall', 'Spearman'

% Output:
%       pval_pos:     p_value of permutation (positive)
%       pval_neg:     p value of permutation (negative)
%       pval_both:    p value of permutation (combined features)
%       prediction_r: a list of all r values in the permutation

%%%%% The script is adopted from Xinlin Shen (2017, Nature Protocols). 
%%%%% Written by Mengxia Gao, PhD student of Dept. of Psychology, the University of Hong Kong.
%%%%% mengxia.gao@gmail.com, 20200310
%%%%% updated on 20201022

no_sub = size(all_mats,3);
prediction_r = zeros(no_iterations+1, 3);

fprintf('\n Performing true prediction');
[true_r_pos, true_r_neg, true_r_both, ~, ~, ~] = CPM_LOOCV_diagonal_regress(thresh_pos, thresh_neg, all_mats, all_behav, all_cov, corr_type);
prediction_r(1,1) = true_r_pos;
prediction_r(1,2) = true_r_neg;
prediction_r(1,3) = true_r_both;

    
    % begin permutation
    for it = 1:no_iterations
        fprintf('\n Performing interation %d out of %d', it, no_iterations);
        new_behav = all_behav (randperm(no_sub));    % random the behavior
        [prediction_r(it+1,1), prediction_r(it+1,2), prediction_r(it+1,3), ~, ~, ~] = CPM_LOOCV_diagonal_regress(thresh_pos, thresh_neg, all_mats, new_behav, all_cov, corr_type);
    end

    % calculate the true value
    true_prediction_r_pos = prediction_r(1,1);
    true_prediction_r_neg = prediction_r(1,2);
    true_prediction_r_both = prediction_r(1,3);

    sorted_prediction_r_pos = sort (prediction_r(:,1),'descend');
    position_pos            = find (sorted_prediction_r_pos == true_prediction_r_pos);
    pval_pos                = position_pos(1)/(no_iterations+1)

    sorted_prediction_r_neg = sort (prediction_r(:,2),'descend');
    position_neg            = find (sorted_prediction_r_neg == true_prediction_r_neg);
    pval_neg                = position_neg(1)/(no_iterations+1)
    
    sorted_prediction_r_both = sort (prediction_r(:,3),'descend');
    position_both            = find (sorted_prediction_r_both == true_prediction_r_both);
    pval_both                = position_both(1)/(no_iterations+1)

end