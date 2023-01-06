function [thresh_results,thresh_pos,thresh_neg] = CPM_select_optimal_p_diagonal_regress(thresh_range, all_mats, all_behav, all_cov, corr_type)
%% This function is to select the optimal p threshold values by trying different thresholds
% Input:
%       thresh_range:   a vector contains a p value range (e.g. [0.001:0.001:0.05], [0.0001:0.0001:0.05])
%       all_mats:       an X*Y*Z 3D matrix
%       all_behav:      a Z*1 vector
%       all_cov:        a Z*M vector contains different covariables
%       corr_type:      correlation type: 'Pearson', 'Kendall', 'Spearman'

% Output:
%       thresh_results: r and p values in all testing thresholds
%       thresh_pos:     the optimal r and p values of the postive network
%       thresh_neg:     the optimal r and p values of the negative network

%%%%% The script is adopted from Xinlin Shen (2017, Nature Protocols). 
%%%%% Written by Mengxia Gao, PhD student of Dept. of Psychology, the University of Hong Kong.
%%%%% mengxia.gao@gmail.com, 20200309
%%%%% updated on 20201022

value = thresh_range;
thresh_results = zeros(length(thresh_range),5);

for thresh = thresh_range
tol = 1.e-6;
position=find(abs(value-thresh)<tol);

fprintf('\n Calculating p threshold # %6.3f. Try %d out of %d',thresh, position,length(thresh_range));

[R_pos,R_neg,~,P_pos,P_neg,~] = CPM_LOOCV_diagonal_regress(thresh, thresh, all_mats, all_behav, all_cov, corr_type);

thresh_results(position,1) = thresh;
thresh_results(position,2) = R_pos;
thresh_results(position,3) = P_pos;
thresh_results(position,4) = R_neg;
thresh_results(position,5) = P_neg;

end

    try
       [r, ind] = max(thresh_results(:,2));
       thresh_pos = [thresh_results(ind,1),r,thresh_results(ind,3)];
       [r, ind] = max(thresh_results(:,4));
       thresh_neg = [thresh_results(ind,1),r,thresh_results(ind,5)];
    catch
        warning('Problem detecting optimal thresholds. Assign a value of 0');
        thresh_pos= [0,0,0];
        thresh_neg= [0,0,0];
    end
    
end
