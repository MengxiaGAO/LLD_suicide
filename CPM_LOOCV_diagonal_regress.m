function [R_pos,R_neg,R_both,P_pos,P_neg,P_both] = CPM_LOOCV_diagonal_regress(thresh_pos, thresh_neg, all_mats, all_behav, all_cov, corr_type)
%% This function is to use the optimal threshold and select the features in the positive and negative networks (no normalization)
% Input:
%       thresh_pos:     the optimal p threshold of positive network (e.g., thresh_pos(1,1))
%       thresh_neg:     the optimal p threshold of negative network (e.g., thresh_neg(1,1))
%       all_mats:       an X*Y*Z 3D matrix
%       all_behav:      a Z*1 vector
%       all_cov:        a Z*M vector contains different covariables (can be empty)
%       corr_type:      correlation type: 'Pearson', 'Kendall', 'Spearman'

% Output:
%       R_pos:  r value of the positive network
%       R_neg:  r value of the negative network
%       R_both: r value of the combined network
%       P_pos:  p value of the positive network
%       P_neg:  p value of the negative network
%       P_both: p value of the combined network

%%%%% The script is adopted from Xinlin Shen (2017, Nature Protocols). 
%%%%% Written by Mengxia Gao, PhD student of Dept. of Psychology, the University of Hong Kong.
%%%%% mengxia.gao@gmail.com, 20200310
%%%%% updated on 20201022

no_sub = size(all_mats,3);
no_node1 = size(all_mats,1);
no_node2 = size(all_mats,2);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred_both = zeros(no_sub,1);

network = ones(no_node1,no_node2);
net_upp = triu(network,1);
upp_id = find(net_upp);
no_edge = length(upp_id);

for leftout = 1:no_sub
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    %-------(Step 2) Divide data into training and testing sets
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats, no_node1*no_node2, no_sub-1)';
    train_vcts_upp = train_vcts(:,upp_id);

    train_behav = all_behav;
    train_behav(leftout) = [];
    
    train_cov = all_cov;
    train_cov(leftout,:) = [];
    
    %-------(Step 3) Relate connectivity to behavior
    % correlate all edges with behavior

    [r_mat, p_mat] = partialcorr(train_vcts_upp, train_behav, train_cov, 'type', corr_type);

    pos_edges = zeros(1, no_edge);
    neg_edges = zeros(1, no_edge);
    
    pos_mask = zeros(no_node1, no_node2);
    neg_mask = zeros(no_node1, no_node2);
    %-------(Step 4) Edge selection
    % set threshold and define masks
    
    pos_edges(r_mat > 0 & p_mat < thresh_pos) = 1;
    neg_edges(r_mat < 0 & p_mat < thresh_neg) = 1;
    
    pos_mask(upp_id) = pos_edges;
    neg_mask(upp_id) = neg_edges;
    
    %-------(Step 5) Calculate single-subject summary values
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:no_sub-1
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask, 'omitnan'));
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask, 'omitnan'));
    end
    train_sum_all = train_sumpos - train_sumneg;
    
    %-------(Step 6) Model fitting
    % build model on TRAIN subs

    fit_pos = regress(train_behav, [ones(no_sub-1,1), train_sumpos]);
    fit_neg = regress(train_behav, [ones(no_sub-1,1), train_sumneg]);
    fit_both = regress(train_behav, [ones(no_sub-1,1), train_sum_all]);
      
    % run model on TEST sub
    test_sumpos = sum(sum(all_mats(:,:,leftout).*pos_mask, 'omitnan'));
    test_sumneg = sum(sum(all_mats(:,:,leftout).*neg_mask, 'omitnan'));
    test_sum_all = test_sumpos - test_sumneg;
    
    behav_pred_pos(leftout) = fit_pos(1) + fit_pos(2)*test_sumpos;
    behav_pred_neg(leftout) = fit_neg(1) + fit_neg(2)*test_sumneg;
    behav_pred_both(leftout) = fit_both(1)+ fit_both(2)*test_sum_all(1);
    
end

%------(Step 7) Prediction in novel subjects
% compare predicted and observed scores
[R_pos, P_pos] = corr(behav_pred_pos,all_behav,'type', corr_type);
[R_neg, P_neg] = corr(behav_pred_neg,all_behav,'type', corr_type);
[R_both, P_both] = corr(behav_pred_both,all_behav,'type', corr_type);

end