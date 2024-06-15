%--------------------------------------------------------------------------
% This function takes the groups resulted from spectral clutsering and the
% ground truth to compute the misclassification rate.
% groups: [grp1,grp2,grp3] for three different forms of Spectral Clustering
% s: ground truth vector
% Missrate: 3x1 vector with misclassification rates of three forms of
% spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2010
%--------------------------------------------------------------------------


function [missrate, index] = Misclassification(groups,s)

n = max(s);             %numero di gruppi

[missrate, index] = missclassGroups(groups(:,1), s, n);
missrate = missrate ./ length(s);