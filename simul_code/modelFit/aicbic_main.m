% AICBIC
clear;
clc;
load('SBJ_structure_each_exp_BESTcollectionJan29.mat');

AIC = zeros(22,1);
BIC = zeros(22,1);
NLL = zeros(22,1);
NUM = zeros(22,1);
for i = 1 : 22
    NUM(i) = SBJ{1, i}.num_data;
    NLL(i) = SBJ{1,i}.model_BayesArb.val;
    AIC(i) = 6*2 + 2 * SBJ{1,i}.model_BayesArb.val;
    BIC(i) = SBJ{1,i}.model_BayesArb.val * 2 + 6*log(SBJ{1, i}.num_data);  
end

% wt= subject number
% val_wt = value that you want to compute
subject_NUM= 2
val_wt = 189;

wt=subject_NUM;



BICori = BIC(wt)
BICnew = val_wt * 2 + 4 * log(SBJ{1, wt}.num_data)