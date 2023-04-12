clear all
close all
clc

% Script for obtaining the result which are gonna be used in the CLHU paper

addpath('methods')

input = load("test\input.mat");
endmember_estimation = load("test\endmember_estimation.mat");

Xim = input.X;
X = transpose(reshape(Xim, [input.nRow*input.nCol, input.nBand]));
n_endmembers = input.n_endmembers;

vca_ee = transpose(endmember_estimation.VCA);
nfindr_ee = transpose(endmember_estimation.NFINDR);

%% FCLS
[A_FCLS,M_FCLS,time_fcls,Yhat_FCLS] = adaptor_FCLS(Xim,vca_ee,[],[]);
A_VCA = A_FCLS;
Xhat_VCA = Yhat_FCLS;
%A_VCA = reshape(transpose(A_FCLS), [input.nRow, input.nCol, n_endmembers]);

A_init = A_FCLS; % with VCA its better

[A_FCLS,M_FCLS,time_fcls,Yhat_FCLS] = adaptor_FCLS(Xim,nfindr_ee,[],[]);
%A_NFINDR = reshape(transpose(A_FCLS), [input.nRow, input.nCol, n_endmembers]);
A_NFINDR = A_FCLS;
Xhat_NFINDR = Yhat_FCLS;

save("test\FCLS.mat", 'A_NFINDR', 'Xhat_NFINDR', 'A_VCA', 'Xhat_VCA')
%% ELMM
% spectral library extraction
bundle_nbr = 5; % number of VCA runs
percent = 20; % percentage of pixels considered in each run
Lib_vca = extractbundles_batchVCA(X, vca_ee, bundle_nbr, percent);
Lib_nfindr = extractbundles_batchVCA(X, nfindr_ee, bundle_nbr, percent); 

opt_elmm.lambda_s = 1;
opt_elmm.lambda_a = 0.05;
opt_elmm.lambda_psi = 0.01;

% VCA + ELMM
[A_ELMM,M_ELMM,time_elmm,Yhat_ELMM] = adaptor_ELMM(Xim,vca_ee,Lib_vca,A_init,opt_elmm);
A_VCA = A_ELMM;
M_VCA = M_ELMM;
Xhat_VCA = Yhat_ELMM;

% NFINDR + ELMM
[A_ELMM,M_ELMM,time_elmm,Yhat_ELMM] = adaptor_ELMM(Xim,nfindr_ee,Lib_nfindr,A_init,opt_elmm);
A_NFINDR = A_ELMM;
M_NFINDR = M_ELMM;
Xhat_NFINDR = Yhat_ELMM;

save('test/ELMM.mat', 'A_NFINDR', 'M_NFINDR', "Xhat_NFINDR", 'A_VCA', 'M_VCA', 'Xhat_VCA')