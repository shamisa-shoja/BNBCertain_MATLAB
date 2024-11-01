% certification of MILP/MIQP

%% upload required folders
upload_folders; % in main folder
%%
clc; clear; close all;
yalmip('clear')

%% Load system
N_ex = 6; % Random MIQP Example
[prob, P_theta, AS0, x0, vartype] = example_MILP(N_ex);

%% Load settings
[options,opts_default] = miqpsettings;
opts_cert = certsettings;

%% Depth-first Strategy for MIQP
%_____________
options.method = 'depth';
tic;
[final_part_QP, expl_sol_QP, bnb_info_off_QP,~] = Bnb_mpMILP(prob,P_theta,AS0,x0,options,opts_default,opts_cert); 
time_MILP_QP = toc;

%% Online MIP for depth first
%%__________Analyze online results______________
options.gen_chebyCenters = false; 
options.gen_rand_samples = true;  
N_sample = 10;                   % Number of samples in parameter set
[detail_sol_QP] = MIP_online(prob,P_theta,AS0,x0,final_part_QP,expl_sol_QP,bnb_info_off_QP,options,opts_default,opts_cert,N_sample);

