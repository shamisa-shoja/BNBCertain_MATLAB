%online mpMIQP

function [iter_on,info_on,iter_on2,info_on2] = compute_online_miqp(mpQP,th_sample,options,default_opts,cert_opts)

options.solver = 'cert_qp';
Nx = size(mpQP.H,1); 
x0_1 = zeros(Nx,1);
%__________________________________________________________________________
%initialzie variables to store data
nr_sample = size(th_sample,2);
x_online2 = zeros(Nx,nr_sample);      J_online2 = zeros(nr_sample,1); 
iter_on2 = zeros(nr_sample,1);        sol_track =cell(0);
x_online = zeros(Nx,nr_sample);       J_online = zeros(nr_sample,1);
iter_on = zeros(nr_sample,1);         bnb_info = cell(0);

info_on = cell(0);                    info_on2 = cell(0);
seq_cut_cond = [];
%__________________________________________________________________________
% define problem formulation
subprob = mpQP;
H = subprob.H; 
A = subprob.A;
Aeq = []; beq = []; %Aeq = subprob.Aeq; beq = subprob.beq;
lb = subprob.lb;   ub = subprob.ub;
vartype = subprob.vartype;

% Solve online miqp for each parameter
%__________________________________________________________________________
for k = 1 : nr_sample
    
    f = subprob.f + subprob.f_theta * th_sample(:,k);
    b = subprob.b + subprob.W * th_sample(:,k);
    %______________________________________________________________________
    %written online certified miqp
    [xmin, Jmin, iter_num, bnb_info_on] = onlineMIQP_exp(H, f, A, b, Aeq, beq, vartype, lb, ub, x0_1, options,default_opts,cert_opts);
    x_online(:,k) = xmin;
    J_online(k,1) = Jmin;
    iter_on(k,1) = iter_num;
    bnb_info{k} = bnb_info_on;
    %seq_cut_cond = [seq_cut_cond ; bnb_info_on.seq_cut_cond];
    %__________________________________________________________________
    %original miqp with certified QP solver
    if options.onlineMiqp_trad
        [xmin2, Jmin2, iter_num2, sol, ~, ~] = online_miqp(H, f, A, b, Aeq, beq, vartype, lb, ub, x0_1, options); %flag, Extendedflag
        x_online2(:,k) = xmin2;
        J_online2(k,1) = Jmin2;
        sol_track{k} = sol; 
        iter_on2(k,1) = iter_num2;
    end
    
end
%__________________________________________________________________________
% final results in detail 
if options.detail_online
    info_on.x = x_online;
    info_on.J = J_online;
    info_on.bnb_info = bnb_info;
    %info_on.seq_cut_cond = seq_cut_cond; %track cut conditions
    
    
    if options.onlineMiqp_trad
        info_on2.x = x_online2;
        info_on2.J = J_online2;
        info_on2.sol = sol_track;
    end
end
%__________________________________________________________________________
%[max_iter_online, ind_max_iter_online] = max(iter_on);
%[min_iter_online, ind_min_iter_online] = min(iter_on);

end