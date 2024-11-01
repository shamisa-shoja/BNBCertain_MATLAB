function [options,default_opts] = miqpsettings

% settings for MILP and MIQP
%_________________________________________________________________________
options.epsi = 1e-4;                % tolerance for zero

options.lowmemory = true;           % false: get many information of BnB
options.find_max = true;            % find point with maximum iteration number in offline case
options.find_sol_off = false;       % find solution and value function at each sample point
options.int_tol = 1e-4;             % integer tolerance, use in online case, to round small number
options.check_cut_cond = true;      % to enable cut conditions, in flase: enumerate all nodes
options.search_type= 'lifo';        % Strategy:  Last in- First Out   
options.compareQP = true;           % if compare QPs in offline case which is time consuming
options.solver_QPcompare ='gurobi'; % solver to compare QPs: 'bmibnb', 'gurobi', 'cplex'
                                    % cplex: faster than gur but less reliable. 'gurobi' is precise than cplex,but slower;
                                    %'bmibnb' is a bit more precise than gurobi but gurobi is so faster
options.all_compareQPs = false;     % compare QP with all 3solvers: bmibnb,gurobi,cplex
options.reduceOverlap = false;      % reduce overlao in final integer solution (upper bound)
options.detail_online = true;       % get all info of online solver
options.onlineMiqp_trad = false;    % apply tradition online_miqp code as well
options.reduceOverlap_eachstep = false;    % reduce overlap at each step after solving QP certification
options.save_both_Red_notRed_sol = false;  % save both reduceOverlap and original integer solutions
options.enum_check_inf = true;             % enumeratin case, cut whent node is infeasible
options.store_all_upperbound = false;      % store all upper bounds
options.compare_upp_low = false;           % compare QPs in just upper bounds. True: compare lower bound with upperbound
options.final_red_overlap = false;         % final overlap removal of explicit solution in bnb
options.analyze2online = false;            % find diffs in offline and online results
options.check_dominance_again = false;     % check cut condition once more if reduceOverlap says to cut the node
options.gen_rand_samples = true;           % generate some random sample points and apply online miqp in this points
options.gen_chebyCenters = false;          % find chebychev centers and apply online miqp in this points
options.existed_point = false;             % the existed points for the online analysis
options.apply_red_overlap_func = false;    % Algorithm 3: use already existed reduceOverlap function in QPcomparision
options.try_enum_miqp = false;             % apply enumeration in BnB
options.QPcompare_debug = false;           % to do minor debugging in QP comparision algoritms
options.disregard_lowerdim_reg = true;     % disregard lower dimentional regions that mpQP generates 
options.plot_innerQP_region = true;        % plot inner regions after certifying QP
options.append_to_S = true;                % method for pushing to S: append or different search strategies
options.push_proc_to_S = true;             % to push to S in the inserting procedure  
options.store_all_upperbound = true;

options.warm_start = false;                % for warm-starting, default is cold start
options.MILP = false;                      % true for MILP, false in MIQP
options.computeAS1 = false; % use modified AS0 for p1 as well (not compute AS1)
options.AS0toAS1 = true; % when AS_parent for p1 is empty, consider AS_paret of p0

% certified bnb_miqp options, SS
options.algorithm1 = false;         % Algorithm 1: solve certified mpMIQP problem: use recent int_sol as upperbound
options.algorithm2 = true;          % Algorithm 2: solve certified mpMIQP problem: Store all upper bounds locally
options.algorithm3 = false;         % Algorithm 3: solve certified mpMIQP problem: Store all upper bounds globally

%_________________________________________________________________________
% tradition bnb_miqp options
%-------------------------- Options ---------------------------------------
options.solver     = 'Explicit'; % specify here the default solver for QPs
options.method     = 'depth';    % specify here the default method
options.branchrule = 'first';    % specify here the default branchrule
options.order      = 0;          % specify here the default order
options.verbose    = 1;%0        % specify here the default verbose
options.maxqp      = inf;        % specify here the default maxqp
options.inftol     = 1e8;        % specify here the default inftol
options.matrixtol  = 1e-6;       % specify here the default matrixtol
options.postol     = 1e-6;       % specify here the default postol
options.integtol   = 1e-4;       % specify here the default integtol
options.maxQPiter  = 1000;       % specify here the default maxQPiter
options.solverlp   = 'dsimplex'; %'linprog';  % specify here the default solver for LPs
options.deletefile = 1;          % specify here the default deletefile-flag
options.round      = 1;          % specify here the default round-flag
options.xxmax      = 1e7;        % specify here the default xmax for qp_dantz
options.eps_zero   = 1e-10;

% #____________Suboptimality settings___________
  options.subopt_sol = false; 
  %___epsilont-suboptimality___
  options.subopt_epsilon = false; %relax dom cut
  options.tol_domcut_abs = 1e-6;
  options.tol_domcut_rel = 1e-6;
  options.tol_abs_param = false; %absoulte allowance parameter is parametric?
  options.tol_domcut_abs_param = [1e-6.*eye(2) 1e-6.*ones(2,1)]; %eps = aJ: [Beps Ceps],dim: 1*(nth + 1) for MILP || [Aeps Beps' Ceps], dim:(nth)*(nth + 2) for MIQ
  %___T-cut method___
  options.subopt_Tcut = false; % %generated-node cut
  options.bound_tot_gen_nodes = 10000; % bounds on total generated nodes
  options.bound_decomp_nodes = 10000; % bounds on total generated nodes
  %___M-cut method___
  options.subopt_Mcut = false; % %active-node cut
  options.bound_active_nodes = 10000; % bounds on stored active nodes
  options.find_subopt_err_MILP = true; % find suoptimal error by solving an MILP in each region, If false:solve relaxation instead
%__________________________________________________________________________
%---------------------- Default Options -----------------------------------

% specify defaults for all undefined parameters in Options
default_opts.solver     = 'quadprog'; % specify here the default solver for QPs
default_opts.method     = 'depth';    % specify here the default method
default_opts.branchrule = 'first';    % specify here the default branchrule
default_opts.mp_solver = 'MPT3';
default_opts.order      = 0;          % specify here the default order
default_opts.verbose    = 0;          % specify here the default verbose
default_opts.maxqp      = inf;        % specify here the default maxqp
default_opts.inftol     = 1e8;        % specify here the default inftol
default_opts.matrixtol  = 1e-6;       % specify here the default matrixtol
default_opts.postol     = 1e-6;       % specify here the default postol
default_opts.integtol   = 1e-4;       % specify here the default integtol
default_opts.maxQPiter  = 1000;       % specify here the default maxQPiter
default_opts.solverlp   = 'linprog';  % specify here the default solver for LPs
default_opts.deletefile = 1;          % specify here the default deletefile-flag
default_opts.round      = 1;          % specify here the default round-flag
default_opts.xxmax      = 1e7;        % specify here the default xmax for qp_dantz


end