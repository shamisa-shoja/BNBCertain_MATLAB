% certification of mpMILP 

function [final_partition, explicit_sol, bnb_info, cert_info] = Bnb_mpMILP(prob,P_theta,AS0,x0,options,default_opts,cert_opts)

tstart = cputime;
%__________________________________________________________________________
%%Setup
  % Define the problem  
  main_prob = prob;
  main_prob.Nx = length(prob.f); Nx = main_prob.Nx;
  [main_prob.N_constr, main_prob.N_theta] = size(prob.W); N_theta = main_prob.N_theta;
  main_prob.Nbin = length(main_prob.vartype); 
      
  if isfield(P_theta,'lb')
       main_prob.bndA = [-eye(main_prob.N_theta); eye(main_prob.N_theta)];
       main_prob.bndb = [-P_theta.lb; P_theta.ub];      
  else
       bnd_Theta = P_theta; bnd_Theta = bnd_Theta.minHRep();
       main_prob.bndA = bnd_Theta.A; 
       main_prob.bndb = bnd_Theta.b; 
  end
   
%     if ~isfield(P_theta,'A')                 
%         P_theta.A = [eye(N_theta); -eye(N_theta)];
%         P_theta.b = [P_theta.ub; -P_theta.lb];        
%     end
    
    %%%%%%%%%%%%%%%%%%%%% MILP %%%%%%%%%%%%%%%%%%%%%
    % In case of LP: define H and f_theta as zero
    if (options.MILP  || ~isfield(main_prob,'H'))                
        main_prob.H = zeros(Nx); 
        options.store_all_upperbound = false; 
    end
    if (options.MILP  || ~isfield(main_prob,'f_theta'))                 
        main_prob.f_theta = zeros(Nx,N_theta);      
    end
    time_compare = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change boundery constraints to inequality constraints in the begenning
    options.boundery_to_ineq = false;
    if options.boundery_to_ineq
        %change boundery constraints to inequalities in the begenning
        main_prob.lb(main_prob.vartype) = zeros(main_prob.Nbin,1);
        main_prob.ub(main_prob.vartype) = ones(main_prob.Nbin,1);
        [main_prob.A,main_prob.b,main_prob.W,main_prob.lb,main_prob.ub,bin_constr,constr_type] = add_bin_cons(main_prob.A,main_prob.b,main_prob.W,main_prob.lb,main_prob.ub,main_prob.vartype);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Do some extra settings
  options = miqpsettings_extra(options,default_opts,main_prob);
  [main_prob,x0] = probsettings(main_prob,x0,options);
  
  %options.check_cut_cond = true; % flase: enumerate all nodes
  
  % define some loop variables
  prob_num = 0;           % counts how many times the QP or LP algorithm is invoked
  optP  = 0;		      % number of QP or LP where optimum is found
  

  % initilazie some variables to store BnB information
  seq_ivalues = [];              %track ivalues
  seq_levels = [];               %track some information of bnb
  seq_N_reg = [];                %track number of regions in each solved relaxed QP
  seq_N_reg_int = [];            %track number of regions with integer sol in each solved relaxed QP
  seq_N_reg_int_ratio = [];
  seq_int_var = [];              % accumulate integer solution
  seq_all_regions = cell(0);     % accumulate all solutions in each solved QP
  seq_int_feas_sol = [];         % save if whole regions are integer feasile
  seq_optP = [];
  seq_int_var_acc = [];          % sequences of integer solution
  seq_feas_whole_node = [];
  seq_P_num = [];
  seq_cut_cond = [];          % to keep track of cut conditions
                              % 0: not cut, 1:infeasibility cut, 2:integer
                              % feasibility, 3: dominance cut, 4: partial infeas cut
                              % in enumeraion case: -2: integer-feasibility
                              % cut(leaf-node), -4: no cut (branch further)
  node_info = [];             % information of each node: [ivalues; -9; level; -9; cut cond.], -9: to distinguesh
                              % cut cond.= 1:infeas, 2:int feas, 3: dominance, 0: branch further(not cut)


  
  % to store certification information at each iteration
  cert_info.cert_info_indetail = cell(0);
  cert_info.part_detail = cell(0);     
  cert_info.part_stats_dual = cell(0); 
 
  % to update upper bound and integer subproblem
  explicit_sol = cell(0); % to save and accumulate integer solution
  compareprob_info = [];    % to store solution of QP comparision from different solvers    
  time_compQP = 0;        % to find max time taken to compare QPs
  seq_time_compQP = [];    % in case of more than 1 QP comparision in each iterate
  time_redOver_upp = 0;       % time taken to reduce overlap at each time step
  upper_bound_redOver = cell(0);   % store reducedOverlap solution
  %%initilazie final certified partition
   final_partition = cell(0);
   
  %node_inform = [];                % store node information :[ivalues; -9; level; -9; cut cond.], -9: to distinguesh
                                   % cut cond.= 1:infeas, 2:int feas, 3: dominance, 0: branch further(not cut)
  seq_ivalues_node = [];           % store integer values in final partition for each node
  F_stack = cell(0);               % Final stack
  upper_bound = cell(0);           % to accumulate upper bound solution (integer feasible sol)
  options.compare_upp_low = true;  % to compare upper bound with relaxed sol at each step
  
   %_________________________________________________________________________
  %%Define starting condition
  % initialzie STACK  
  stack{1}    = struct('H',main_prob.H, 'f',main_prob.f,  'f_theta',main_prob.f_theta,...                             % objective function
                'A',main_prob.A,   'b',main_prob.b,  'W',main_prob.W, 'Aeq',main_prob.Aeq, 'beq',main_prob.beq, ...   % constraints
                'lb',main_prob.lb, 'ub',main_prob.ub,'bndA',main_prob.bndA,'bndb',main_prob.bndb,...                  %lb and ub on variables, bndA,b are on parameters
                'vartype',main_prob.vartype, 'ivalues',-ones(main_prob.Nbin,1), 'level',0, 'x0',x0, 'e',0, 'e_param',zeros(1,N_theta),'xb',-ones(main_prob.Nbin,1),...      %settings for binary variables
                'box_ind' ,main_prob.box_ind, 'upper_ind',main_prob.upper_ind , 'lower_ind',main_prob.lower_ind,...    %for certification purpose
                'AS0', AS0,'x0_case',1,'normal',0,'iter',0,...          %for certification purpose 
                'iter_stack',0,'parent_stack',[],'AS_stack',[],'AS_parent',[],... %for warm-starting purpose          
                'J_parent', inf, 'Theta', P_theta);       % For best-first startegy                                                 
                % 'branching_priority', options.branching_priority, 'explicit_eval', options.explicit_eval,'integer_times', options.integer_times, ... % for bnb purpose

  
  % initilazie number of nodes to be proceed 
  N_nodes_total = 2^(length(main_prob.vartype)+1)-1;
  N_nodes_processed = 0; N_nodes = 0;
  
   % Warm start: save active sets, initilazie it as empty
  stack{1}.AS = [];
  stack{1}.AS_tot = [];
  
  % initialize internal stack T (state of branch and bound)
  T_stack{1} =  stack{1};                          % store the state of tree
  %T_stack_cost = 0;
  
  % initialize spatial information
  S_stack{1}.Theta = P_theta;                       % initialize parameter set                         
  S_stack{1}.iter = 0;                              % initialize complexity number
  S_stack{1}.Tree = stack;                          % initialize internal stack
  J.Ai =zeros(N_theta); J.Bi =zeros(N_theta,1)';    % initialize upper bound
  J.Ci =inf;                                        % initialize upper bound to be inf
  x.Fi = zeros(Nx,N_theta); x.Gi =zeros(Nx,1);      % initialize solution
  if options.store_all_upperbound                   % if storing all upper bounds
      S_stack{1}.J_upper{1} = J;                        
      S_stack{1}.x_upper{1} = x;
  else
      S_stack{1}.J_upper = J;
      S_stack{1}.x_upper = x;
  end
  
  % to find out which node we are in:
  prob_info.level = 0;
  prob_info.ivalue = -ones(main_prob.Nbin,1);
  prob_info.prob_num = 0;
  prob_info.jreg = 0;
  prob_info.allreg = 0;
  S_stack{1}.prob_info = prob_info;
  % to keep track of working set changes
  WS.parent = [];                  % store working set info: parent and active sets
  WS.AS = [];
  WS.AS_logic = [];
  S_stack{1}.WS_info = WS;
  S_stack{1}.exp_int_sol = cell(0);% integer solution in each region
  S_stack{1}.info.node_info = [];       % info of solution
  S_stack{1}.info.ivalues = [];         % info of solution
  S_stack{1}.info.current_node = [-ones(main_prob.Nbin,1);-9;0;-9;0];        %info of cut conditions
                                        %current_node = [subprob.ivalues; -9; subprob.level; -9; 1];
                                        %-9 to distinguesh between values, last num=0: branch
  S_stack{1}.info.current_cut_cond = 'root node';
  S_stack{1}.info.level = 0;
  S_stack{1}.info.cut_cond = [];        %info of cut conditions
  S_stack{1}.info.sizeTree = 0;         %info of size of the tree (number of nodes)

  % Warm start: save active sets, initilazie it as empty
  S_stack{1}.AS_par = [];
  AS_parent = []; %warm start
  
  S_stack{1}.J_lower = inf;  % For Best first strategy
  S_stack_out = cell(0);
  info_copy= cell(0); 
% =========================================================================
%__________________________________________________________________________
% Main Loop
  starttime = cputime;

  k = 0; loop_iter = 2; % test loop
  % while (k < loop_iter)
  
  while(~isempty(S_stack))
	k = k + 1;	
    
    if (prob_num > options.maxqp)|| (cputime-starttime > options.maxtime) %||(flag == -1) 
	  if(cert_opts.break_if_above)
		S_stack = [];
		final_partition= cell(0);
        final_partition{1}.Region = [];
		final_partition{1}.iter = nan;
		continue	
	   end
	  cert_info.unfin_regs = cert_info.unfin_regs+1;
	  continue
    end
    
    %______________________________________________________________________
     % pop BnB state (T_stack) from spatial stack (S_stack)
     % [S_stack,T_stack,Theta,iter,J_upper,~] = pop_BnB(S_stack,cert_opts.search_type,'spatial',[]); 
     [S_stack,S_tuple] = pop(S_stack,cert_opts.search_type,[]);
     T_stack = S_tuple.Tree;
     Theta = S_tuple.Theta;
     iter = S_tuple.iter;
     J_upper = S_tuple.J_upper;
     x_upper = S_tuple.x_upper; 
     prob_info = S_tuple.prob_info;
     WS = S_tuple.WS_info;
     exp_int_sol = S_tuple.exp_int_sol;
     info = S_tuple.info; 
     AS_par = S_tuple.AS_par; % for warm start
     
    T_stack_cost = 0;
    
    if ~isempty(T_stack)
        % Get the next subproblem from tree (T_stack)
        [T_stack,subprob,T_stack_cost] = pop(T_stack, cert_opts.search_type,T_stack_cost); 

    %__________________________________________________________________________
        % Solve relaxed problem
        try
            if options.MILP
                [all_regions,N_regions, int_regions,N_int_regions, non_int_regions, feasible, int_feas_sol, int_sol_ind, part,cert_info_,~,~] = ...
                    LP_cert(subprob,main_prob,Theta,cert_opts,seq_int_var_acc,AS_parent,options);  
            else
                %[all_regions,N_regions, int_regions,N_int_regions, non_int_regions, feasible, int_feas_sol, int_sol_ind, part,cert_info_,~,~] = ...
                %   QPcertification_WS2(subprob,main_prob,Theta,cert_opts,seq_int_var_acc,subprob.AS,options); %expl_int_seq_sols,seq_int_var_acc
                [all_regions,N_regions, int_regions,N_int_regions, non_int_regions, feasible, int_feas_sol, int_sol_ind, part,cert_info_,~,~] = ...
                    QPcertification(subprob,main_prob,Theta,cert_opts,seq_int_var_acc,options); %expl_int_seq_sols,seq_int_var_acc
            end
        catch
            keyboard
        end
        prob_num = prob_num + 1; 
     %__________________________________________________________________________
        % expand dimension of sobproblem to main dimension
        if feasible        % all or at least some parts of parameter set are feasible

            if(options.verbose > 0)
               fprintf(['#%d | prob_num = %d | level = %d | iter = %d | int/all regions : %d/%d | T_Stack: %d | S_Stack: %d\n'],k,prob_num,subprob.level,iter,N_int_regions,N_regions,length(T_stack),length(S_stack))  %
            end
            % Store some useful info of certification            
            if ~isempty(cert_info_)
                cert_info.cert_info_indetail{end+1} = cert_info_;  
                cert_info.part_detail{end+1} = part;     
                cert_info.part_stats_dual{end+1} = cert_info_.part_stats_dual; 
            end

            %__________________________________________________________________________
            % Store some useful info of BnB
            seq_N_reg = [seq_N_reg; subprob.level N_regions];
            seq_N_reg_int = [seq_N_reg_int; subprob.level N_int_regions];
            seq_N_reg_int_ratio = [seq_N_reg_int_ratio; subprob.level (N_int_regions/N_regions)];
            seq_int_var = [seq_int_var, seq_int_var_acc];
            seq_int_feas_sol = [seq_int_feas_sol; subprob.level int_feas_sol];
            seq_P_num = [seq_P_num;  subprob.level prob_num];
            %__________________________________________________________________________
            % expand reduced subprob dimention to real dimenstion                    
            all_regions = subsol2entiresol_sh(main_prob.H,main_prob.f,main_prob.f_theta,all_regions,main_prob.vartype,subprob.ivalues);
            
            seq_all_regions{end+1} = all_regions; % store all solutions

             if(N_int_regions > 0)
                 int_regions = subsol2entiresol_sh(main_prob.H,main_prob.f,main_prob.f_theta,int_regions,main_prob.vartype,subprob.ivalues); 
                 int_sol_exist = true;                                     

                 if(options.add_complete_int_seq_sols)
                        for j = 1:length(expl_int_seq_sols)
                            expl_int_seq_sols{j} = subsol2entiresol_sh(main_prob.H,main_prob.f,main_prob.f_theta,expl_int_seq_sols{j},main_prob.vartype,subprob.ivalues); 
                        end
                  end
             else
                    int_sol_exist = false;
             end

              if( N_regions - N_int_regions > 0)
                  non_int_regions = subsol2entiresol_sh(main_prob.H,main_prob.f,main_prob.f_theta,non_int_regions,main_prob.vartype,subprob.ivalues);
              else  
                  non_int_regions =[]; 
              end
                            
             %__________________________________________________________________________
        else % whole node is infeasible
            
            % Since we are here, the QP was feasible, and hence, also the parametric solution should be feasible.
            %warning('check_for_possible_int_sols: Numerical problems, or low dimentional region');
            
            all_regions =[]; 
            int_regions =[]; int_sol_exist = false;
            non_int_regions = [];
            
            seq_N_reg = [seq_N_reg; subprob.level -inf];
            seq_N_reg_int = [seq_N_reg_int; subprob.level -inf];
            seq_N_reg_int_ratio = [seq_N_reg_int_ratio; subprob.level NaN];
            
            % Infeasible region: push to S like dominance cut
            disp('whole infeasibility cut');
            infeas_cut_whole = true;
            
            % update info of the node
            seq_cut_cond= [seq_cut_cond; subprob.level 4];
            node = [subprob.ivalues; -9; subprob.level; -9; 4];%-9 to distinguesh between values, last num=4: whole infeas
            node_info = [node_info, node];
            info_copy = cell(0);
            info_copy = info;
            %info_copy.node_info = node_info;    % info of solution
            info_copy.node_info = [info_copy.node_info, node];    % info of solution
            info_copy.ivalues = [info_copy.ivalues, subprob.ivalues];% info of solution
            info_copy.current_node = node;      % info of solution
            info_copy.current_cut_cond = 'whole infeas';
            info_copy.level = [info_copy.level; subprob.level];
            info_copy.cut_cond = [info_copy.cut_cond; subprob.level 4];
            info_copy.sizeTree = info_copy.sizeTree + 1;
          %__________________________________________________________________________
            % update S_stack by pushing new tree to 
            new_S = [];
            new_S.Theta = Theta;
            new_S.iter = iter + 1; %1: if number of nodes %all_regions.iter; %exp_sol_part{jreg}.iter; 
            new_S.Tree = T_stack;         % update tree
            new_S.J_upper = J_upper;     % update J_upper!
            new_S.x_upper = x_upper;
            new_S.exp_int_sol = [];
            new_S.info = info_copy;
            prob_info = [];
            prob_info.level = subprob.level;
            prob_info.ivalue = subprob.ivalues;
            prob_info.prob_num = prob_num;
            prob_info.jreg = 0;
            prob_info.allreg = N_regions;
            new_S.prob_info = prob_info;
            %%store parent and active set info
            %%precise but time consuming
            %%[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
            %if ~iscell(WS.parent)
            %    WS.parent = {WS.parent};
            %end
            %WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
            if ~iscell(WS.AS)
                WS.AS = {WS.AS};
            end
            %WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
            if ~iscell(WS.AS_logic)
                WS.AS_logic = {WS.AS_logic};
            end
            %WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
            new_S.WS_info = WS;
            %new_S.exp_int_sol = exp_int_sol;

            new_S.AS_par = AS_parent;  % WARM start
            
            new_S.J_lower.Ci = inf; %exp_sol_part{jreg}.J_lower; % For Best first strategy
                        
            if options.append_to_S %usual appending to S
                 S_stack = push(S_stack,new_S);
              else 
                 switch options.method
                    case 'depth'
                        S_stack = push(S_stack,new_S);
                    case 'breadth'
                        S_stack = push(S_stack,new_S); % FIX it
                    case 'best'
                        S_stack = push_BF(S_stack,new_S);
                    case 'bestdepth'
                        S_stack = push(S_stack,new_S); % FIX it
                 end
             end         

        end %feasible
        %__________________________________________________________________
        % update some info of BnB
        seq_feas_whole_node = [seq_feas_whole_node; subprob.level feasible];
        seq_ivalues = [seq_ivalues, subprob.ivalues];
        seq_levels = [seq_levels; subprob.level];
        seq_ivalues_node = [seq_ivalues_node, subprob.ivalues]; %store in each final region

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________
        exp_sol = all_regions;
        leaf_node =  ~sum(subprob.ivalues==-1); %all(subprob.ivalues ~= -1); 
        J_lower = cell(0); %initialize lower bound
        
        
     if feasible % ~isempty(exp_sol) % whole node is feasible in some parts
         
        exp_sol_part = reg2part(exp_sol);

        for jreg = 1 : length(exp_sol_part)
            
            T_stack_copy = cell(0);
            T_stack_copy = T_stack;         % copy internal BnB stack 
            
            J_upper_copy = cell(0);
            J_upper_copy = J_upper;         % copy internal upper bound
            
            x_upper_copy = cell(0);
            x_upper_copy = x_upper;         % copy internal upper bound
            
            exp_int_sol_copy = cell(0);     % clean it up from previous iteration
            exp_int_sol_copy = exp_int_sol; % copy internal integer feasible solution 
            
            % info_copy = cell(0);            % clean it up from previous iteration
            % info_copy = info;               % copy internal integer feasible solution 

            exp_sol_part{jreg}.Pfinal = exp_sol_part{jreg}.Region; 
            %__________________________________________________________________________
            %%WARM start: update subproblems active set
            % subprob.AS =  exp_sol_part{jreg}.AS;
            AS_copy = [];                   % clean it up from previous iteration
            AS_copy = AS_par;               % copy internal active set 
            %__________________________________________________________________________            
             % update lower bound (solution of relaxed subprob)
            if exp_sol_part{jreg}.state > 2
                %subproblem in this region is infeasible
                J_lower{jreg}.Ai = exp_sol_part{jreg}.Ai;
                J_lower{jreg}.Bi = exp_sol_part{jreg}.Bi;
                J_lower{jreg}.Ci = inf;                   %infeasible
                infeas_cut = true;                
            else
                J_lower{jreg}.Ai = exp_sol_part{jreg}.Ai;
                J_lower{jreg}.Bi = exp_sol_part{jreg}.Bi;
                J_lower{jreg}.Ci = exp_sol_part{jreg}.Ci;  %feasible
                infeas_cut = false;
            end
                x_lower{jreg}.Fi = exp_sol_part{jreg}.Fi;
                x_lower{jreg}.Gi = exp_sol_part{jreg}.Gi;
        %__________________________________________________________________________
        %                      Check 3 Cut Conditions
        %_________________________________________________________________________           
          J_upper_i = cell(0);  x_upper_i = cell(0);
%           %_____________________________________________________________
%           %check infeasibility
%           if infeas_cut
%             seq_cut_cond= [seq_cut_cond; subprob.level 1];
%             node = [subprob.ivalues; -9; subprob.level; -9; 1];%-9 to distinguesh between values, last num=1: infeas
%             node_info = [node_info, node];
%           end
         %_____________________________________________________________
          % check dominance cut condition 
          dominance_cut = false;
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%% MILP %%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if options.MILP
              
             if infeas_cut
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % infeasibility cut
                disp('infeas cut');
                % update info of the node
                seq_cut_cond= [seq_cut_cond; subprob.level 1];
                node = [subprob.ivalues; -9; subprob.level; -9; 1];%-9 to distinguesh between values, last num=4: whole infeas
                node_info = [node_info, node];                
                info_copy = cell(0);
                info_copy = info;
                info_copy.node_info = [info_copy.node_info, node];    % node_info; 
                info_copy.ivalues = [info_copy.ivalues, subprob.ivalues];
                info_copy.current_node = node;     
                info_copy.current_cut_cond = 'infeas cut';
                info_copy.level = [info_copy.level; subprob.level];
                info_copy.cut_cond = [info_copy.cut_cond; subprob.level 1]; %seq_cut_cond;
                info_copy.sizeTree = info_copy.sizeTree + 1;
                %__________________________________________________________________________
                % update S_stack by pushing new tree to 
                new_S = [];
                new_S.Theta = exp_sol_part{jreg}.Region;
                new_S.iter = iter + exp_sol_part{jreg}.iter; 
                new_S.Tree = T_stack_copy;         % update tree
                new_S.J_upper = J_upper_copy;     % update J_upper!
                new_S.x_upper = x_upper_copy;
                new_S.exp_int_sol = exp_sol_part{jreg}; %exp_int_sol_copy;
                new_S.info = info_copy;
                prob_info = [];
                prob_info.level = subprob.level;
                prob_info.ivalue = subprob.ivalues;
                prob_info.prob_num = prob_num;
                prob_info.jreg = jreg;
                prob_info.allreg = N_regions;
                new_S.prob_info = prob_info;
                if ~iscell(WS.parent)
                    WS.parent = {WS.parent};
                end
                WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
                if ~iscell(WS.AS)
                    WS.AS = {WS.AS};
                end
                WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
                if ~iscell(WS.AS_logic)
                    WS.AS_logic = {WS.AS_logic};
                end
                WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
                new_S.WS_info = WS;
                %new_S.exp_int_sol = exp_int_sol;
                new_S.AS_par = AS_parent;  % WARM start
                new_S.J_lower = J_lower{jreg}; %exp_sol_part{jreg}.J_lower; % For Best first strategy

                S_stack = push(S_stack,new_S);

             else
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                            % dominance cut condition
              
              %if ( ~infeas_cut && ~( J_lower{jreg}.Ci - J_upper.Ci == -inf )) %&& ~leaf_node 
                  t1 = cputime;
                  [~,J_diff_min,~,~,Theta1,Theta2,flag_reg] = spatialSeparate(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region, options);
                  t_com = cputime -t1;
                  time_compare = time_compare + t_com;
                  
                  % Theta1: to be pruned, Theta2: continue further
                  %Theta3 = Union([Theta1,Theta2]); figure; Theta3.plot; % for debugging purposes
                  % figure; exp_sol_part{jreg}.Region.plot
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  if strcmp(flag_reg,'partition') || strcmp(flag_reg,'prune')
                      % Partial Region to be pruned: push it to S
                       disp('dominance cut');
                       
                       % update info of the node
                        seq_cut_cond= [seq_cut_cond; subprob.level 3];
                        node = [subprob.ivalues; -9; subprob.level; -9; 3];%-9 to distinguesh between values, last num=4: whole infeas
                        node_info = [node_info, node];
                        info_copy = cell(0);
                        info_copy = info;
                        info_copy.node_info = [info_copy.node_info, node];    
                        info_copy.ivalues = [info_copy.ivalues, subprob.ivalues];
                        info_copy.current_node = node;      
                        info_copy.current_cut_cond = 'dominance cut';
                        info_copy.level = [info_copy.level; subprob.level];
                        info_copy.cut_cond = [info_copy.cut_cond; subprob.level 3];
                        info_copy.sizeTree = info_copy.sizeTree + 1;
                      %__________________________________________________________________________
                        % update S_stack by pushing new tree to 
                        new_S = [];
                        new_S.Theta = Theta1;
                        new_S.iter = iter + exp_sol_part{jreg}.iter; 
                        new_S.Tree = T_stack_copy;         % update tree
                        new_S.J_upper = J_upper_copy;     % update J_upper!
                        new_S.x_upper = x_upper_copy;
                        new_S.exp_int_sol = exp_sol_part{jreg}; %exp_int_sol_copy;
                        new_S.info = info_copy;
                        prob_info = [];
                        prob_info.level = subprob.level;
                        prob_info.ivalue = subprob.ivalues;
                        prob_info.prob_num = prob_num;
                        prob_info.jreg = jreg;
                        prob_info.allreg = N_regions;
                        new_S.prob_info = prob_info;
                        %store parent and active set info
                        % precise but time consuming
                        %[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
                        if ~iscell(WS.parent)
                            WS.parent = {WS.parent};
                        end
                        WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
                        if ~iscell(WS.AS)
                            WS.AS = {WS.AS};
                        end
                        WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
                        if ~iscell(WS.AS_logic)
                            WS.AS_logic = {WS.AS_logic};
                        end
                        WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
                        new_S.WS_info = WS;

                        new_S.AS_par = AS_parent;  % WARM start
                        new_S.J_lower = J_lower{jreg}; %exp_sol_part{jreg}.J_lower; % For Best first strategy
                        
                        S_stack = push(S_stack,new_S);
                        
                  end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  if strcmp(flag_reg,'keep') || strcmp(flag_reg,'partition')
                      % Partial Region to be considered further

                    if (int_sol_ind(jreg) || leaf_node || subprob.level == length(main_prob.vartype)) % second and third conditions are equal
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        % integer feas: subproblem is solved
                        disp('integer feasible cut');
                        
                        % update info of the node
                        seq_cut_cond= [seq_cut_cond; subprob.level 2];
                        node = [subprob.ivalues; -9; subprob.level; -9; 2];%-9 to distinguesh between values, last num=4: whole infeas
                        node_info = [node_info, node];
                        info_copy = cell(0);
                        info_copy = info;
                        info_copy.node_info = [info_copy.node_info, node];    
                        info_copy.ivalues = [info_copy.ivalues, subprob.ivalues];
                        info_copy.current_node = node;      
                        info_copy.current_cut_cond = 'int feas';
                        info_copy.level = [info_copy.level; subprob.level];
                        info_copy.cut_cond = [info_copy.cut_cond; subprob.level 2]; 
                        info_copy.sizeTree = info_copy.sizeTree + 1;

                        J_upper_i = J_lower{jreg}; % update (local) upper bound
                        x_upper_i = x_lower{jreg}; % update (local) upper bound
                        if ~iscell(exp_int_sol)
                            exp_int_sol = {exp_int_sol};
                        end
                        exp_int_sol{end+1} = exp_sol_part{jreg};
                        
                        %__________________________________________________________________________
                        % final check that J_upper_i isn't empty
                        if isempty(J_upper_i)
                            %disp('J_upper_i is empty!')
                            J_upper_i = J_upper;
                            x_upper_i = x_upper;
                        end
                        %__________________________________________________________________________
                        % update S_stack by pushing new tree to 
                        new_S = [];
                        new_S.Theta = Theta2;
                        new_S.iter = iter + exp_sol_part{jreg}.iter; 
                        new_S.Tree = T_stack_copy;     % update tree
                        new_S.J_upper = J_upper_i;     % update J_upper!
                        new_S.x_upper = x_upper_i;
                        new_S.exp_int_sol = exp_int_sol;
                        new_S.info = info_copy;
                        prob_info = [];
                        prob_info.level = subprob.level;
                        prob_info.ivalue = subprob.ivalues;
                        prob_info.prob_num = prob_num;
                        prob_info.jreg = jreg;
                        prob_info.allreg = N_regions;
                        new_S.prob_info = prob_info;
                        %store parent and active set info
                        % precise but time consuming
                        %[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
                        if ~iscell(WS.parent)
                            WS.parent = {WS.parent};
                        end
                        WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
                        if ~iscell(WS.AS)
                            WS.AS = {WS.AS};
                        end
                        WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
                        if ~iscell(WS.AS_logic)
                            WS.AS_logic = {WS.AS_logic};
                        end
                        WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
                        new_S.WS_info = WS;
                        %new_S.exp_int_sol = exp_int_sol;

                        new_S.AS_par = AS_parent;  % WARM start

                        new_S.J_lower = J_lower{jreg}; %exp_sol_part{jreg}.J_lower; % For Best first strategy

                        S_stack = push(S_stack,new_S);
                        %__________________________________________________________________________ 
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        % branch further
                        % subproblem is feasible, but not integer feasible, branch further             
                    else

                        %__________________________________________________________________________
                        % split subproblem to two new subproblems: P0 and P1
                        disp(['Splitting feasible space (' num2str(jreg) '/' num2str(N_regions) ')']);                    
                        % branchvar:position of branching integer variable within the set of integer variables in present subprob  
                        branchvar            = 1; %decision(x, subprob.vartype, options.branchrule);
                        % separate subproblem
                        options.update_form = false; % do not change problem formulation for subproblems
                        if options.update_form
                            [p0,p1,zeroOK,oneOK] = separate(subprob, branchvar,options);  
                        else
                            [p0,p1,zeroOK,oneOK] = separate_in2e(subprob, branchvar,options,AS_parent);
                        end

                        % push two subproblems to T_stack
                        cost = inf;
                        switch options.method
                            case 'depth'
                                cost = 1/(subprob.level+1);
                            case 'breadth'
                                cost = subprob.level+1;
                            case 'best'
                                cost = J_lower{jreg};
                            case 'bestdepth'
                                if isstruct(J_lower{jreg})
                                    cost =cell(0);
                                    constt = (1/(subprob.level+1));
                                    if isfield(J_lower{jreg},'Ai')
                                       cost.Ai = constt.*J_lower{jreg}.Ai; 
                                    end
                                    cost.Bi = constt.*J_lower{jreg}.Bi;
                                    cost.Ci = constt.*J_lower{jreg}.Ci;
                                else
                                    cost = (1/(subprob.level+1))*J_lower{jreg};  %This privilegiates deep nodes
                                end
                        end
                        % node info for best-first search strategy
                        p0.cost = cost;  
                        p1.cost = cost;
                        p0.Theta = Theta2;
                        p1.Theta = Theta2;
                        
                        % update info of the node
                        seq_cut_cond= [seq_cut_cond; subprob.level 0];
                        node = [subprob.ivalues; -9; subprob.level; -9; 0];%-9 to distinguesh between values, last num=0: branch
                        node_info = [node_info, node];
                        info_copy = cell(0);
                        info_copy = info;
                        info_copy.node_info = [info_copy.node_info, node];    % info of solution
                        info_copy.ivalues = [info_copy.ivalues, subprob.ivalues];% info of solution
                        info_copy.current_node = node;      % info of solution
                        info_copy.current_cut_cond = 'branch';
                        info_copy.level = [info_copy.level; subprob.level];
                        info_copy.cut_cond = [info_copy.cut_cond; subprob.level 2];
                        info_copy.sizeTree = info_copy.sizeTree + 1;
                        
                        if options.push_proc_to_S
                            % update S_stack by pushing new tree to
                            new_S = [];
                            new_S.Theta = Theta2;
                            new_S.iter = iter + exp_sol_part{jreg}.iter; 
                            new_S.Tree = T_stack_copy;     % update tree
                            new_S.J_upper = J_upper_copy;  % not update J_upper
                            new_S.x_upper = x_upper_copy;
                            new_S.exp_int_sol = exp_int_sol;
                            new_S.info = info_copy;
                            prob_info = [];
                            prob_info.level = subprob.level;
                            prob_info.ivalue = subprob.ivalues;
                            prob_info.prob_num = prob_num;
                            prob_info.jreg = jreg;
                            prob_info.allreg = N_regions;
                            new_S.prob_info = prob_info;
                            %store parent and active set info
                            % precise but time consuming
                            %[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
                            if ~iscell(WS.parent)
                                WS.parent = {WS.parent};
                            end
                            WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
                            if ~iscell(WS.AS)
                                WS.AS = {WS.AS};
                            end
                            WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
                            if ~iscell(WS.AS_logic)
                                WS.AS_logic = {WS.AS_logic};
                            end
                            WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
                            new_S.WS_info = WS;
                            %new_S.exp_int_sol = exp_int_sol;

                            new_S.AS_par = AS_parent;  % WARM start
                            new_S.J_lower = J_lower{jreg}; %exp_sol_part{jreg}.J_lower; % For Best first strategy                                 

                            S_stack = push_proc_CS(p0,p1,T_stack_copy,Theta2, new_S,S_stack,options); %all search strategies %[T_stack_copy, S_stack]                            
                            %S_stack = S_stack_out{end};
                            
                        else
                            if options.order == 0
                                if oneOK                                   
                                     T_stack_copy = push_proc(T_stack_copy,p1,options); %all search strategies

                                     %%[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p1,cost);
                                     %%T_stack_copy = push(T_stack_copy,p1);
                                     %
%                                     switch options.method
%                                         case 'depth'
%                                             T_stack_copy = push(T_stack_copy,p0);
%                                         case 'breadth'
%                                             T_stack_copy = push(T_stack_copy,p0); % FIX it
%                                         case 'best'
%                                             S_stack = push_BF(S_stack,new_S);
%                                             T_stack_copy = push_BF(T_stack_copy,p0);
%                                         case 'bestdepth'
%                                             T_stack_copy = push(T_stack_copy,p0); % FIX it
%                                     end
                                end
                                if zeroOK                               
                                    T_stack_copy = push_proc(T_stack_copy, p0, options);
                                end
                            else
                                if zeroOK
                                    %push_BnB(p0,cost);
                                    %[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p0,cost);
                                    T_stack_copy = push_proc(T_stack_copy,p0,options);
                                end
                                if oneOK
                                    %push_BnB(p1,cost);
                                    T_stack_copy = push_proc(T_stack_copy,p1,options);
                                end
                            end  %if options.order 
                                %_________________________________________________________
                                % final check that J_upper_i isn't empty
                                if isempty(J_upper_i)
                                    %disp('J_upper_i is empty!')
                                    J_upper_i = J_upper;
                                    x_upper_i = x_upper;
                                end
                                %__________________________________________________________________________
                                % update S_stack by pushing new tree to 
                                new_S = [];
                                new_S.Theta = Theta2;
                                new_S.iter = iter + exp_sol_part{jreg}.iter; 
                                new_S.Tree = T_stack_copy;     % update tree
                                new_S.J_upper = J_upper_i;     % update J_upper!
                                new_S.x_upper = x_upper_i;
                                new_S.exp_int_sol = exp_int_sol;
                                new_S.info = info_copy;
                                prob_info = [];
                                prob_info.level = subprob.level;
                                prob_info.ivalue = subprob.ivalues;
                                prob_info.prob_num = prob_num;
                                prob_info.jreg = jreg;
                                prob_info.allreg = N_regions;
                                new_S.prob_info = prob_info;
                                %store parent and active set info
                                % precise but time consuming
                                %[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
                                if ~iscell(WS.parent)
                                    WS.parent = {WS.parent};
                                end
                                WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
                                if ~iscell(WS.AS)
                                    WS.AS = {WS.AS};
                                end
                                WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
                                if ~iscell(WS.AS_logic)
                                    WS.AS_logic = {WS.AS_logic};
                                end
                                WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
                                new_S.WS_info = WS;
                                %new_S.exp_int_sol = exp_int_sol;

                                new_S.AS_par = AS_parent;  % WARM start

                                new_S.J_lower = J_lower{jreg}; %exp_sol_part{jreg}.J_lower; % For Best first strategy

                              if options.append_to_S  %usual appending to S
                                 S_stack = push(S_stack,new_S);
                              else 
                                 switch options.method
                                    case 'depth'
                                        S_stack = push(S_stack,new_S);
                                    case 'breadth'
                                        S_stack = push(S_stack,new_S); % FIX it
                                    case 'best'
                                        S_stack = push_BF(S_stack,new_S);
                                    case 'bestdepth'
                                        S_stack = push(S_stack,new_S); % FIX it
                                 end
                              end % if push to S
                       end %push to S

                            if(jreg>=2)
                                N_nodes_total = N_nodes_total + 2*(2^length(subprob.vartype-1+1)-1); % Adding the extra nodes from the _two_ branches (hence 2*...) _below_ the current node (which should not be added, hence, ...-1...).
                            end
                            %______________________________________________________                    
                     end % integer feasibility
                        
                  end %end keep or partition
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              end % cut conditions
          %_______________________________________________________________
        N_nodes = N_nodes + 1;    
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%% MIQP %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            if options.check_cut_cond % check out 3cut conditions
              if options.compareQP
                 if ~iscell(J_upper)
                     %______________________________________________________
                      % do not store several upper bound, use the last upper bound locally
                    if ( ~infeas_cut && ~( J_lower{jreg}.Ci - J_upper.Ci == -inf )) %&& ~leaf_node 
                        %compare QP functions, using either bmibnb, gurobi, cplex!                  
                        switch options.solver_QPcompare
                            case 'bmibnb'
                                t1 = cputime;
                                [~,J_diff_min,~,exitflag_diff] = compareQPs_bmibnb(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region);  %[sol_diff,.,th_min_diff,.] = ..
                                t_bmi = cputime -t1;
                                time_compQP = time_compQP + t_bmi;
                            case 'cplex'
                                t1 = cputime;
                                [~,J_diff_min,~,exitflag_diff] = compareQPs_cplex(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region);
                                t_c = cputime -t1;
                                time_compQP = time_compQP + t_c;                        
                            case 'gurobi'
                                t1 = cputime;
                                [~,J_diff_min,~,exitflag_diff] = compareQPs_gurobi(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region);%Theta
                                t_g = cputime -t1;
                                time_compQP = time_compQP + t_g;
                        end
                        %______________________________________________________
                        % compare QPs sol with other solvers
                        if (strcmp(options.solver_QPcompare,'bmibnb') && options.all_compareQPs )
                            t2 = cputime; 
                            [~,J_diff_min_g,~,exitflag_diff_g] = compareQPs_gurobi(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region);
                            t_g = cputime -t2;
                            %compareprob_info = [compareprob_info, [-1; -2], [J_diff_min exitflag_diff t_bmi; J_diff_min_g exitflag_diff_g t_g]];
                            t3 = cputime; 
                            [~,J_diff_min_c,~,exitflag_diff_c] = compareQPs_cplex(J_lower{jreg}, J_upper, exp_sol_part{jreg}.Region); 
                            t_c = cputime -t3;
                            compareprob_info = [compareprob_info, [-1; -2; -3], [J_diff_min exitflag_diff t_bmi; J_diff_min_g exitflag_diff_g t_g; J_diff_min_c exitflag_diff_c t_c]];
                        end     
                        %______________________________________________________
                        %%Dubugging: having one sample point
                        if options.QPcompare_debug
                            th = exp_sol_part{jreg}.Pn.chebyCenter().x;
                            J_th = th'*J_lower{jreg}.Ai*th + J_lower{jreg}.Bi*th + J_lower{jreg}.Ci;
                            J_upp_th = th'*J_upper.Ai*th + J_upper.Bi*th + J_upper.Ci;
                        end
                        %______________________________________________________
                        % Final check dominance criteria
                        if exitflag_diff && (J_diff_min >= 0) 
                            % min (J_lower-J_upper) >= 0 --> J_lower > J_upper
                            seq_cut_cond= [seq_cut_cond; subprob.level 3]; 
                            node = [subprob.ivalues; -9; subprob.level; -9; 3];%-9 to distinguesh between values, last num=3: dominance
                            node_info = [node_info, node];
                            dominance_cut = true;
                        end
                        %__________________________________________________
                    end %if not infeasible and not inf 
                    %______________________________________________________
                    %______________________________________________________
                      % store upper bound globally
                 elseif options.apply_red_overlap_func
                       % compare lower bound with all upper bounds
                    if  ~infeas_cut && ~isempty(upper_bound) 
                        %__________________________________________________
                        %%Debugging: having one sample point
                        if options.QPcompare_debug
                            th = exp_sol_part{jreg}.Pn.chebyCenter().x;
                            J_th = th'*J_lower{jreg}.Ai*th + J_lower{jreg}.Bi*th + J_lower{jreg}.Ci;
                            J_upp1 = th'*upper_bound{1}.Ai{1}*th + upper_bound{1}.Bi{1}*th + upper_bound{1}.Ci{1};
                        end
                        %______________________________________________________
                              current_exp_sol = part2reg(exp_sol_part{jreg});  % convert for the sake of reduce_overlap_func
                             [upper_bound_tmp,time_overlap_rem,dominance_cut] = reduce_overlap_sol(upper_bound,current_exp_sol,true);                     
                             time_compQP = time_compQP + time_overlap_rem;
                             seq_time_compQP = [seq_time_compQP, time_overlap_rem];
                             
                             if dominance_cut
                                seq_cut_cond= [seq_cut_cond; subprob.level 3]; 
                                node = [subprob.ivalues; -9; subprob.level; -9; 3];%-9 to distinguesh between values, last num=3: dominance
                                node_info = [node_info, node];
                             end
                     end  %if not infeasible and globally
                     %_____________________________________________________
                     %_____________________________________________________
                     % store upper bound locally: J_upper is a cell of possibly more than one upper bounds
                 else
                         % Case: make copy of upper bounds
                     if  ~infeas_cut && J_upper_copy{end}.Ci ~= inf %%( J_lower{jreg}.Ci - J_upper.Ci == -inf )) %&& ~leaf_node 
                        % compare QP functions, using either bmibnb, gurobi, cplex!   
                        N_upper = length(J_upper_copy);
                        J_diff_min = zeros(N_upper,1); exitflag_diff= zeros(N_upper,1);
                        
                        %__________________________________________________
                        %%Debugging: having one sample point
                        if options.QPcompare_debug
                            th = exp_sol_part{jreg}.Pn.chebyCenter().x;
                            J_th = th'*J_lower{jreg}.Ai*th + J_lower{jreg}.Bi*th + J_lower{jreg}.Ci;
                            %J_upp1 = th'*upper_bound{1}.Ai{1}*th + upper_bound{1}.Bi{1}*th + upper_bound{1}.Ci{1};
                            J_upp1 = [];
                        end
                        %__________________________________________________
                        
                        switch options.solver_QPcompare
                            case 'bmibnb'
                                t1 = cputime;
                                for kk = 1 : N_upper
                                    if J_upper_copy{kk}.Ci == inf
                                        continue
                                    end
                                    [~,J_diff_min(kk),~,exitflag_diff(kk)] = compareQPs_bmibnb(J_lower{jreg}, J_upper_copy{kk}, exp_sol_part{jreg}.Region);  %[sol_diff,.,th_min_diff,.] = ..
                                end
                                t_bmi = cputime -t1;
                                time_compQP = time_compQP + t_bmi;
                            case 'cplex'
                                t1 = cputime;
                                for kk = 1 : N_upper
                                    if J_upper_copy{kk}.Ci == inf
                                        continue
                                    end
                                    [~,J_diff_min(kk),~,exitflag_diff(kk)] = compareQPs_cplex(J_lower{jreg}, J_upper_copy{kk}, exp_sol_part{jreg}.Region);
                                end
                                t_c = cputime -t1;
                                time_compQP = time_compQP + t_c;                        
                            case 'gurobi'
                                t1 = cputime;
                                for kk = 1 : N_upper
                                    if J_upper_copy{kk}.Ci == inf
                                        continue
                                    end
                                    [~,J_diff_min(kk),~,exitflag_diff(kk)] = compareQPs_gurobi(J_lower{jreg}, J_upper_copy{kk}, exp_sol_part{jreg}.Region);%Theta
                                    if options.QPcompare_debug && J_diff_min(kk) >=0
                                        J_upp1(end+1,1) = th'*J_upper_copy{kk}.Ai*th + J_upper_copy{kk}.Bi*th + J_upper_copy{kk}.Ci;
                                    end
                                end
                                t_g = cputime -t1;
                                time_compQP = time_compQP + t_g;
                        end   
                        %______________________________________________________
                        %check to see if there is any difference with min difference greater than zero
                        flag_diff = find(exitflag_diff);
                        pos_min_diff = any(J_diff_min(flag_diff) >= 0);                            
                        %______________________________________________________
                        % Final check dominance criteria
                        if any(flag_diff) && (pos_min_diff) %if exitflag_diff && (J_diff_min >= 0) 
                            % min (J_lower-J_upper) >= 0 --> J_lower > J_upper
                            seq_cut_cond= [seq_cut_cond; subprob.level 3]; 
                            node = [subprob.ivalues; -9; subprob.level; -9; 3];%-9 to distinguesh between values, last num=3: dominance
                            node_info = [node_info, node];
                            dominance_cut = true;
                        end
                        %______________________________________________________
                     end %if not infeasible and not inf locally
                     %_____________________________________________________
                 end % if J_upper iscell
                %__________________________________________________________
              end %if compares
            %______________________________________________________________
            %%%%%%%%%%%%%%%%%%%%%%  Cut BnB: MIQP %%%%%%%%%%%%%%%%%%%%%%%%%
            %______________________________________________________________
            % subproblem is infeasible or dominance criteria hold, cut
            if  infeas_cut || dominance_cut 
                    %upper bound is already updated
                    if infeas_cut
                        disp('infeasible')   % subprob is infeasible
                    else
                        disp('dominance cut') %subprob is dominance cut
                    end
                     if ~iscell(J_upper)
                         J_upper_i = J_upper; % not update upper bound
                         x_upper_i = x_upper;                      
                     end
            %__________________________________________________________________________ 
              % subproblem is integer feasible, cut              
            elseif (int_sol_ind(jreg) || leaf_node || subprob.level == length(main_prob.vartype)) % second and third conditions are equal
                    
                    % subproblem is solved 
                    seq_cut_cond= [seq_cut_cond; subprob.level 2];  
                    node = [subprob.ivalues; -9; subprob.level; -9; 2];%-9 to distinguesh between values, last num=2: int feas
                    node_info = [node_info, node];

                     if ~iscell(J_upper)
                        J_upper_i = J_lower{jreg}; % update (local) upper bound
                        x_upper_i = x_lower{jreg}; % update (local) upper bound
                        exp_int_sol{end+1} = exp_sol_part{jreg};
                     else   
                        J_upper_copy{end+1} = J_lower{jreg}; % store (local) upper bound                        
                        exp_int_sol_copy{end+1} = exp_sol_part{jreg}; 
                        x_upper_copy{end+1} = x_lower{jreg}; 
%                         %________________________
%                         %WARM start : update AS
%                             % approach 1: just consider current AS
%                         %AS_parent =  exp_sol_part{jreg}.AS;  
%                             % approach 2: append current AS to the parent AS
%                             options.add_AS = true;
%                             if ~options.add_AS
%                                 % just consider the current parent active sets
%                                 AS_parent =  exp_sol_part{jreg}.AS;
%                             else
%                                % append all parent's active sets
%                                 AS_mix = [ AS_copy, exp_sol_part{jreg}.AS];
%                                 AS_parent = sort(AS_mix);
%                             end
%                             %________________________
                     end

            %__________________________________________________________________________ 
            % subproblem is feasible, but not integer feasible, branch further             
            else
                    seq_cut_cond= [seq_cut_cond; subprob.level 0];
                    node = [subprob.ivalues; -9; subprob.level; -9; 0];%-9 to distinguesh between values, last num=0: branch further
                    node_info = [node_info, node];
                
                     if ~iscell(J_upper)
                         J_upper_i = J_upper; % not update upper bound
                         x_upper_i = x_upper;
                     %else
                     %    J_upper_i = [J_upper_i, J_upper]; % not update upper bound
                     end
                     %__________________________________________________________________________
                    %________________________
                    %WARM start : update AS
                        % approach 1: just consider current AS
                    %AS_parent =  exp_sol_part{jreg}.AS;  
                        % approach 2: append current AS to the parent AS
                        options.add_AS = true;
                        if ~options.add_AS
                            % just consider the current parent active sets
                            AS_parent =  exp_sol_part{jreg}.AS;
                        else
                           % append all parent's active sets
                            AS_mix = [ AS_copy, exp_sol_part{jreg}.AS];
                            AS_parent = sort(AS_mix);
                        end
                        %________________________
                    % WARM start: update subproblems active set
                    subprob.AS =  exp_sol_part{jreg}.AS;
                    subprob.AS_tot =  exp_sol_part{jreg}.AS;
                    %__________________________________________________________________________
                    % split subproblem to two new subproblems: P0 and P1
                    disp(['Splitting feasible space (' num2str(jreg) '/' num2str(N_regions) ')']);                    
                    % branchvar:position of branching integer variable within the set of integer variables in present subprob  
                    branchvar            = 1; %decision(x, subprob.vartype, options.branchrule);
                    % separate subproblem
                    options.update_form = false; % do not change problem formulation for subproblems
                    if options.update_form
                        [p0,p1,zeroOK,oneOK] = separate(subprob, branchvar,options);  
                    else
                        [p0,p1,zeroOK,oneOK] = separate_in2e(subprob, branchvar,options);
                    end
                    % push two subproblems to T_stack
                    %T_stack_copy = push(T_stack_copy,p1);
                    %T_stack_copy = push(T_stack_copy,p0);
                    
                    switch options.method
                        case 'depth'
                            cost = 1/(subprob.level+1);
                        case 'breadth'
                            cost = subprob.level+1;
                        case 'best'
                            cost = subprob.e; % Best-first. This tends to go breadth-first
                        case 'bestdepth'
                            cost = subprob.e/(subprob.level+1);  %zpi %This privilegiates deep nodes
                    end
                    if options.order == 0
                        if oneOK
                            %[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p1,cost);
                             T_stack_copy = push(T_stack_copy,p1);
                        end
                        if zeroOK
                            %[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p0,cost);
                            T_stack_copy = push(T_stack_copy,p0);
                        end
                    else
                        if zeroOK
                            %push_BnB(p0,cost);
                            [T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p0,cost);
                        end
                        if oneOK
                            %push_BnB(p1,cost);
                            [T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p1,cost);
                        end
                    end  %if options.order ...
                   
                    if(jreg>=2)
                        N_nodes_total = N_nodes_total + 2*(2^length(subprob.vartype-1+1)-1); % Adding the extra nodes from the _two_ branches (hence 2*...) _below_ the current node (which should not be added, hence, ...-1...).
                    end
                    %______________________________________________________                    
            end % cut conditions            
            %_________________________________________________________
            % save some node info for final partition
            info_copy.node_info = node; %[info_copy.node_info, node];
            info_copy.ivalues = subprob.ivalues; %[info_copy.ivalues, subprob.ivalues];
            info_copy.sizeTree = prob_num;
         %__________________________________________________________________________
         %__________________________________________________________________________
               % enumeration case 
          else
               %Don't check cut conditions and explicitly enumerate all nodes
               if (leaf_node || subprob.level == length(main_prob.vartype)) % second and third conditions are equal
                    % subproblem is integer feasible, cut 
                    seq_cut_cond= [seq_cut_cond; subprob.level -2];                     
                    J_upper_i = J_lower{jreg}; % update upper bound
                  %__________________________________________________________________________ 
                    %check feasibility, % check: cut here or continue?
               elseif infeas_cut && options.enum_check_inf
                    seq_cut_cond= [seq_cut_cond; subprob.level -1];
                %__________________________________________________________________________ 
                    % subproblem is feasible, but not integer feasible, branch further            
                else
                    seq_cut_cond= [seq_cut_cond; subprob.level -4];
                    J_upper_i = J_upper; % not update upper bound
                    %__________________________________________________________________________
                    % WARM start: update subproblems active set
                    subprob.AS =  exp_sol_part{jreg}.AS;
                    %__________________________________________________________________________
                    % split subproblem to two new subproblems: P0 and P1
                    disp(['Splitting feasible space (' num2str(jreg) '/' num2str(N_regions) ')']);
                    branchvar = 1; %          = decision(x, subprob.vartype, options.branchrule,subprob.branching_priority);
                    [p0,p1,zeroOK,oneOK] = separate(subprob, branchvar,options); 
                    
                    switch options.method
                        case 'depth'
                            cost = 1/(subprob.level+1);
                        case 'breadth'
                            cost = subprob.level+1;
                        case 'best'
                            cost = zpi; % Best-first. This tends to go breadth-first
                        case 'bestdepth'
                            cost = zpi/(subprob.level+1); % This privilegiates deep nodes
                    end
                    if options.order == 0
                        if oneOK
                            %[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p1,cost);
                             T_stack_copy = push(T_stack_copy,p1);
                        end
                        if zeroOK
                            %[T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p0,cost);
                            T_stack_copy = push(T_stack_copy,p0);
                        end
                    else
                        if zeroOK
                            %push_BnB(p0,cost);
                            [T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p0,cost);
                        end
                        if oneOK
                            %push_BnB(p1,cost);
                            [T_stack_copy,T_stack_cost] = push_BnB(T_stack_copy,T_stack_cost,p1,cost);
                        end
                    end  %if options.order ...
                   
                    if(jreg>=2)
                        N_nodes_total = N_nodes_total + 2*(2^length(subprob.vartype-1+1)-1); % Adding the extra nodes from the _two_ branches (hence 2*...) _below_ the current node (which should not be added, hence, ...-1...).
                    end
               end % if leafnode
                   %__________________________________________________________________________
           end % if check cut cond.
           %_______________________________________________________________
           N_nodes_processed = N_nodes_processed + 1;
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%    Update Stacks    %%%%%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %__________________________________________________________________________
            % final check that J_upper_i isn't empty
            if isempty(J_upper_i)
                %disp('J_upper_i is empty!')
                J_upper_i = J_upper;
                x_upper_i = x_upper;
            end
            %______________________________________________________________
            % make sure J_upper_i is cell when there are more than one upper bound
            if iscell(J_upper) && ~iscell(J_upper_i)
                J_upper_i = {J_upper_i};
            end
            %__________________________________________________________________________
            % update S_stack by pushing new tree to 
            new_S = [];
            new_S.Theta = exp_sol_part{jreg}.Region;
            new_S.iter = iter + exp_sol_part{jreg}.iter; 
            new_S.Tree = T_stack_copy;         % update tree
            if ~iscell(J_upper)
                new_S.J_upper = J_upper_i;     % update J_upper!
                new_S.x_upper = x_upper_i;
                new_S.exp_int_sol = exp_int_sol;
            else
                new_S.J_upper = J_upper_copy;  % update J_upper!
                new_S.x_upper = x_upper_copy;
                new_S.exp_int_sol = exp_int_sol_copy;
            end
            new_S.info = info_copy;
            prob_info = [];
            prob_info.level = subprob.level;
            prob_info.ivalue = subprob.ivalues;
            prob_info.prob_num = prob_num;
            prob_info.jreg = jreg;
            prob_info.allreg = N_regions;
            new_S.prob_info = prob_info;
            %store parent and active set info
            % precise but time consuming
            %[WS.parent,WS.AS,WS.AS_logic] = parent_AS_stack(WS.parent,exp_sol_part{jreg}.parent,WS.AS,exp_sol_part{jreg}.AS,WS.AS_logic,exp_sol_part{jreg}.AS_log);
            if ~iscell(WS.parent)
                WS.parent = {WS.parent};
            end
            WS.parent = [WS.parent,{exp_sol_part{jreg}.parent}];
            if ~iscell(WS.AS)
                WS.AS = {WS.AS};
            end
            WS.AS = [WS.AS,{exp_sol_part{jreg}.AS}];
            if ~iscell(WS.AS_logic)
                WS.AS_logic = {WS.AS_logic};
            end
            WS.AS_logic = [WS.AS_logic,{exp_sol_part{jreg}.AS_log}];
            new_S.WS_info = WS;
            %new_S.exp_int_sol = exp_int_sol;
            
            
            new_S.AS_par = AS_parent;  % WARM start
            
            S_stack = push(S_stack,new_S); 
    %__________________________________________________________________________
          end % MILP or MIQP
        end %for jreg
        %_______________________________________________________________
        N_nodes_processed = N_nodes_processed + 1;   
        %__________________________________________________________________
        % update upper bound
              if int_sol_exist % update explicit solution
                  upper_bound =[ upper_bound, int_regions];
                  %________________________________________________________
                  %reduce overlap at each QP iteration
                  if options.reduceOverlap_eachstep
                      %options.compare_upp_low = false;  %compare upper bounds, nou with lower bound
                      if options.save_both_Red_notRed_sol
                          %save both reducedOverlap and original integer solutions
                          upper_bound_redOver = [upper_bound_redOver, int_regions];
                          [upper_bound_redOver,t_redOv] = reduce_overlap_sol(upper_bound_redOver,[],false);
                          time_redOver_upp = time_redOver_upp + t_redOv;
                      else
                          [upper_bound,t_redOv] = reduce_overlap_sol(upper_bound,[],false);
                          time_redOver_upp = time_redOver_upp + t_redOv;
                      end
                  end                                    
                  %________________________________________________________
              end
      end %if feasible
      %seq_J_lower{end+1} = J_lower; % store lower bound for analysis
      %____________________________________________________________________
    else 
        % Tree is empty, update final stack
        N_final_part = length(final_partition);
        final_partition{N_final_part+1}.Region = Theta;
        final_partition{N_final_part+1}.iter = iter;        
        if ~iscell(J_upper) 
            final_partition{N_final_part+1}.J = J_upper;
            final_partition{N_final_part+1}.x = x_upper;
        else
            J_uppr = cell(0);    x_uppr = cell(0); % to skip initial inf 
            for kk = 1: length(J_upper) % erase iinitial infinite upper bound
                if J_upper{kk}.Ci ~= inf
                    J_uppr{end+1} = J_upper{kk};
                    x_uppr{end+1} = x_upper{kk};
                end
            end
            final_partition{N_final_part+1}.J = J_uppr;
            final_partition{N_final_part+1}.x = x_uppr;
         end
        final_partition{N_final_part+1}.node_info = info;
        %final_partition{N_final_part+1}.node_info.exp_sol = exp_int_sol;
        %final_partition{N_final_part+1}.prob_info = prob_info;
        final_partition{N_final_part+1}.sizeTree = info.sizeTree;
        final_partition{N_final_part+1}.node_seq = info.ivalues;
        %final_partition{N_final_part+1}.WS_info = WS;
        %__________________________________________________________________   
        % store sequence of integer values       
        %final_partition{N_final_part+1}.int_values =  seq_ivalues_node;
        seq_ivalues_node = []; %clean to take next node info
        
        %__________________________________________________________________
        F_stack = push(F_stack,final_partition{N_final_part+1});
        optP = prob_num;
        seq_optP = [seq_optP; optP];
        
%         figure;
%         for ii= 1: length(F_stack)
%             F_stack{ii}.Region.plot; 
%             hold on;
%         end
        
    end % isempty(T_stack)
%__________________________________________________________________________
  end %end while
  % Final sol
  explicit_sol = upper_bound; 
  time_bnb = cputime - tstart;
  %________________________________________________________
  % Final reduce overlap in explicit solution
  if options.final_red_overlap
     options.compare_upp_low = false;  % to reduce overlap in final sol
     [explicit_sol_remOverlap,time_final_overlap,~] = reduce_overlap_sol(explicit_sol,[],true);
  end  
%__________________________________________________________________________  
%   % Store some information of BnB 
if ~options.lowmemory
  bnb_info.ivalues = seq_ivalues;  
  bnb_info.levels = seq_levels;  
  bnb_info.Nr_reg_int_ratio = seq_N_reg_int_ratio;
  bnb_info.seq_int_var = seq_int_var;
  bnb_info.seq_cut_cond = seq_cut_cond;
  bnb_info.seq_whole_node_feas = seq_feas_whole_node;
  bnb_info.seq_whole_node_int_feas = seq_int_feas_sol;
  bnb_info.maxprob_num = prob_num;
  bnb_info.maxQPopt_num = optP;
  if options.all_compareQPs
    bnb_info.compareQPs_info = compareprob_info;
    bnb_info.time_compQP = time_compQP;
  end
  if options.reduceOverlap_eachstep && options.save_both_Red_notRed_sol
      bnb_info.explicit_sol_redOver = upper_bound_redOver;
  end
  if options.final_red_overlap
     bnb_info.explicit_sol_redOver = explicit_sol_remOverlap;
     bnb_info.time_final_overlap = time_final_overlap;
  end
  %bnb_info.seq_time_compQP = seq_time_compQP;
  bnb_info.node_info = node_info;
  if options.reduceOverlap_eachstep
    bnb_info.time_upp_redOver = time_redOver_upp;
  end
  bnb_info.time_compQP = time_compQP;
  bnb_info.time_bnb = time_bnb;
  %________________________________________________________________________
   % store less info
else
    bnb_info.int_values = seq_ivalues;
    bnb_info.maxprob_num = prob_num;
    if options.reduceOverlap_eachstep && options.save_both_Red_notRed_sol
       bnb_info.explicit_sol_remOver = upper_bound_redOver;
    end
    if options.final_red_overlap
         bnb_info.explicit_sol_redOver = explicit_sol_remOverlap;
         bnb_info.time_final_overlap = time_final_overlap;
    end
    bnb_info.node_info = node_info;
    if options.reduceOverlap_eachstep
        bnb_info.time_upp_redOver = time_redOver_upp;
    end
    bnb_info.time_compQP = time_compQP;
    bnb_info.time_bnb = time_bnb;
end
%________________________________________________________________________

%    % reduce overlap in explicit (mixed-integer) solution 
%     explicit_sol_remOverlap = cell(0);
%
%      if options.reduceOverlap
%          % Running one last time with redoverlapopts.mark_only = false to get a acc_int_regions without inf.
%         if(~isempty(explicit_sol)) 
%             disp('Final overlap removal...')
%             redoverlapopts.chk_part = 0;
%             redoverlapopts.verbose = 0; % 2 gives lots of information
%             redoverlapopts.mark_only = false;
%             redoverlapopts.complex = 0;
%             fin_overlap_rem_t_start = cputime;
%             explicit_sol_remOverlap = reduceOverlaps_ext_Sh_cert2(explicit_sol,redoverlapopts);
%             final_overlap_rem_time = cputime - fin_overlap_rem_t_start;
%             if(~iscell(explicit_sol_remOverlap))
%                 explicit_sol_remOverlap = {explicit_sol_remOverlap};
%             end
%             disp('Done!')
%         end 
%         bnb_info.final_overlap_rem_time = final_overlap_rem_time;
%         
%      end % reduceOverlap     
% %_________________________________________________________________________ 
 
end %main function
 
% =========================================================================
%_________________________________________________________________________
% Auxiliary functions (Stack)
%_________________________________________________________________________
function  [stack_out,popped,stack_cost] = pop(stack_in,pop_rule,stack_cost)
    switch pop_rule
	case 'lifo'
	  popped = stack_in{1};
	  stack_out = stack_in(2:end);
      if ~isempty(stack_cost)
        stack_cost(1) = [];
      end
	case 'fifo'
	  popped = stack_in{end};
	  stack_out = stack_in(1:end-1);
      if ~isempty(stack_cost)
        stack_cost(end) = [];
      end
    end
end
%_________________________________________________________________________
function stack_out = push(stack_in,to_add)
     stack_out = {to_add stack_in{:}}; 
end

