% Solve relaxed mpLP and certify it
% Find out integer feasible and noninteger solutions
% Warm Start: use parent Active set as initial active set rather than an empty set?

function [all_regions,N_regions, int_regions,N_int_regions, non_int_regions, feasible, int_feas_sol,int_sol, part,cert_info,expl_sols,known_int_seqs] = ...
            LP_cert(subprob,main_prob,P_theta,cert_opts,known_int_seqs,AS_par,options) 

    N_x = length(subprob.f);
    N_theta = size(subprob.W,2);
    N_cons = size(subprob.A,1);%number of all inequality constraints;
    param_vars=(1:N_theta)'; sol = cell(0);

    if all(subprob.ivalues== -1)
        rootnode = true;
    end
    %__________________________________________________________________________
    % Add lowerbound and upperbound to inequality constraints
    %prob = subprob; % for test
    if ~options.boundery_to_ineq
        % changing boundery to inequality has already been done
        [subprob.A,subprob.b,subprob.W,subprob.lb,subprob.ub,bin_constr,constr_type] = add_bin_cons(subprob.A,subprob.b,subprob.W,subprob.lb,subprob.ub,subprob.vartype);
    else
         N_ineq = size(subprob.A,1);
         bin_constr = sort(N_ineq + [subprob.vartype ;(length(vartype)+ vartype)]);
    end
    %________________________________________________________________
    sol_part = cell(0); part = cell(0); cert_info = cell(0);
    %________________________________________________________________
    %%%%%%%%%%%%%%%%%%%%%% mp-LP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.MILP
        if options.MILP_solver
            options.MILP_mpt = false;   % To solve with mpt3
            options.MILP_yalmip = true; % To solve with yalmip
            [sol_mpt, sol_yalm] = parametricSolver(subprob,P_theta,options);
            
            options.MILP_figures = false;
            if options.MILP_figures
                for i = 1:N_x
                    figure;
                    sol_mpt.xopt.fplot('primal', 'position', i);
                    xlabel('t');
                    ylabel(sprintf('x_%d(t)', i));
                end
                figure;
                sol_mpt.xopt.fplot('obj');
                xlabel('t');
                ylabel('J(t)');
            end
            if ~isempty(sol_yalm)
                N_region = size(sol_yalm.Fi,2);
                feasible = true;               
                sol.Fi = sol_yalm.Fi;
                sol.Gi = sol_yalm.Gi;
                sol.Ai = sol_yalm.Ai;
                sol.Bi = sol_yalm.Bi;
                sol.Ci = sol_yalm.Ci;
                sol.activeConstraints = get_AS_from_mp_sol(subprob.A,subprob.b,subprob.W,sol.Fi,sol.Gi);
                sol.AS = sol.activeConstraints;
                %Pnn = sol_mpt.xopt.Domain;
                
                H_stack = []; K_stack = [];
                for i = 1:N_region
                    Pn = double(sol_yalm.Pn(i));
                    H = Pn(:,1:end-1); K = Pn(:,end); 
                    H_stack = [H_stack; H]; K_stack = [K_stack; K];  % for union of all created regions
                    Poly(i) = Polyhedron(H,K);
                    sol.state{i} = 2;  % >2 if it is infeasible
                    sol.iter{i} = 1;   % complexity measure, iter=1: number of nodes
                    sol.parent{i} = false(size(subprob.A,1));   
                    sol.parent{i}(sol.AS{i}) = true; % dummy parent
                end
                Poly.minHRep;
                sol.Pn = Poly;
                sol.Region = sol.Pn;
                
                %%Do some debugging for N_exm = 2
                options.debugging = false;
                if options.debugging
                    %th_samp = [4.5, 0.5, 0.5]'; P_ind = find(sol.Pn.contains(th_samp));
                    th_ch = P_theta.chebyCenter.x;
                    ch_ind = sol.Pn.contains(th_ch);
                    th_samp = sol.Pn(1).chebyCenter.x;
                    P_ind = 1;
                    x_samp = sol.Fi{P_ind}*th_samp + sol.Gi{P_ind};
                    J_samp = subprob.f'*x_samp;
                    % J_samp = sol.Bi{P_ind}*th_samp + sol.Ci{P_ind}; % or
                end
                
        %%_Different alternatives to find the partial infeasible region
            options.try_diff_approach_restreg = false;
            if options.try_diff_approach_restreg
                if length(sol.Pn) == 1
                    Theta_tot = sol.Pn;
                else
                    Theta_tot = Polyhedron(H_stack,K_stack);
                    Theta_tot.minHRep;
                end                
                Theta_rest1 = P_theta\Theta_tot; %P_theta - Theta_tot;
                Theta_rest1.minHRep;
                figure; Theta_rest1.plot;
                
                H_tot = Theta_tot.H(:,1:end-1); K_tot = Theta_tot.H(:,end);
                H_theta = P_theta.H(:,1:end-1); K_theta = P_theta.H(:,end); 
                Theta_rest2 = Polyhedron([H_theta; -H_tot], [K_theta; K_tot]);
                figure; Theta_rest2.plot;
                
                Theta_rest3 = P_theta\Theta_tot;
                figure; Theta_rest3.plot;
                 
                Theta_minus = uminus(sol.Pn);
                figure; Theta_minus.plot;
                % %_________________________________________________
               % see difference of bounding boxes
                options.bounding_box = false;
                if options.bounding_box
                    Box_tot = Theta_tot.outerApprox(); 
                    Box_theta = P_theta.outerApprox();

                    lb_th = Box_theta.Internal.lb; ub_th = Box_theta.Internal.ub;

                    figure; Box_tot.plot; axis([lb_th ub_th, lb_th ub_th, lb_th ub_th]);
                    figure; Box_theta.plot; axis([lb_th ub_th, lb_th ub_th, lb_th ub_th]);
                    Box_rest = Box_theta - Box_tot;
                    figure; Box_rest.plot; axis([lb_th ub_th, lb_th ub_th, lb_th ub_th]);
                    
                    %th_ch_tot = Theta_tot.chebyCenter.x; ch_tot_ind = P_theta.contains(th_ch_tot);                                       
                    figure; Theta_rest.plot; axis([lb_th ub_th, lb_th ub_th, lb_th ub_th]);
                    figure; sol.Pn.plot; axis([lb_th ub_th, lb_th ub_th,lb_th ub_th]); %Theta_tot
                    figure; P_theta.plot; axis([lb_th ub_th, lb_th ub_th, lb_th ub_th]);
                end
            end
            %________________________________________________________
            
                Theta_rest = mldivide(P_theta,sol.Pn); % difference of two polyhedra
                                                  
                ind_part_infs = find(Theta_rest.isFullDim);
               
                if  ind_part_infs  %if ~isempty(Theta_rest) %if N_part_infeas == 1 && Theta_rest.isFullDim                                                             
                    %________________________________
%                    %%check if the leftover region is convex
%                         empty_poly = Polyhedron;
%                         Theta_un = PolyUnion([Theta_rest,empty_poly]);
%                         if ~Theta_un.isConvex()
%                             keyboard;
%                             figure; Theta_rest.plot;                    
%                             % %Theta_rest_convex = Theta_rest.affineHull;
%                             % Theta_rest_convex = Theta_rest.forEach(@affineHull, 'UniformOutput', false);
%                             % %Theta_un_convex = Theta_un.affineHull
%                             % 
%                             % %%outer approximation of the Theta_rest
%                             % Theta_outApprox = Theta_rest.outerApprox;
%                             % figure; Theta_outApprox.plot;                   
%                         end   
                    %________________________________
                    ind_convex_reg = []; %Theta_union = Polyhedron;
                    for ii = 1: length(Theta_rest)
                        % find non-convex region
                        empty_poly = Polyhedron; 
                        Theta_uni = empty_poly; 
                        Theta_uni = PolyUnion([Theta_rest(ii),empty_poly]);
                        %Theta_union(ii) = Theta_uni;
                        ind_convex_reg(ii) = Theta_uni.isConvex();
                    end
                    ind_nonconvex_reg = ~ind_convex_reg;
                    if any(ind_nonconvex_reg)
                        keyboard;
                    end                    
                    %________________________________

                    for kk = 1 : length(ind_part_infs)
                        sol.Fi{N_region+ind_part_infs(kk)} = zeros(N_x,N_theta);
                        sol.Gi{N_region+ind_part_infs(kk)} = nan(N_x,1);
                        sol.Ai{N_region+ind_part_infs(kk)} = zeros(N_theta);
                        sol.Bi{N_region+ind_part_infs(kk)} = zeros(1,N_theta);
                        sol.Ci{N_region+ind_part_infs(kk)} = inf;
                        sol.activeConstraints{N_region+ind_part_infs(kk)} = [];
                        sol.AS{N_region+ind_part_infs(kk)} = [];                   
                        sol.state{N_region+ind_part_infs(kk)} = 3;  % infeasible region
                        sol.iter{N_region+ind_part_infs(kk)} = 1;   % complexity measure, iter=1: number of nodes
                        sol.parent{N_region+ind_part_infs(kk)} = false(size(subprob.A,1));
                        sol.Pn(N_region+ind_part_infs(kk)) = Theta_rest(kk).minHRep; 
                        sol.Region(N_region+ind_part_infs(kk)) = Theta_rest(kk);                         
                    end % for
                    sol.partial_infeas = true;
                end % if ind_partial_infeas
                %___________________________________________________
                    
            else
                 feasible = false;
                 sol.Pn = [];
                 sol.activeConstraints = [];
                 sol.Fi = [];
                 sol.Gi = [];
                 sol.Ai = [];
                 sol.Bi = [];
                 sol.Ci = inf;
                 sol.state = 3;  % >2 if it is infeasible
                 sol.iter = 1;   % complexity measure, iter=1: number of nodes
                 sol.AS = sol.activeConstraints;
                 sol.parent = []; %false(size(subprob.A,1));   
                 %sol.parent(sol.AS) = true; % dummy parent
            end
            %___________________________________________________________
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        else
            %_____________    Certification of LP     ____________
            % use certification algorithm instead
            % use QP certification temporarily
            %%Cold start
            stack_ = [];
            %[part,cert_info] = AScertificationDualV3(subprob,P_theta,[],cert_opts,stack_);
            AS_parent = []; % cold start
            %[part,cert_info] = LPCertification(subprob,P_theta,AS_parent,cert_opts,stack_);
            if options.warm_start 
               AS_parent = AS_par;
               %[part_WS,cert_info_WS] = LPCertification(subprob,P_theta,AS_parent,cert_opts,stack_);
            end
            [part,cert_info] = LPCertification(subprob,P_theta,AS_parent,cert_opts,stack_);
            %__________________________________________________
            %%plot data for debug
            % max_iter = 0; for i= 1: length(part); max_iter = max(max_iter, part{i}.iter); end
            % colors = distinguishable_colors(max_iter+2); inds= [1,2]; values = zeros(size(subprob.W,2)-length(inds),1); %for debug
            % figure; plot_partition(part,'iter',inds,colors,values,1,P_theta,[0,max_iter])
            %___________________________________________________
            %%cert_info.certification_time
            %part_stats_dual = partition_statistics(part,{'iter','FLOPs','GI_FLOPs','GI_sqrts'},false);
            %cert_info.part_stats_dual = part_stats_dual;
             %_________________________________________________________
%              % check for lambda
%              for kk = 1:length(part)
%                  lam{kk} = -subprob.A(part{kk}.AS,:)'\subprob.f; % lam should be positive
%                  lam_p{kk} = part{kk}.lam;
%                  lam_diff{kk} = lam{kk}-lam_p{kk};
%              end
             %%_________________________________________________________
             % disregard low dimentional regions 
             if ~isempty(part) 
                 if options.disregard_lowerdim_reg
                     for jj = 1: length(part)
                         if part{jj}.Region.isFullDim()
                            sol_part{end+1} = part{jj};
                         end
                     end 
                 else
                  % do not disregard lower dimentional regions    
                  sol_part = part;
                 end
                 %_________________________________________________________
                 %change partitions to regions 
                 if ~isempty(sol_part)
                     sol = part2reg(sol_part); %change partitions to regions        
                 end  
             end
             if ~isempty(sol)       
                 feasible = true;
             else
                 feasible = false;
                 sol.Pn = [];
                 sol.activeConstraints = [];
                 sol.Fi = [];
                 sol.Gi = [];
                 sol.Ai = [];
                 sol.Bi = [];
                 sol.Ci = inf;
                 sol.state = 3;  % >2 if it is infeasible
                 sol.iter = 1;   % complexity measure, iter=1: number of nodes
                 sol.AS = sol.activeConstraints;
                 sol.parent = []; %false(size(subprob.A,1));   
                 %sol.parent(sol.AS) = true; % dummy parent
             end 
                        
        end
        
    %%%%%%%%%%%%%%%%%%%%%% mp-QP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %___________________________________________________
        sol = cell(0);
        %_______________Solve QP_____________
        %%Solve mpQP with yalmip
        % [sol,exitflag]=mp_yalmip(mpQP,Pnn);
        %___________________________________________________
        %________________Certify QP_______________
        %%Cold start
        stack_ = [];
        [part,cert_info] = AScertificationDualV3(subprob,P_theta,[],cert_opts,stack_);
        %___________________________________________________
        %%Warm start: use parent Active set as initial active set rather than an empty set
        %[part,cert_info] = AScertificationDualV3(subprob,P_theta,AS_par,cert_opts,stack_); %warm satrt
        %[part,cert_info] = AScertificationDualV3(subprob,P_theta,[],cert_opts,stack_); %cold start
        %___________________________________________________
        %cert_info.certification_time
        part_stats_dual = partition_statistics(part,{'iter','FLOPs','GI_FLOPs','GI_sqrts'},false);
        cert_info.part_stats_dual = part_stats_dual;
         %_________________________________________________________
         % disregard low dimentional regions 
         if ~isempty(part) 
             if options.disregard_lowerdim_reg
                 for jj = 1: length(part)
                     if part{jj}.Region.isFullDim()
                        sol_part{end+1} = part{jj};                        
                     end
                 end 
             else
              % do not disregard lower dimentional regions    
              sol_part = part;
              %th_ch1 = [2.5;2.5;2.5]; (sol_part{1}.RJ*th_ch1)+sol_part{1}.SJ; 
              %subprob.f'*((sol_part{1}.Fx*th_ch1)+sol_part{1}.Gx)
              %[x1,f_val]  = linprog(subprob.f,subprob.A,subprob.b + subprob.W*th_ch1)
             end
             %_________________________________________________________
             %change partitions to regions 
             if ~isempty(sol_part)
                 sol = part2reg(sol_part); %change partitions to regions        
             end  
         end
         if ~isempty(sol)       
             feasible = true;
         else
             sol.Pn = [];
             sol.activeConstraints = [];
             feasible = false;
         end  
    end
    %_____________________________________________________________________
    % TODO: check below: necessary?
%     if(isempty(sol.Pn)) % Can be infeasible, but also because of SEMIdefiniteness, let's try to regularize and resolve.     
%          
%          
%          subprob.H = subprob.H + 1e-5*eye(size(subprob.H,1));
%          disp('Warning: Subproblem perturbed because of bad numerics');
% 
%          %[sol,exitflag]=mp_yalmip(mpQP,Pnn);
%          [part,cert_info] = AScertificationDualV3(subprob,P_theta,[],cert_opts,stack);
%          cert_info.certification_time
%          part_stats_dual = partition_statistics(part,{'iter','FLOPs','GI_FLOPs','GI_sqrts'},false);
%          
%          %part=part;
% 
%         if ~isempty(part)
%              sol = part2reg(part);                          
%              %activeConstraints = get_AS_from_mp_sol(mpQP.A,mpQP.b,mpQP.W,sol.Fi,sol.Gi);            
%         else
%               sol.Pn = [];
%               sol.activeConstraints = [];
%         end
%     end % end if sol is empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %_________________________________________________________________
    % store Active set in logical
   % if ~islogical(sol.AS{1})
       for k = 1: length(sol.Pn)
           AS = false(length(subprob.b),1);
           if ~isempty(sol.AS{k})
               AS(sol.AS{k})=true;
           end
           sol.AS_log{k} = AS; %logial active constraints
       end
    %end
    %_________________________________________________________________
    %find integer region index
    N_regions = length(sol.Pn); 
    leaf_node = ~sum(subprob.ivalues==-1);
    active_const = sol.activeConstraints; %AS
    
    int_sol = false(length(sol.Pn),1);
    if(~isempty(active_const) && ~leaf_node)
        for i = 1: N_regions
            activeconstr_tmp = false(length(subprob.b),1);
            activeconstr_tmp(active_const{i}) = true;
            active_bin_constr = activeconstr_tmp(bin_constr);
            if(sum(active_bin_constr) >= length(subprob.vartype)) % Note that it cannot happen that xi == 0 and xi == 1 simultaneously.
                int_sol(i) = true;
            end
        end
    elseif(leaf_node)
        int_sol = true(N_regions,1);
    end

    N_int_regions = sum(int_sol);
    if(all(int_sol)) % Only integer solutions possible
        int_feas_sol = true;
    else
        int_feas_sol = false;
    end

    all_regions = sol;
    if feasible
        all_regions.Pfinal = PolyUnion(all_regions.Pn);
        %all_regions.activeConstraints = activeConstraints;
    end
    %___________________________________________________________
    % specify integer regions (sol)
    int_regions = all_regions;
    non_int_sol_ind = bin2ind(~int_sol);
    nr_of_non_int_sol_regions = sum(~int_sol);
    for i = 1:nr_of_non_int_sol_regions
        int_regions.Ai{non_int_sol_ind(i)} = zeros(N_theta);
        int_regions.Bi{non_int_sol_ind(i)} = zeros(1,N_theta);
        int_regions.Ci{non_int_sol_ind(i)} = inf;
    end
    
    %___________________________________________________________
    % specify non integer regions (sol)
    non_int_regions = all_regions;
    int_sol_ind = bin2ind(int_sol);
    for i = 1:N_int_regions
        non_int_regions.Ai{int_sol_ind(i)} = zeros(N_theta);
        non_int_regions.Bi{int_sol_ind(i)} = zeros(1,N_theta);
        non_int_regions.Ci{int_sol_ind(i)} = inf;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do some extra calculation % Not-required!
    %_____________________________________________________________________
    if((sum(int_sol)>0) && (nr_of_non_int_sol_regions < options.tighten_feasible_set_limit))
        %[bndA_new,bndb_new] = double(hull(sol.Pn(~int_sol)));
        sol.Pn(~int_sol) = sol.Pn(~int_sol).minHRep();
%         bndA_new = sol.Pn(~int_sol).A;  % SH
%         bndb_new = sol.Pn(~int_sol).b;  % SH
%     else
%         bndA_new = subprob.bndA;
%         bndb_new = subprob.bndb;
    end
    %________________________________________________________________
    expl_sols = [];
    if(options.add_complete_int_seq_sols)
        int_constr_info = zeros(length(bin_constr),2);
        for i = 1:length(bin_constr)
            ind = find(subprob.A(bin_constr(i),:) == -1);

            if(~isempty(ind))
                int_constr_info(i,:) = [ind 0];
            else
                ind = find(subprob.A(bin_constr(i),:) == 1);
                int_constr_info(i,:) = [ind 1];
            end
        end
    %_______________________________________________________________
    
        [int_seqs,int_seqs_regs] = find_relax_int_seqs(subprob,bin_constr,active_const,int_constr_info); 
        if (size(subprob.H,1)>=2) 
        if(~isempty(int_seqs))
            [expl_sols,known_int_seqs,expl_sols_feas] = fixed_bin_seq_expl_sol_comp2(int_seqs,subprob,main_prob.bndA,main_prob.bndb,param_vars,main_prob.vartype,known_int_seqs); 
            if(sum(expl_sols_feas) ~= length(expl_sols_feas))
                expl_sols = [expl_sols int_regions];
            end
        elseif(leaf_node)
            expl_sols = cell(1,1);
            expl_sols{1} = all_regions;
        end
        end
    end
    %_______________________________________________________________

end %function