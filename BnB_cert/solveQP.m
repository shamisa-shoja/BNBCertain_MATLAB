% Solve relaxed mpQP (not certification)
% Find out integer feasible and noninteger solutions

function [all_regions,N_regions, int_regions,N_int_regions, non_int_regions, feasible, int_feas_sol,int_sol,expl_sols,known_int_seqs] = ...
            solveQP(subprob,main_prob,P_theta,known_int_seqs,options) 

    Nx = length(subprob.f);
    N_theta = size(subprob.W,2);
    param_vars=(1:N_theta)';

    if all(subprob.ivalues== -1)
        rootnode = true;
    end
    
    %__________________________________________________________________________
    % Add lowerbound and upperbound to noninequality constraints
    [subprob.A,subprob.b,subprob.W,subprob.lb,subprob.ub,bin_constr,constr_type] = add_bin_cons(subprob.A,subprob.b,subprob.W,subprob.lb,subprob.ub,subprob.vartype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %__________________________________________________________________________
    % solve mpQP certification
%     stack_ = [];
%     [part,cert_info] = AScertificationDualV3(subprob,P_theta,[],cert_opts,stack_);
%     %cert_info.certification_time
%     part_stats_dual = partition_statistics(part,{'iter','FLOPs','GI_FLOPs','GI_sqrts'},false);
%     cert_info.part_stats_dual = part_stats_dual;
%     
%      sol_part = cell(0);
%      %_________________________________________________________
%      % disregard low dimentional regions 
%      if ~isempty(part) 
%          if options.disregard_lowerdim_reg
%              for jj = 1: length(part)
%                  if part{jj}.Region.isFullDim()
%                     sol_part{end+1} = part{jj};
%                  end
%              end 
%          else
%           % do not disregard lower dimentional regions    
%           sol_part = part;
%          end
%      end
%     %sol = part2reg(part); %change partitions to regions
%%___________________________________________
%      %change partitions to regions 
%      if ~isempty(sol_part)
%          sol = part2reg(sol_part);         
%          feasible = true;
%      else
%          sol.Pn = [];
%          sol.activeConstraints = [];
%          feasible = false;
%      end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%__________________________________________________________________________
    %Solve mpQP
 %%----YALMIP-----
 %[sol,exitflag]=mp_yalmip(subprob,P_theta);
 solver = 'yalmip';
 %%-----MPT3-----
 %solver = 'mpt';
 [sol,exitflag] = mpQP_solution(subprob,P_theta,solver);
 
 %_________________________________________________________ 
     %if isempty(sol.Pn)  % polyhedron
     if ( exitflag == -1 || isempty(double(sol.Pn)) )% polytope % isempty(sol)
        feasible = false;
     else
        feasible = true;
     end
%___________________________________________
 if ~isempty(sol)
     %if feasible 
         if ~isfield(sol,'Pfinal')
            sol.Pfinal = PolyUnion(sol.Pn);
         else
            [H_f,K_f] = double(sol.Pfinal);
            sol.Pfinal = Polyhedron('A', H_f, 'b', K_f); % convert Pn to polyhedron
         end
     %end

     for i = 1: length(sol.Pn) % convert polytope to polyhedron
         H = []; K= [];
        try
            [H,K] = double(sol.Pn(i)); % polytope
         catch
             H = sol.Pn(i).A;     % Polyhedron
             K = sol.Pn(i).b;
         end
         sol.Pn(i) = Polyhedron('A', H, 'b', K); % convert Pn to polyhedron
     end
     sol.Region = sol.Pn;
     activeConstraints = get_AS_from_mp_sol(subprob.A,subprob.b,subprob.W,sol.Fi,sol.Gi);
     sol.activeConstraints = activeConstraints;
     sol.AS = activeConstraints;
 else 
     sol.Pn=[];
     sol.Fi=[];
     sol.Gi=[];
     activeConstraints = [];
     sol.activeConstraints = activeConstraints;
     sol.Pfinal = [];
     sol.dynamics = [];
     sol.overlaps= [];
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____________________________________________________________________
    % TODO: check below: necessary?
%     if(isempty(sol.Pn)) % Can be infeasible, but also because of SEMIdefiniteness, let's try to regularize and resolve.     
%          
%          
%          subprob.H = subprob.H + 1e-5*eye(size(subprob.H,1));
%          disp('Warning: Subproblem perturbed because of bad numerics');
% 
%          %[sol,exitflag]=mp_yalmip(subprob,Pnn);
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
%     if(isempty(sol.Pn))
%         feasible = false;
%     else
%         feasible = true;
%     end
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
    if feasible && ~isfield(sol,'Pfinal')
        all_regions.Pfinal = PolyUnion(all_regions.Pn);
        all_regions.activeConstraints = activeConstraints;
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