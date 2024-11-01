
function [sol_diff,J_diff_min,th_min_diff,exitflag] = compareQPs_cplex(J_low, J_upp, Theta)

    A_low = J_low.Ai;
    B_low =J_low.Bi; 
    C_low =J_low.Ci; 

    A_upp = J_upp.Ai;
    B_upp =J_upp.Bi; 
    C_upp =J_upp.Ci; 

    A_diff = A_low - A_upp;
    B_diff = B_low - B_upp;
    C_diff = C_low - C_upp;

    N_theta = size(A_low,1);

    if N_theta ~= size(A_upp,1)
        error ('incompatible objective functions')
    end
    
    if ~issymmetric(A_diff)
        A_diff = (A_diff + A_diff') /2;
    end
    H = Theta.A;
    K = Theta.b;%-(1E-6)
    if isfield(Theta,'lb') && ~isempty(Theta.lb)
        lb = Theta.lb; %+(1E-6); 
        ub = Theta.ub; %-(1E-6);
    else
        lb = -inf(N_theta,1); 
        ub = inf(N_theta,1); 
    end
    
    absolute_mip_tolerance = 0; %1e-6
    relative_mip_tolerance = 0; %1e-6
    integer_tolerance = 1e-09; %0  (minimum of gurobi is 1e-09)
    max_sol_time = 1200; % 1e15 % maximum solution time in sec
    cplex_options = cplexoptimset('optimalitytarget',3,'output.clonelog',0);
    cplex_options.mip.tolerances.absmipgap = absolute_mip_tolerance; % Identifier: 2008
    cplex_options.mip.tolerances.mipgap = relative_mip_tolerance; % Identifier: 2009
    cplex_options.mip.tolerances.integrality = integer_tolerance; % Identifier: 2010
    cplex_options.timelimit = max_sol_time;  % Identifier: 1039
    
    % case indefinite QP
    if min(eig(A_diff)) < 0
        cplex_options.mip.cuts.disjunctive = 3; % Identifier: 2053
        cplex_options.mip.strategy.variableselect = 4; % Identifier: 2028
        cplex_options.emphasis.mip = 1; % Identifier: 2058
        cplex_options.mip.strategy.presolvenode = 1; % Identifier: 2037
    end
    
    try        
        [th_min_diff,J_diff_min,exitflag,sol_diff] = cplexqp(A_diff,B_diff,H,K,[],[],lb,ub,[],cplex_options);
        %J_diff_min = J_diff_min + C_diff; % apply also C_diff %TODO: check? 
        J_diff_min = th_min_diff'*A_diff* th_min_diff + B_diff*th_min_diff + C_diff;
        %exitflag = 1;
    catch
        sol_diff = [];
        th_min_diff =[];
        exitflag = 0;
        J_diff_min = nan;
    end
