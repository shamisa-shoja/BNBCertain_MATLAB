
function [sol_diff,J_diff_min,th_min_diff,exitflag] = compareQPs_gurobi(J_low, J_upp, Theta)

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
    if ~issymmetric(A_diff) % necessary?
        A_diff = (A_diff + A_diff') /2;
    end
    
    H = Theta.A;
    K = Theta.b;
    if isfield(Theta,'lb') && ~isempty(Theta.lb)
        lb = Theta.lb; %+(1E-6); 
        ub = Theta.ub; %-(1E-6);
    else
        lb = -inf(N_theta,1); 
        ub = inf(N_theta,1); 
    end
    
    absolute_mip_tolerance = 0; %1e-6  %In Gurobi: default: 1e-10
    relative_mip_tolerance = 0; %1e-6  %In Gurobi: default: 1e-4
    integer_tolerance = 1e-09; %0  (In Gurobi: minimum is 1e-09 not zero)
    
    params = struct();
    params.IntFeasTol = integer_tolerance;     % Integer feasibility tolerance
    params.MIPGap = relative_mip_tolerance;    % Relative MIP optimality gap
    params.MIPGapAbs = absolute_mip_tolerance; % Absolute MIP optimality gap
    params.OutputFlag = 0;
	params.NonConvex = 2;
	params.Method= 0;
	params.TimeLimit= 2;
	params.FeasibilityTol= 0.5*1e-6;
    
    model.Q = sparse(A_diff);
	model.obj = B_diff;
	model.modelsense = 'min';

	model.A   = sparse(H);
	model.rhs = K;%-(1E-6)
	model.sense = repmat('<',size(H,1),1); 
    
	model.lb= lb; %+(1E-6)
	model.ub= ub; %-(1E-6)
%	N_quad = size(A_diff,1)/N_theta; % 1 in this case!
% 	for n = 1:N_quad
% 	  model.quadcon(n).Qc=sparse(A_diff(((n-1)*N_theta+1):n*N_theta,:));
% 	  model.quadcon(n).q=B_diff(n,:);
% 	  model.quadcon(n).rhs = -C_diff(n); %-(1E-6);
%     end	
    
    try
        result = gurobi(model,params);
        empty=strcmp(result.status,'INFEASIBLE')|strcmp(result.status,'INF_OR_UNBD');

        if(~empty && isfield(result,'x')) % Sometimes return 'NUMERIC' %TODO look into this      
            sol_diff = result;
            th_min_diff = result.x;
            J_diff_min = result.objval + C_diff;
            %J_diff_min = th_min_diff'*A_diff* th_min_diff + B_diff*th_min_diff + C_diff;           
            exitflag = 1;
        else
            sol_diff = [];
            J_diff_min = nan;
            th_min_diff =[];
            exitflag = 0;
        end
    catch        
	    sol_diff = [];
        J_diff_min = nan;
        th_min_diff =[];
        exitflag = 0;
    end

end % function
