
function [sol_diff,J_diff_min,th_min_diff,exitflag,Theta_1,Theta_2,flag_reg] = spatialSeparate(J_low, J_upp, Theta, options)
    
    N_theta = size(J_low.Bi,2);
    
    if isfield(J_low, 'Ai') 
        A_low = J_low.Ai;
    else %if ~isfield(J_low, 'Ai') || isempty(J_low.Ai)
        A_low = zeros(N_theta);
    end
    B_low =J_low.Bi; 
    C_low =J_low.Ci; 
    
    
    if isfield(J_upp, 'Ai') 
        A_upp = J_upp.Ai;
    else
        A_upp = zeros(N_theta);
    end
    B_upp =J_upp.Bi; 
    C_upp = J_upp.Ci; 
    
    A_diff = A_low - A_upp;
    B_diff = B_low - B_upp;
    C_diff = C_low - C_upp;

    
    if N_theta ~= size(B_upp,2)
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
        
    % initialize outputs
    sol_diff = [];
    J_diff_min = nan;
    th_min_diff =[];
    exitflag = 0;
    epsi = 1e-4; % epsilon to be consider as zero
    
    % initialize output regions
    Theta1 = Polyhedron; Theta2 = Polyhedron; 
    Theta_1 = Polyhedron; % region to be pruined
    Theta_2 = Theta;      % region to be kept
    flag_reg = 'keep'; % By default keep whole region
           
%__________________________________________________________            
if options.MILP
    %%%%%%%%%%% MILP problem %%%%%%%%%%%%%
    dom_cut = false;     % for Theta_1 to be pruined
    Theta1_empty = false;
    not_dom_cut = false; % for Theta_2 to be kept
    Theta2_empty = false;
      
    % J is affine % TODO: add obvoisly better and worse also
    % Equal_worse_J = (norm(B_diff) <= epsi) && (C_diff >= -epsi)
    % Better_J = (norm(B_diff) <= epsi) && (C_diff < -epsi)
    if (norm(B_diff) <= epsi) && (abs(C_diff) <= epsi)
        % if deltaJ ==0 --> cut it as previous is better
        dom_cut = true;
        flag_reg = 'prune';
        % completely prune region
        Theta_1 = Theta;
        Theta_2 = Polyhedron;
        
    else        
        % J_diff = B_diff*theta+C_diff >= 0
        Theta1 = Polyhedron([H;-B_diff],[K;C_diff]);

        % J_diff = B_diff*theta+C_diff <= 0
        Theta2 = Polyhedron([H;B_diff],[K;-C_diff]);
        %Theta2 = Theta - Theta1;

          if Theta1.isEmptySet()
              Theta1_empty = true;                      
          end
          if ~Theta1_empty && Theta1.isFullDim
              dom_cut = true;
          end

          if Theta2.isEmptySet()
              Theta2_empty = true;
          end
          if ~Theta2_empty && Theta2.isFullDim
              not_dom_cut = true;
          end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if dom_cut && not_dom_cut
          % both finer regions are full dim polyhedra,
          % partition region to Theta_1 & Theta_2
          flag_reg = 'partition';
          % devide region to two
          Theta_1 = Theta1;
          Theta_2 = Theta2;
          
      elseif dom_cut && Theta2_empty 
          flag_reg = 'prune';
          % completely prune region
          Theta_1 = Theta;
          Theta_2 = Polyhedron;
          
      %elseif not_dom_cut && Theta1_empty % no need
      %    flag_reg = 'keep';
          % No dominance cut: keep whole region
       %   Theta_1 = Polyhedron; 
       %   Theta_2 = Theta;
      end
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%       if strcmp(flag_reg,'partition')
%           partial_cut = true;
%       elseif strcmp(flag_reg,'prune')
% 
%           dominance_cut = true;
%       else
%           not_dominance_cut = true;
%       end
%     %____________________________________________
else
    %%%%%%%%%%% MIQP problem %%%%%%%%%%%%%
    % J is quadratic
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
        end
    catch        
	    sol_diff = [];
        J_diff_min = nan;
        th_min_diff =[];
        exitflag = 0;
    end
end
end % function
