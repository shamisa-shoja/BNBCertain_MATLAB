%% Desicion and Seperate
% The problem formulation is not changed in this function
% Just inequality constraint is changed to equality related to fixed binary variable

function [p0, p1, zeroOK, oneOK] = separate_in2e(prob, branchvar,options,AS_parent)
if (length(prob.vartype) >= 1)
    Nx = size(prob.f,1);
    
    this  = prob.vartype(branchvar);
    others= [1:this-1, this+1:Nx];
    %others= [1:this-1, this+1:n]';  %SH
    % extract the values of the box bounds for the binary branching variable
    % this is used, to check, whether there are box bounds that do not allow 
    % to set one variable to a particular value
    
    lbbranch = prob.lb(this);
    ubbranch = prob.ub(this);
    
    if (lbbranch <= 0) & (ubbranch >= 0)
        zeroOK = 1;
    else
        zeroOK = 0;
    end
    
    if (lbbranch <= 1) & (ubbranch >= 1)
        oneOK = 1;
    else
        oneOK = 0;
    end   
    
    if (zeroOK == 0) & (oneOK == 0)
        error('box constraints on the binary variables are infeasible')
    end    
    
    update_form = options.update_form;
    % Generate new H                                                  
    % Partition old H into 4 blocks, some of which are possibly empty 
    % Note that the old H itself is not empty                         
    
    H11 = prob.H(others, others);
    H12 = prob.H(others, this);
    H22 = prob.H(this, this);
    
    if update_form                  % SH
        p0.H = H11;
        p1.H = H11;
    else                            % SH
        p0.H = prob.H;              % SH % not change problem formulation
        p1.H = prob.H;              % SH
    end                             % SH
  
    if isfield(prob,'H_inv')         % SH
        %p0.H_inv = inv(H11);         % SH
        %p1.H_inv = inv(H11);         % SH
        
        if update_form                  % SH
            p0.H_inv = inv(H11);        % SH
            p0.H_inv = inv(H11);        % SH
        else                            % SH
            p0.H_inv = H_inv;           % SH % not change problem formulation
            p1.H_inv = H_inv;           % SH
        end                             % SH
    
    end
    
    
    % Generate new f                                                  
    % Partition old f into 2 blocks, some of which are possibly empty 
    % Note that a contribution from the partitioning of H is present  
    
    b1  = prob.f(others);
    b2  = prob.f(this);
    
    p0.f = b1;
    p1.f = b1(:)+H12(:);
    
     % Generate new f_theta 
        
     if ~isempty(prob.f_theta)              %SH  
        f_theta = prob.f_theta(others,:);   %SH 
     else
        f_theta = [];                       %SH
     end 
       
     p0.f_theta = f_theta;                    %SH
     p1.f_theta = f_theta;                    %SH
   
     %TODO: theta'*f_theta(this,:)' need to be added to objectove func. or cost of p1 (p1.e)
    
    % Generate new A                                                  
    
    if ~isempty(prob.A)       
        A = prob.A(:,others);
    else
        A = [];
    end                      
    
    p0.A = A;
    p1.A = A;
    
    % Generate new b                                                  
    % The only modification is a contribution from the matrix A       
    
    if ~isempty(prob.A)      
        p0.b = prob.b;
        p1.b = prob.b - prob.A(:,this);
    else                       
        p0.b = prob.b;          
        p1.b = prob.b;          
    end                         

   
  % Generate new W 
        
   W = prob.W; %(:,others); %No..
   
   p0.W = W;                               %SH
   p1.W = W;                               %SH
 
    
    % Generate new Aeq                                                  
    
    if size(prob.Aeq,2) > 0      
        Aeq    = prob.Aeq(:,others);
        p0.Aeq = Aeq;
        p1.Aeq = Aeq;
    else
        p0.Aeq	= prob.Aeq;
        p1.Aeq =  prob.Aeq;
    end	 
    
    % Generate new beq                                                 
    % The only modification is a contribution from the matrix A       
    
    if ~isempty(prob.beq)
        p0.beq = prob.beq;
        p1.beq = prob.beq - prob.Aeq(:,this);
    else
        p0.beq = prob.beq;
        p1.beq = prob.beq;	 
    end
   
    
    % Generate new lb,ub,x0                                         
    
    if ~isempty(prob.lb),
        lb = prob.lb(others);
    else
        lb = [];
    end
    p0.lb = lb;
    p1.lb = lb;
    
    if ~isempty(prob.ub),
        ub = prob.ub(others);
    else
        ub = [];
    end
    p0.ub = ub;
    p1.ub = ub;
  
    if isfield(prob,'bndA')
        p0.bndA = prob.bndA;
        p1.bndA = prob.bndA;
    end
    if isfield(prob,'bndA')
        p0.bndb = prob.bndb;
        p1.bndb = prob.bndb;
    end
  
    
    % Generate new vartype
    
    %EX:              1 2 3 4 5 6 7 8 9 
    %       old_vartype=[    3 4 5     8]'
    %       branchvar=5
    %       newvartype= [    3 4     7]'
    vartype = [prob.vartype(1:branchvar-1); ...
            prob.vartype(branchvar+1:length(prob.vartype))-1];
    
    % Collect the terms for the new subproblems                     
    
    p0.vartype = vartype;
    p1.vartype = vartype;
       
    % Find the absolute index of the branching variable
    ifree   = find(prob.ivalues==-1); % Collect free integer variables
    ibranch = ifree(branchvar);       % Pick up the branch variable 
    
    aux         = prob.ivalues;
    aux(ibranch)= 0;
    p0.ivalues  = aux;   
    aux(ibranch)= 1;
    p1.ivalues  = aux;
    
    p0.level = prob.level+1;
    p1.level = prob.level+1;
    
      % Generate new x0
  if isstruct(prob.x0)                     
      x0.F=prob.x0.F(others,:);                   
      x0.G=prob.x0.G(others,:);            
  elseif ~isempty(prob.x0)
        x0 = prob.x0(others);
    else
        x0 = [];
   end
    
    p0.x0 = x0;
    p1.x0 = x0;
    
    % Generate new e                                                  
    % The only modification is a contribution from H22, b2            
    
    p0.e = prob.e;
    p1.e = prob.e+ .5*H22 + b2; 
    %p1.e = prob.e+ .5*H22 + b2+ theta'*f_theta(this,:)'; %TODO add to final answers in acc_int    %SH

    
    p0.xb(this)=0;                  
    p1.xb(this)=1;                    
    
    if isfield(prob,'branching_priority')
         if(~isempty(prob.branching_priority))
            p0.branching_priority = prob.branching_priority(2:end)-(prob.branching_priority(2:end)>prob.branching_priority(1));
            p1.branching_priority = prob.branching_priority(2:end)-(prob.branching_priority(2:end)>prob.branching_priority(1));
        else
            p0.branching_priority = [];
            p1.branching_priority = [];
         end
    end
    
    if isfield(prob,'explicit_eval')
        p0.explicit_eval = prob.explicit_eval;
        p1.explicit_eval = prob.explicit_eval;
    end
    
   if isfield(prob,'integer_times')
       integer_times = [prob.integer_times(1:branchvar-1); prob.integer_times(branchvar+1:length(prob.integer_times))];     

       p0.integer_times = integer_times; 
       p1.integer_times = integer_times; 
   end
      
   if isfield(prob,'box_ind')
       p0.box_ind = prob.box_ind;
       p1.box_ind = prob.box_ind;
   end
      
   if isfield(prob,'upper_ind')
       p0.upper_ind = prob.upper_ind;
       p1.upper_ind = prob.upper_ind;
   end
   
   if isfield(prob,'lower_ind')
       p0.lower_ind = prob.lower_ind;
       p1.lower_ind = prob.lower_ind;
   end
   
   
   p0.AS0=prob.AS0;
   p1.AS0=prob.AS0;
   
   p0.x0_case= prob.x0_case;
   p0.normal = prob.normal;
%         p0.Fx= prob.Fx;
%         p0.Gx= prob.Gx;
 

    p1.x0_case= prob.x0_case;
    p1.normal= prob.normal;
%         p1.Fx= prob.Fx;
%         p1.Gx= prob.Gx;

   
    if isfield(prob,'iter')
         p0.iter = prob.iter;
         p1.iter = prob.iter;
    end
        % save parents iter number and add up for each its child node
        if isfield(prob,'iter_stack')
             p0.iter_stack = prob.iter_stack;
             p1.iter_stack = prob.iter_stack;
        end
        
        if isfield(prob,'parent_stack')
             p0.parent_stack = prob.parent_stack;
             p1.parent_stack = prob.parent_stack;
        end
        
        if isfield(prob,'AS_stack')
             p0.AS_stack = prob.AS_parent;
             p1.AS_stack = prob.AS_parent;
        end
        
    if isfield(prob, 'xxmax')
        if ~isempty(prob.xxmax),
            xxmax = prob.xxmax(others);
        else
            xxmax = [];
        end
        p0.xxmax = xxmax;
        p1.xxmax = xxmax;
    end

   % Keep parameter space info
   if isfield(prob,'Theta')
       p0.Theta=prob.Theta;
       p1.Theta=prob.Theta;
   end
   
   if isfield(prob, 'lb_th')
       p0.lb_th= prob.lb_th;
       p0.ub_th= prob.ub_th;

       p1.lb_th= prob.lb_th;
       p1.ub_th= prob.ub_th;
       
%       p1.A_th= [];  % Polyhedron of parameter space
%     	p1.b_th= []; 

   end
   
   if isfield(prob, 'Q_th')
       p0.Q_th=prob.Q_th; p1.Q_th=prob.Q_th;
       p0.R_th=prob.R_th; p1.R_th=prob.R_th;
       p0.S_th=prob.S_th; p1.S_th=prob.S_th;
   end
 
   % Warm start 
    %Define starting cond. for active set in warm start
        p0.AS_parent = []; 
        p1.AS_parent = [];
        
        %m = length(prob.b); %number of constraints %TODO
        m = length(p0.b); %number of constraints %TODO
        nb = length(prob.vartype); %%number of binary variables
        
        if ~isempty(AS_parent)
            options.notkeep_AS_order = true;
            %___________________________
            % Approach 1: order does not remain
            if options.notkeep_AS_order
            % ignore AS for current fixed binary value: ex. AS = [1 3 5 8] --> AS = [1 3 6];
                AS_ineq = AS_parent(AS_parent <= m); %ActiveSet of original inequality constraints
                %AS_b = AS_parent(AS_parent > m);     %AS of left binary constraints
                AS_b_lb = AS_parent(AS_parent < m+nb+1 & AS_parent > m);    %AS of lower bound of left binary constraints
                AS_b_ub = AS_parent(AS_parent > m+nb);     %AS of upper bound of left binary constraints

                %AS_par = [AS_ineq, AS_b(AS_b <= m+branchvar-1), AS_b(AS_b > m+branchvar+1)-2];
        %         AS_par = [AS_ineq, AS_b_lb(1:branchvar-1), AS_b_lb(branchvar+1:length(AS_b_lb))-1,...
        %                   AS_b_ub(1:branchvar-1), AS_b_ub(branchvar+1:length(AS_b_ub))-2];
                 AS_par = [AS_ineq, AS_b_lb(AS_b_lb < m+branchvar), AS_b_lb(AS_b_lb > m+branchvar)-1,...
                          AS_b_ub(AS_b_ub < m+nb+branchvar), AS_b_ub(AS_b_ub > m+nb+branchvar)-2];
            else
                %____________________________________
                % Approach 2: order of AS is kept
                  AS_par = prob.AS;
                  ind1 = find(AS_par <= m);
                  ind2 = find(AS_par < m+nb+1 & AS_par > m);              

                  for j = 1:length(ind2)
                      if AS_par(ind2(j)) == m+branchvar
                          AS_par(ind2(j)) = 0;
                      else
                          AS_par(ind2(j)) = AS_par(ind2(j)) -1;
                      end
                  end
                  %AS_par(ind2) = AS_par(ind2) -1;

                  ind3 = find(AS_par > m+nb);
                  for j = 1:length(ind3)
                      if AS_par(ind3(j)) == m+nb+branchvar
                         AS_par(ind3(j)) = 0;
                      else 
                          AS_par(ind3(j)) = AS_par(ind3(j)) - 2;
                      end
                  end
                  %AS_par(ind3) = AS_par(ind3) -2;

                  ind = find(AS_par);
                  AS_par = AS_par(ind);
        %___________________________
        %         % Approach 3: order of AS is kept
        %        AS_par = [];
        %        for j = 1 : length(AS_parent)
        %            %if (AS_parent(j) > m) && (AS_parent(j) < m+nb+1)
        %            switch AS_parent(j) 
        %                case (AS_parent(j) < m) 
        %                    AS_par = [AS_par, AS_parent(j)];
        %                case (AS_parent(j) > m) && (AS_parent(j) < m+nb+1)
        %                    AS_par = [AS_par, AS_parent(j)-1];
        %                case (AS_parent(j) > m+nb+1)
        %                    AS_par = [AS_par, AS_parent(j)-2];
        %            end
        %        end
            end
            %______________________________   
            % Collect the terms of AS for the new subproblems
           p0.AS_parent = AS_par; %AS(ind); 
           p1.AS_parent = AS_par; %AS(ind);
        end
           
             % update the entire Active set if there is: Warm start
%         if isfield(prob, 'AS_tot') 
%             m = length(prob.b); %number of constraints
%             integers = sort(prob.vartype(prob.ivalues ~= -1)); %fixed integer variables
%             %AS_tot = unique([prob.AS_tot, m+integers]);
%             %p0.AS_tot = AS_tot; %[AS_tot, m+branchvar]
%             %p1.AS_tot = AS_tot; %[AS_tot, m+branchvar+1]
%             
%             p0.AS_tot = unique([prob.AS_tot, m + branchvar]);   %prob.vartype(branchvar)         
%             %p1.AS_tot = unique([prob.AS_tot, m + prob.vartype(p1.ivalues ~= -1)]);            
%             p1.AS_tot = unique([prob.AS_tot, m + branchvar+1]);
%         end
               
else
%     error('no more integer variables to branch on') %TODO: change to warning
    warning('no more integer variables to branch on') %TODO: change to warning
end

end %function
% ------------------------------------------------------------------------------
%
% decision: when a problem has to be separated, this function decides   
%           which will be the next branching variable                   
%           input: x       = present value of the solution of the qp            
%                  vartype = indices of the free integer variables relative to x
%                  br      = parameter denoting the branching rule that has to  
%                            be adopted     
%          output: branchvar = next branching variable position within vartype

% function branchvar = decision(x,vartype,br,branching_priority)
% switch br
%     case 'first'
%         % first free variable is chosen as branching variable        
%         branchvar = 1;
%     case 'last'
%         % DA: last free variable is chosen as branching variable
%         branchvar = size(vartype,1);
%     case 'max'
%         % integer free variable with maximal frac part is 
%         % chosen as branching variable          
%         xi          = x(vartype);
%         [aux1,aux2] = max(abs(xi-round(xi)));
%         branchvar   = aux2(1); % pick the first of the variables with max value
%     case 'min'
%         % integer free variable with minimal frac part is 
%         % chosen as branching variable          
%         xi          = x(vartype);
%         [aux1,aux2] = min(abs(xi-round(xi)));
%         branchvar   = aux2(1); % pick the first of the variables with max value
%     case 'predefined'
%         branchvar = branching_priority(1);
%     otherwise
%         % decision not implemented                                   
%         warning('decision not implemented: switch to first free');
%         branchvar = 1;
% end
% 
% end %function

