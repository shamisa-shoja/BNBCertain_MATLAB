
function [sol_diff,J_diff_min,th_min_diff,exitflag] = compareQPs_bmibnb(J_low, J_upp, Theta)

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
    
    th = sdpvar(N_theta,1);
    
    cons = Theta.A * th <= Theta.b;
    if isfield(Theta,'lb') && ~isempty(Theta.lb)
        cons2 = Theta.lb <= th <= Theta.ub;
        cons = [cons, cons2];
    end    
    
    obj = th'* A_diff* th + B_diff * th + C_diff;
    
    %ops = sdpsettings('solver','moment','moment.order',3)
    ops = sdpsettings('solver','bmibnb');
    
    try
        sol_diff = optimize(cons,obj,ops);
        J_diff_min = value(obj);
        th_min_diff = value(th);
        %J_diff_min = th_min_diff'*A_diff* th_min_diff + B_diff*th_min_diff + C_diff;           
        exitflag = 1;        
    catch
        sol_diff = [];
        J_diff_min = nan;
        th_min_diff =[];
        exitflag = 0;
    end
   