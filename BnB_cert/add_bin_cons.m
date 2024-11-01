% *Finds implicitly defined equality constraints defined as inequality
% constraints and redefines them as explicit equality constraints
% *Removes all zero, but feasible, constraints

function [A_ext,b_ext,W_ext,lb,ub,bin_constr,constr_type] = add_bin_cons(A,b,W,lb,ub,vartype)

tol = 1e-10;

Nx = size(A,2);
N_ineq = size(A,1);
N_th = size(W,2);

lbgeqninf = lb > -inf;
ubleqinf = ub < inf;
A_lb = -eye(Nx);
A_lb = A_lb(lbgeqninf,:);
A_ub = eye(Nx);
A_ub = A_ub(ubleqinf,:);
A_ext = [A; A_lb; A_ub];
b_ext = [b; -lb(lbgeqninf); ub(ubleqinf)];
W_ext = [W; zeros(size(A_lb,1),N_th); zeros(size(A_ub,1),N_th)];

constr_type = sort(N_ineq + [bin2ind(lbgeqninf(vartype));(length(lbgeqninf(vartype))+bin2ind(ubleqinf(vartype)))]);


 lbgeqninf2 = find(lb > -inf);   %SS
% ubleqinf2 = find(ub < inf);    %SS
% bin_constr = sort(N_ineq + [vartype ;(length(lbgeqninf2)+ vartype)]);  %SS

if all(lbgeqninf)
   bin_constr = sort(N_ineq + [vartype ;(length(lbgeqninf2)+ vartype)]);  %SS
else
   bin_constr = sort(N_ineq + [bin2ind(lbgeqninf(vartype));(length(lbgeqninf(vartype))+bin2ind(ubleqinf(vartype)))]); 
end

lb = -inf*ones(Nx,1);
ub = inf*ones(Nx,1);


%__________________________________________________________________________
% nr_of_implicit_eq_constr = 0;
% implicit_eq_constr = cell(1);

% for i = 1:nineq_ext
%     eq_constr_found = false;
%     for j = i+1:nineq_ext
%         if(norm([A(i,:) b(i)]/norm([A(i,:) b(i)]) + [A(j,:) b(j)]/norm([A(j,:) b(j)])) < tol)
%             if(~eq_constr_found) % First eq constraint found
%                 nr_of_implicit_eq_constr = nr_of_implicit_eq_constr + 1;
%                 implicit_eq_constr{nr_of_implicit_eq_constr,1} = i;
%                 eq_constr_found = true;
%             end
%             
%             implicit_eq_constr{nr_of_implicit_eq_constr} = [implicit_eq_constr{nr_of_implicit_eq_constr};j];
%         end
%     end
% end
% 
% conv_eq_constr = [];
% 
% for i = 1:nr_of_implicit_eq_constr
%     Aeq = [Aeq;A(implicit_eq_constr{i}(1),:)];
%     beq = [beq;b(implicit_eq_constr{i}(1))];
%     conv_eq_constr = [conv_eq_constr ;length(beq)];
% end
% 
% implicit_eq_constr_mat = unique(cell2mat(implicit_eq_constr));
% 
% A(implicit_eq_constr_mat,:) = [];
% b(implicit_eq_constr_mat) = [];
% 
% % Note that both implicit_eq_constr_mat (by unique) and constr_type (by sort at construction) are sorted.
% while(~isempty(implicit_eq_constr_mat))
%     constr_type(constr_type == implicit_eq_constr_mat(1)) = [];
%     ind_tmp = (constr_type >= implicit_eq_constr_mat(1));
%     constr_type(ind_tmp) = constr_type(ind_tmp) - 1;
%     constr_type = constr_type(constr_type ~= implicit_eq_constr_mat(1));
%     implicit_eq_constr_mat = implicit_eq_constr_mat(2:end)-1;
% end
% 
% % Removal of all zero rows that are feasible (if they are infeasible the solver has to detect has).
% if(~isempty(b))
%     all_zero_feas_ineq_constr = bin2ind((sum(A.^2,2) == 0) & (b >= 0))';
% else
%     all_zero_feas_ineq_constr = [];
% end
% 
% if(~isempty(beq))
%     all_zero_feas_eq_constr = bin2ind((sum(Aeq.^2,2) == 0) & (beq == 0))';
% else
%     all_zero_feas_eq_constr = [];
% end
% 
% A(all_zero_feas_ineq_constr,:) = [];
% b(all_zero_feas_ineq_constr,:) = [];
% 
% % Note that both implicit_eq_constr_mat (by unique) and constr_type (by sort at construction) are sorted.
% while(~isempty(all_zero_feas_ineq_constr))
%     constr_type = constr_type(constr_type ~= all_zero_feas_ineq_constr(1));
%     ind_tmp = (constr_type >= all_zero_feas_ineq_constr(1));
%     constr_type(ind_tmp) = constr_type(ind_tmp) - 1;
%     all_zero_feas_ineq_constr = all_zero_feas_ineq_constr(2:end)-1;
% end
% 
% int_constr_info = zeros(length(constr_type),2);
% 
% for i = 1:length(constr_type)
%     ind = find(A(constr_type(i),:) == -1);
%     
%     if(~isempty(ind))
%         int_constr_info(i,:) = [ind 0];
%     else
%         ind = find(A(constr_type(i),:) == 1);
%         int_constr_info(i,:) = [ind 1];
%     end
% end
% 
% Aeq(all_zero_feas_eq_constr,:) = [];
% beq(all_zero_feas_eq_constr,:) = [];