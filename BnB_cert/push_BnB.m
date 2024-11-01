
function  [stack_out,stack_cost] = push_BnB(stack_in,stack_cost,to_add,cost)

  % Determine position in STACK where problem is inserted, according to a best first strategy
    ii = find(stack_cost>=cost);  % EX: STACKCOST=[100 80 33 22 ^ 5 3 2], cost=10
    if isempty(ii)
        i=1;
    else
        i=ii(end)+1;
    end
    
    if length(stack_in) > 1
        stack_out = { stack_in{1:i}, to_add, stack_in{i+1:end}};
    else
        stack_out = {to_add stack_in{:}};
    end
    
    stack_size = length(stack_in);
    
    for j = stack_size:-1:i
        %stack_out(j+1)= stack_in(j);
        %stack_out{j+1}= stack_in{j};
        if isempty(stack_cost)
             stack_cost(j)=0;
        end
         stack_cost(j+1)= stack_cost(j);    
    end
%     
%     %if isempty(stack)
%     %    stack = {};
%     %end
%     
%     %stack(i) = to_add; 
%     stack{i}  = to_add;
      stack_cost(i) = cost;
%     %stack_size    = stack_size+1;
end




%___________________________________________________________________________
% function  [stack_out,stack_cost] = push_BnB(stack_in,stack_cost,to_add,cost)
%     
%     stack_out = cell(0);
%   % Determine position in STACK where problem is inserted, according to a best first strategy
%     ii = find(stack_cost>=cost);  % EX: STACKCOST=[100 80 33 22 ^ 5 3 2], cost=10
%     if isempty(ii)
%         i=1;
%     else
%         i=ii(end)+1;
%     end
%      
%     stack_size = length(stack_in);
%     
%     for j = stack_size:-1:i
%         %stack_out(j+1)= stack_in(j);
%         stack_out{j+1}= stack_in{j};
%         if isempty(stack_cost)
%              stack_cost(j)=0;
%         end
%          stack_cost(j+1)= stack_cost(j);    
%     end
% 
%     %stack_out(i) = to_add; 
%     stack_out{i}  = to_add;
%     stack_cost(i) = cost;
%     %stack_size    = stack_size+1;
% end

% function stack_out = push(stack_in,to_add)
%   stack_out = {to_add stack_in{:}};
% end
