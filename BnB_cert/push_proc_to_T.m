
function  [stack_out] = push_proc_to_T(stack_in,to_add,opts)    
  
  epsi = 1e-4; % tolerance for zero
  stack_size = length(stack_in);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if stack_size == 0
      % stack is empty--> push it
      stack_out = {to_add};
  else
      % stack is nonempty--> find the right place in the prioritery order list
      switch opts.method
            case 'depth'
                % parameter-independent node search strategy
                stack_out = {to_add stack_in{:}}; 
            case 'breadth'
                % parameter-independent node search strategy
                for k = 1: stack_size
                    cost_in(k) = stack_in{k}.cost;                    
                end %for
                
                jj = find(to_add.cost <= cost_in);
                if isempty(jj)
                   kk = stack_size;
                else
                   kk = jj(1)-1; 
                end

                 if kk >= stack_size  
                    stack_out = { stack_in{:}, to_add}; 
                 else
                    stack_out = { stack_in{1:kk}, to_add, stack_in{kk+1:end}};
                 end               

            case {'best','bestdepth'}
                % parameter-dependent node search strategy
              flag_larger = false(stack_size,1);
              flag_smaller = false(stack_size,1);
              flag_devide= false(stack_size,1);

              if to_add.cost.Ci == inf 
                 % TODO: check: put in the end?
                 % J_add is larger (infeasible prob) --> no need to compare, put in in the end of the list
                stack_out = {stack_in{:} to_add};
              
              else
                  % Determine position in S_stack where problem is inserted, according to a best first strategy      
                   cheby1= cell(0); cheby2= cell(0);      
                  for k = 1: stack_size
                      if stack_in{k}.cost.Ci == inf
                          % J_add is smaller, no need to compare
                          flag_smaller(k) = true;
%                           if to_add.cost.Ci == inf 
%                              %put in before the inf one
%                              % J_add is larger (infeasible prob) --> no need to compare, put in in the end of the list
%                              %stack_out = {stack_in{:} to_add};
%                              kk = k-1;
%                              stack_out = { stack_in{1:kk}, to_add, stack_in{kk+1:end}};
%                              return
%                           end
                      else

                      %[~,~,~,~,Theta_1(k),Theta_2(k)] = spatialSeparate(J_add, stack_in{k}.J_lower, Theta_add, opts);
                        B_diff{k} = to_add.cost.Bi - stack_in{k}.cost.Bi;
                        C_diff{k} = to_add.cost.Ci - stack_in{k}.cost.Ci;

                        if (norm(B_diff{k}) <= epsi) && (abs(C_diff{k}) <= epsi)
                            % J_add and Ji are equal
                            kk = k-1;
                             if kk >= stack_size  
                                stack_out = { stack_in{:}, to_add}; 
                             else
                                stack_out = { stack_in{1:kk}, to_add, stack_in{kk+1:end}};
                             end
                           return
                        end

                        H = to_add.Theta.A;
                        K = to_add.Theta.b;

                        Theta1 = Polyhedron; Theta2 = Polyhedron; 

                        Theta2(k) = Polyhedron('A', [H;B_diff{k}],'b', [K;-C_diff{k}]); % J_diff <= 0 : B_diff*theta <= -C_diff
                        %Theta1(k) = to_add.Theta - Theta2(k);      % J_diff >= 0 : -B_diff*theta <= C_diff
                        Theta1(k) = Polyhedron('A', [H;-B_diff{k}],'b', [K;C_diff{k}]); % J_diff >= 0

                        flag_larger(k) = ~Theta1(k).isEmptySet();
                        flag_smaller(k) = ~Theta2(k).isEmptySet();

                        if Theta1(k).isFullDim() && Theta2(k).isFullDim()
                            % test if they are not just lines but polyhedra
                           cheby1{k} = Theta1(k).chebyCenter(); 
                           cheby2{k} = Theta2(k).chebyCenter();
                           % based on cheby.r, one can also see if it is not full dim 
                           % (r = 0 --> not full dimential)
                           if ~isinf(cheby1{k}.r) && ~isinf(cheby2{k}.r)
                               flag_devide(k) = true;
                           end
                        end  

                        if flag_devide(k)    
                          % partition region to_add to two new regions 
                           disp('spatial splitting in Pushing to List based on best-first Strategy');
                           to_add_p1 = to_add;
                           to_add_p1.Theta = Polyhedron;
                           to_add_p1.Theta = Theta2;
                           j = k-1; 
                           stack_in = { stack_in{1:j}, to_add_p1, stack_in{j+1:end}};

                           to_add.Theta = Theta1; %continue further with the rest of Theta1
                           stack_size = stack_size + 1;
                        end

                     end
                  end
                  % now put the node (with the left region) into the list
                  if ~to_add.Theta.isEmptySet()&& to_add.Theta.isFullDim()
                      ii = find(~flag_larger & flag_smaller);

                      if isempty(ii)
                          i = stack_size;
                      else
                          i = ii(1)-1; %i=ii(end)+1;
                      end

                     if i >= stack_size  
                        stack_out = { stack_in{:}, to_add}; 
                     else
                        stack_out = { stack_in{1:i}, to_add, stack_in{i+1:end}};
                     end
                  end 
               end    %if cost is inf   
      end %switch
  end %if stack is not empty
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      %%%%%%%%%%%%%%%%%%%%%%%%
      opts.illustrate = false;
      if opts.illustrate
            figure; 
            for j=1:stack_size
                subplot(1,stack_size,j)
                stack_in{j}.Theta.plot
            end
            title('Theta os stack input')

            figure; 
            to_add.Theta.plot
            title('Theta of added in stack')

            Theta_diff = Polyhedron;
            figure; title('Theta os stack input and added stack') 
            for j=1:stack_size
                subplot(2,stack_size,j)
                stack_in{j}.Theta.plot; hold on; to_add.Theta.plot; hold off;
                axis([-2.5 2.5, -2.5 2.5])

                subplot(2,stack_size,stack_size+j)
                Theta_diff(j) = stack_in{j}.Theta - to_add.Theta;
                Theta_diff(j).plot;
                axis([-2.5 2.5, -2.5 2.5])
                title('Diff of Thetas')
            end       
        end
      if opts.illustrate
            figure; 
            for j=1:stack_size
                subplot(1,stack_size,j)
                Theta1(j).plot
            end
            title('Theta1: delta_j >= 0')
            
            figure; 
            for j=1:stack_size
                subplot(1,stack_size,j)
                Theta2(j).plot
            end
            title('Theta2: delta_j <= 0')
            
            figure; 
            for i=1 : stack_size
                subplot(1,stack_size,i);
                Theta2(i).plot;
                %Theta_diff(k) = Theta_1(k)-Theta_2(2);
                %Theta_diff(k).plot
            end
            title('Theta2: delta_j <= 0') %title('difference of Theta: Theta1-Theta2')
         %close all;   
      end
    
end % end function


