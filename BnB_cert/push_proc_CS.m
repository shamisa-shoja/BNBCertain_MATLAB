function  S_stack_out = push_proc_CS(p0,p1,T_stack,Theta, new_S,S_stack,opts)  %[T_stack_out, S_stack_out]   
  
  epsi = opts.epsi; % 1e-4; % tolerance for zero
  stack_size = length(T_stack);
  S_stack_out = cell(0);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if stack_size == 0
      % stack is empty--> push it
      T_stack_out = {p0, p1};
      % update new tuple to put in S
      new_S.Theta = Theta;
      new_S.Tree = cell(0);
      new_S.Tree = T_stack_out;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Push to S  
      S_stack_out = {new_S S_stack{:}};
      return
  else
      % stack is nonempty--> find the right place in the prioritery order list
      switch opts.method
            case 'depth'
                % parameter-independent node search strategy 
                % put in first always (no need to compare)
                  T_stack_out = {p0, p1, T_stack{:}}; 
                % update new tuple to put in S
                  new_S.Theta = Theta;
                  new_S.Tree = cell(0);
                  new_S.Tree = T_stack_out;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Push to S  
                  S_stack_out = {new_S S_stack{:}};
            case 'breadth'
                % parameter-independent node search strategy
                %cost = subprob.level+1;
                %p0.cost = cost; p1.cost = cost;
                for k = 1: stack_size
                    cost_in(k) = T_stack{k}.cost;                    
                end %for
                
                jj = find(p0.cost- cost_in <= epsi);
                if isempty(jj)
                   kk = stack_size;
                else
                   kk = jj(1)-1; 
                end

                 if kk >= stack_size  
                    T_stack_out = { T_stack{:}, p0, p1}; 
                 else
                    T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                 end  
                 % update new tuple to put in S
                  new_S.Theta = Theta;
                  new_S.Tree = cell(0);
                  new_S.Tree = T_stack_out;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Push to S  
                  S_stack_out = {new_S S_stack{:}};

            case {'best','bestdepth'}
                % parameter-dependent node search strategy
              flag_larger = false(stack_size,1);
              flag_smaller = false(stack_size,1);
              flag_devide= false(stack_size,1);
%%
              if p0.cost.Ci == inf 
                 % TODO: check: put in the end?
                 % cost of new nodes is inf (infeasible prob) --> no need to compare, put in in the end of the list
                  T_stack_out = {T_stack{:} p0, p1};
                 % update new tuple to put in S
                  new_S.Theta = Theta;
                  new_S.Tree = cell(0);
                  new_S.Tree = T_stack_out;
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % Push to S  
                  S_stack_out = {new_S S_stack{:}};
                  return
%%              
              else
                  % Determine position in T_stack where problem is inserted, according to a best first strategy                             
                  for k = 1: stack_size
                      %%
                      if T_stack{k}.cost.Ci == inf
                          % cost of new nodes is smaller, no need to compare more
                          flag_smaller(k) = true;
                          kk = k-1;
                          if kk >= stack_size  
                            T_stack_out = { T_stack{:}, p0, p1}; 
                          else
                            T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                          end                          
                          % update new tuple to put in S
                          new_S.Theta = Theta;
                          new_S.Tree = cell(0);
                          new_S.Tree = T_stack_out;
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          % Push to S  
                          %S_stack_out = {new_S S_stack{:}};
                          S_stack_out{end+1} = new_S;
                          S_stack_out = {S_stack_out{:} S_stack{:}};
                          return
                      else
                        %%
                        %[~,~,~,~,Theta_1(k),Theta_2(k)] = spatialSeparate(J_add, stack_in{k}.J_lower, Theta_add, opts);
                        if isfield(T_stack{k}.cost, 'Ai') && isfield(p0.cost, 'Ai')
                           A_diff{k} = p0.cost.Ai - T_stack{k}.cost.Ai; % for QP case 
                        end
                        B_diff{k} = p0.cost.Bi - T_stack{k}.cost.Bi;
                        C_diff{k} = p0.cost.Ci - T_stack{k}.cost.Ci;

                        if (norm(B_diff{k}) <= epsi) && (abs(C_diff{k}) <= epsi) % norm(A_diff{k}) <= epsi %for QP case
                            % J_add and Ji are equal
                            kk = k-1;
                             if kk >= stack_size  
                                T_stack_out = { T_stack{:}, p0, p1}; 
                             else
                                T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                             end
                             % update new tuple to put in S
                             new_S.Theta = Theta;
                             new_S.Tree = cell(0);
                             new_S.Tree = T_stack_out;
                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              % Push to S  
                            %S_stack_out = {new_S S_stack{:}};
                            S_stack_out{end+1} = new_S;
                            S_stack_out = {S_stack_out{:} S_stack{:}};
                            return
                        else
                            %% partition the possible region
                            H = []; K = [];
                            H = Theta.A;
                            K = Theta.b;

                            Theta1 = Polyhedron; Theta2 = Polyhedron; 

                            Theta2(k) = Polyhedron('A', [H;B_diff{k}],'b', [K;-C_diff{k}]); % J_diff <= 0 : B_diff*theta <= -C_diff
                            %Theta1(k) = cost.Theta - Theta2(k);      % J_diff >= 0 : -B_diff*theta <= C_diff
                            Theta1(k) = Polyhedron('A', [H;-B_diff{k}],'b', [K;C_diff{k}]); % J_diff >= 0
                            %Theta_int = Theta1(k) & Theta2(k); Theta_int.isEmptySet
                            %figure; Theta_int.plot;
                            
                            flag_larger(k) = ~Theta1(k).isEmptySet();
                            flag_smaller(k) = ~Theta2(k).isEmptySet();
                            flag_keep_theta1(k) = ~Theta1(k).isEmptySet() && Theta1(k).isFullDim(); %larger cost
                            flag_keep_theta2(k) = ~Theta2(k).isEmptySet() && Theta2(k).isFullDim(); %smaller cost

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
                                %_________________________________________
                            if flag_devide(k)    
                              % partition region to_add to two new regions 
                               disp('spatial splitting in Pushing to List based on the best-first Strategy');
                               p0.Theta = Polyhedron; p1.Theta = Polyhedron;  %empty the polyhedron
                               p0.Theta = Theta2; % update polyhedron
                               p1.Theta = Theta2;
                               
                               kk = k-1;
                                 if kk >= stack_size  
                                    T_stack_out = { T_stack{:}, p0, p1}; 
                                 else
                                    T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                                 end
                                 % update new tuple to put in S
                                 new_S.Theta = Theta2;
                                 new_S.Tree = cell(0);
                                 new_S.Tree = T_stack_out;
                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  % Push to S  
                                  %S_stack_out = {new_S S_stack{:}};
                                  S_stack_out{end+1} = new_S;
                            
                                  % Update remaining region
                                 Theta = Theta1; %continue further with the rest of Theta: Theta1
                                 
                                 %flag_reg_left = ~Theta.isEmptySet && Theta.isFullDim;
                                     % push remaining region in the end if in the end of list
                                     if k == stack_size % ( && flag_reg_left)
                                         %keyboard;
                                         p0.Theta = Theta; % for p0
                                         p1.Theta = Theta; %for p1
                                        %region has the worst cost
                                         T_stack_out = { T_stack{:}, p0, p1};
                                         % update new tuple to put in S
                                         new_S.Theta = Theta;
                                         new_S.Tree = cell(0);
                                         new_S.Tree = T_stack_out;
                                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                          % Push to S  
                                          %S_stack_out = {new_S S_stack{:}};
                                           S_stack_out{end+1} = new_S;
                                           S_stack_out = {S_stack_out{:} S_stack{:}};
                                           return
                                     end
                                 %_________________________________________
                            elseif ~flag_keep_theta1(k) && flag_keep_theta2(k)
                                    % the whole region should be put here
                                   kk = k-1;
                                     if kk >= stack_size  
                                        T_stack_out = { T_stack{:}, p0, p1}; 
                                     else
                                        T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                                     end
                                     % update new tuple to put in S
                                     new_S.Theta = Theta;
                                     new_S.Tree = cell(0);
                                     new_S.Tree = T_stack_out;
                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                      % Push to S  
                                      %S_stack_out = {new_S S_stack{:}};
                                       S_stack_out{end+1} = new_S;
                                       S_stack_out = {S_stack_out{:} S_stack{:}};
                                      return
                            elseif k == stack_size
                                    %region has the worst cost
                                     T_stack_out = { T_stack{:}, p0, p1}; 
                                     %kk = k;
                                     %if kk >= stack_size  
                                     %   T_stack_out = { T_stack{:}, p0, p1}; 
                                     %else
                                     %   T_stack_out = { T_stack{1:kk}, p0, p1, T_stack{kk+1:end}};
                                     %end
                                     % update new tuple to put in S
                                     new_S.Theta = Theta;
                                     new_S.Tree = cell(0);
                                     new_S.Tree = T_stack_out;
                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                      % Push to S  
                                      %S_stack_out = {new_S S_stack{:}};
                                      S_stack_out{end+1} = new_S;
                                      S_stack_out = {S_stack_out{:} S_stack{:}};
                                      return %no need
                                     
                                      %%now put the node (with the rest region) into the list
                                      %if ~p0.Theta.isEmptySet()&& p0.Theta.isFullDim()
                                       %   ii = find(~flag_larger & flag_smaller);

                                      %    if isempty(ii)
                                      %        i = stack_size;
                                      %    else
                                      %        i = ii(1)-1; %i=ii(end)+1;
                                      %    end
        %
                                       %  if i >= stack_size  
                                       %     T_stack_out = { T_stack{:}, p0}; 
                                       %  else
                                       %     T_stack_out = { T_stack{1:i}, p0, T_stack{i+1:end}};
                                      %   end
                                     % end 
                                      %
                                       % end %if costs are equal
                            end %devide
                        end %if cost are equal
                      end % if cost_k is inf
                   end   %for k=1:stack_size  
              end  %if cost is inf
      end  %switch
  end  %if stack is not empty
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Push to S  
  % S_stack_out = {new_S S_stack{:}};
  
  %%%%%%%%%%%%%%%%%%       
      %%%%%%%%%%%%%%%%%%%%%%%%
      opts.illustrate = false;
      if opts.illustrate
            figure; 
            for j=1:stack_size
                subplot(1,stack_size,j)
                T_stack{j}.Theta.plot
            end
            title('Theta os stack input')

            figure; 
            p0.Theta.plot
            title('Theta of added in stack')

            Theta_diff = Polyhedron;
            figure; title('Theta os stack input and added stack') 
            for j=1:stack_size
                subplot(2,stack_size,j)
                T_stack{j}.Theta.plot; hold on; p0.Theta.plot; hold off;
                axis([-2.5 2.5, -2.5 2.5])

                subplot(2,stack_size,stack_size+j)
                Theta_diff(j) = T_stack{j}.Theta - p0.Theta;
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

%   cost = inf;
%     switch options.method
%         case 'depth'
%             cost = 1/(subprob.level+1);
%         case 'breadth'
%             cost = subprob.level+1;
%         case 'best'
%            cost = subprob.e; % Best-first. This tends to go breadth-first
%             cost = J_lower{jreg};
%         case 'bestdepth'
%             cost = subprob.e/(subprob.level+1); %This privilegiates deep nodes
%             if isstruct(J_lower{jreg})
%                 cost =cell(0);
%                 constt = (1/(subprob.level+1));
%                 if isfield(J_lower{jreg},'Ai')
%                    cost.Ai = constt.*J_lower{jreg}.Ai; 
%                 end
%                 cost.Bi = constt.*J_lower{jreg}.Bi;
%                 cost.Ci = constt.*J_lower{jreg}.Ci;
%             else
%                 cost = (1/(subprob.level+1))*J_lower{jreg};  %zpi %This privilegiates deep nodes
%             end
%     end
%     node info for best-first search strategy
%     p0.cost = cost;  
%     p1.cost = cost;
%     p0.Theta = Theta;
%     p1.Theta = Theta;
