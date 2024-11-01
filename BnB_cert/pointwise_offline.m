
function [iter_offline,xmin_off,Jmin_off,inwhich_acc] = pointwise_offline(final_part,th_sample,exp_sol_off,options)

N_sample =  size(th_sample,2);
iter_offline  = zeros(N_sample,1);
xmin_off = cell(0); 
Jmin_off =zeros(N_sample,1);
inwhich_acc = cell(0);

%%%%%%%%%%%%%%%% find iteration number for earch sample %%%%%%%%%%%%%%%%%
    if(iscell(final_part))
          
          for k = 1 : length(final_part)
              Pnn(k) = final_part{k}.Region;
          end
          
          for i_sam = 1: N_sample
    
             inwhich = find(Pnn.contains(th_sample(:,i_sam)));                                        
             N_samp_reg = length(inwhich);
             
             if N_samp_reg
                if  N_samp_reg == 1      % && options.find_max
                    reg_th = inwhich;
                    iter_offline(i_sam) = final_part{reg_th}.iter; 
                    
                elseif options.find_max
                    %find the max one for the worst-case if theta is in more than one region
                    iter_off = zeros(N_samp_reg,1);
                    for jj = 1: N_samp_reg
                        iter_off(jj) = final_part{inwhich(jj)}.iter;
                    end
                    [~,max_iter_ind] = max(iter_off);
                    reg_th = inwhich(max_iter_ind);
                    iter_offline(i_sam) = final_part{reg_th}.iter;
                
                else
                    reg_th = inwhich(end);
                    iter_offline(i_sam) = final_part{reg_th}.iter;
                end
                inwhich_acc{i_sam} = reg_th;
              end
%                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    % find solution and value function at each sample point
%                    if options.find_sol_off       
%                        x_min{i_sam} = exp_sol_off{reg_th}.Fi{reg_th}* th_sample(:,i_sam) + exp_sol_off{reg_th}.Gi{reg_th};
%                        J_inwhich(i_sam,1) = th_sample(:,i_sam)'* exp_sol_off{i_sam}.Ai{i_sam}* th_sample(:,i_sam) + exp_sol_off{reg_th}.Bi{reg_th}* th_sample(:,i_sam) + exp_sol_off{reg_th}.Ci{reg_th};
%                    end
          end % i_sample
    else
        for k = 1 : length(final_part.Region)
              Pnn(k) = final_part.Region(k);
         end
        
        for i_sam = 1: N_sample
            %isin = any(Pnn.contains(th_sample(:,i_sam))); 
            inwhich = find(Pnn.contains(th_sample(:,i_sam))); %final_partition.Pn
            inwhich_acc{i_sam} = inwhich;
            
            if inwhich
               reg_th=inwhich(end); 
               iter_offline(i_sam) = final_part.iter{reg_th}; 
            end
        end % for i_sample
        
   end % if iscell
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % find solution and value function at each sample point
   
   if options.find_sol_off
%        %___________________________________________________________________
%           % explicit sol is already reduced overlap
%        if options.reduceOverlap_eachstep && options.save_both_Red_notRed_sol
%            for k = 1 : length(exp_sol_off)
%               Pnn2(k) = exp_sol_off{k}.Region;
%           end
%           
%           for i_sam = 1: N_sample
%              isin2 = any(Pnn2.contains(th_sample(:,i_sam)));     
%              inwhich2 = find(Pnn2.contains(th_sample(:,i_sam)));              
%              inwhich_acc2{i_sam} = inwhich2;
%              
%                 if ~any(isin2)  %&& (final_partition{region_cert}.exp_sol.Ci == inf)
%                     sol_off{i_sam}.xmin = nan;
%                     sol_off{i_sam}.Jmin = nan;
%                 else
%                     if (length(inwhich2) > 1) %&& options.find_max
%                         inwhich2 = inwhich2(end);
%                         J_inwhich(jj,1) = th_sample(:,k)'* exp_sol_off{isin_remOv(jj)}.Ai{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Bi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Ci{inwhichh_};
%                         x_inwhich{jj} = exp_sol_off{isin_remOv(jj)}.Fi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Gi{inwhichh_};
% 
%                     else
%                         J_min = th_sample(:,k)'* exp_sol_off{isin_remOv(jj)}.Ai{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Bi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Ci{inwhichh_};
%                         x_min= exp_sol_off{isin_remOv(jj)}.Fi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Gi{inwhichh_};
% 
%                     end
%                 end
%           end % i_sample           
%        %___________________________________________________________________
%        % reduced overlap is not already done in explicit solution 
%        else

           for k = 1: N_sample       
                Pn = cell(0); %Polyhedron;
                isin_ = zeros(length(exp_sol_off),1); inwhich_ = cell(0); 

                for ipart = 1 : length(exp_sol_off)
                    Pn{ipart} = exp_sol_off{ipart}.Pn;
                    isin_(ipart) = any(Pn{ipart}.contains(th_sample(:,k))); 
                    inwhich_{ipart} = find(Pn{ipart}.contains(th_sample(:,k))); %this theta is in which region in final int solution
                end

                isin_remOv = find(isin_); 
                J_inwhich = zeros(length(isin_remOv),1); x_inwhich = cell(0);
                for jj = 1: length(isin_remOv)
                    inwhichh_ = inwhich_{isin_remOv(jj)};
                    J_inwhich(jj,1) = th_sample(:,k)'* exp_sol_off{isin_remOv(jj)}.Ai{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Bi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Ci{inwhichh_};
                    x_inwhich{jj} = exp_sol_off{isin_remOv(jj)}.Fi{inwhichh_}* th_sample(:,k) + exp_sol_off{isin_remOv(jj)}.Gi{inwhichh_};
                end

                if ~isempty(J_inwhich)
                    [Jmin,ind_min_J] = min(J_inwhich);
                    xmin = x_inwhich{ind_min_J};
                else
                    Jmin = nan;
                    xmin = nan;
                end

                %sol_off{k}.xmin = xmin;
                %sol_off{k}.Jmin = Jmin;
                xmin_off{k} = xmin; 
                Jmin_off(k,1) = Jmin;
                
           end %for k
       %end %if reduced overlap explicit sol 
   end %if find offline sol at each sample point
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
    
end %main function


                
                                     
             