% Find chebychev center of final partition as parameter samples

function [th_chebyCenter,iter_chebyCenter,xmin_off,Jmin_off,cheby_center] = find_chebycenter(finalpart,P_theta,exp_sol_off,options)

    N_reg = length(finalpart);
    N_theta = finalpart{1}.Region.Dim;
    
    xmin_off = cell(0);
    Jmin_off = zeros(N_reg,1);
    
    cheby_center = cell(0);
    th_chebyCenter = zeros(N_theta, N_reg);
    iter_chebyCenter = zeros(N_reg,1);
    N_poly = zeros(N_reg,1);
    
    for ireg = 1: N_reg
        N_poly(ireg) = length(finalpart{ireg}.Region);
        cheby_center{ireg} = finalpart{ireg}.Region.chebyCenter();
        th_chebyCenter(:,ireg) = cheby_center{ireg}.x;
        iter_chebyCenter(ireg) = finalpart{ireg}.iter;
    end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % find solution and value function at each sample point in explicit_sol
   th_sample = th_chebyCenter;
   
   if options.find_sol_off

           for k = 1: N_reg       
                Pn = cell(0); %Polyhedron;
                isin_ = zeros(length(exp_sol_off),1); inwhich_ = cell(0); 

                for ipart = 1 : length(exp_sol_off)%/1000
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

                [Jmin,ind_min_J] = min(J_inwhich);
                xmin = x_inwhich{ind_min_J};

                %sol_off{k}.xmin = xmin;
                %sol_off{k}.Jmin = Jmin;
                xmin_off{k} = xmin; 
                Jmin_off(k,1) = Jmin;
                
           end %for k
    end %if find offline sol at each sample point


end % main func