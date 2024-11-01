% find difference in different online and offline results

function analyze_diffs(mpQP,P_theta,th_sample,diff_iter,part,inds,AS0,x0,options,default_opts,cert_opts,iter_off,iter_on)

    [max_diff, max_diff_ind] = max(diff_iter);
    ind_diffs = find(diff_iter);
    
    [min_diff, min_diff_ind] = min(diff_iter);
    
   %______________________________________________________________________
   % First deal with negative diff iters!
    % In regions where offline counted less than online solver, solve
    % explicit sol once more and compute the differences of iters
    ind_diffs_neg = find(diff_iter < 0);
    diff_iter_neg = diff_iter(ind_diffs_neg);
    
    solve_once_more = false;
    if solve_once_more
        iter_off2 = iter_off;
        for k = 1: length(ind_diffs_neg)
            [final_part_diff_u2{k}, sol_diff_u2{k}, bnb_info_diff_u2{k},~] = Bnb_mpMIQP_V2(mpQP,part{ind_diffs_neg(k)}.Region,AS0,x0,options,default_opts,cert_opts); 
            iter_off2(ind_diffs_neg(k)) = final_part_diff_u2{k}{1}.iter;
        end
        diff_iter2 = iter_off2 -iter_on; ind_diff_iter2 = find(diff_iter2 < 0);
    end
    %______________________________________________________________________
    % Consider min case, if min is negative (min should be always zero)
    if min_diff < 0
        max_diff_ind = min_diff_ind; 
    end

    if max_diff_ind
        iter_on(max_diff_ind)
        iter_off(max_diff_ind)
        %______________________________________________________________________
        %online
        Nx = size(mpQP.H,1); 
        x0_1 = zeros(Nx,1);
        H = mpQP.H; 
        A = mpQP.A;
        Aeq = []; beq = []; %Aeq = mpQP.Aeq; beq = mpQP.beq;
        lb = mpQP.lb;   ub = mpQP.ub;
        vartype = mpQP.vartype;
        th_max_diff = th_sample(:,max_diff_ind);  %parameter point with max diffs
        f = mpQP.f + mpQP.f_theta * th_max_diff;
        b = mpQP.b + mpQP.W * th_max_diff;
        
        [xmin_on_diff, Jmin_on_diff, iter_on_diff, bnb_on_diff] = onlineMIQP_exp(H, f, A, b, Aeq, beq, vartype, lb, ub, x0_1, options,default_opts,cert_opts);
        %__________________________________________________________________
        %__________________________________________________________________
         if options.analyze2online
             % tradition online miqp
            [xmin_on_diff2, Jmin_on_diff2, iter_on_diff2, sol_on_diff2, ~, ~] = online_miqp(H, f, A, b, Aeq, beq, vartype, lb, ub, x0_1, options); %flag, Extendedflag
          %______________________________________________________________________
         else
            %offline
            if options.store_all_upperbound
                [final_part_diff_u, expl_sol_diff_u, bnb_info_off_diff_u,~] = Bnb_mpMIQP_V2(mpQP,part{max_diff_ind}.Region,AS0,x0,options,default_opts,cert_opts); 
            else
                [final_part_diff, sol_diff, bnb_info_diff,~] = Bnb_mpMIQP(mpQP,part{max_diff_ind}.Region,AS0,x0,options,default_opts,cert_opts); 
            end
            %options.solver_QPcompare =  'cplex'; %'gurobi'; % default: 'bmibnb' is a bit more precise but gurobi is so faster
            %[final_part_diff_2, sol_diff_2, bnb_info_diff_2,cert_info_diff_2] = Bnb_mpMIQP(mpQP,part{max_diff_ind}.Region,AS0,x0,options,default_opts,cert_opts); 
         end
         
    %______________________________________________________________________
        % plot regions with different iterations
        options.plot_diffs = false;
        if options.plot_diffs
            
            % first plot negative diffs           
            if ind_diffs_neg 
                inds = [1 2];               
                alpha = 1;
                values = zeros(size(mpQP.W,2)-length(inds),1);
                colors = [];
                figure; title('certified solution') 
                part_neg = part(ind_diffs_neg);
                plot_partition(part_neg,'iter',inds,colors,values,alpha,P_theta,[0,max(iter_off)]);
            end
            %______________________________________________________________
            %ind_diffs = ind_diffs_neg;
            
            % plot positive diffs
            figure; title('region with nonequal iteration');
            for k = 1: length(ind_diffs)
                Pn_diff(k) = part{ind_diffs(k)}.Region;
                if Pn_diff(k).Dim ==2
                    Pn(k) = Pn_diff(k).minHRep();        %in 2D case
                else
                    Pn(k) = Pn_diff(k).projection(inds); %project in higher dimentional parameter set
                end
                Pn(k).plot(); hold on;
            end
            axis([P_theta.lb(1) P_theta.ub(1) P_theta.lb(2) P_theta.ub(2)]);
        end
        %______________________________________________________________________
  
    end %if there is difference
%_________________________________________________________________________
 
end %function