
function [detail_MIP] = MIP_online(prob,P_theta,AS0,x0,final_part,explicit_sol, bnb_info_off,options,default_opts,cert_opts,nr_sample,th_points)
    
    %%Find lower and upper bound of complexity certification
    options.find_upp_low_iter = true;
    if options.find_upp_low_iter    
        [iter_off,max_iter_off] = struct2vec(final_part,'iter'); 
        detail_MIP.iter_off = iter_off;
        detail_MIP.max_iter_off = max_iter_off;
    end
%% Generate random parameters in parameter set to evaluate results
    if options.gen_rand_samples
        %nr_sample = 500;
        sample_kind = 'deter';   % sampling kinds: 
                                 % 'deter' : deterministic: 
                                 % dimention of theta is (N_theta, N_sample^N_theta)
                                 % 'rand'  : randoms:
                                 % dimention of theta is (N_theta, N_sample)
        if isfield(P_theta,'lb')
            th_sample = sampling(nr_sample,P_theta.lb,P_theta.ub,sample_kind);
        else 
            % use boundong box as lower and upper bounds
            box = P_theta.outerApprox();
            box_lb = box.Internal.lb;
            box_ub = box.Internal.ub;
            th_sample = sampling(nr_sample,box_lb,box_ub,sample_kind);
        end
        

        %% calculate iteration number of explicit (offline) solution
        options.find_sol_off = false;   %find solution at each sample point
        [iter_offline,xmin_off,Jmin_off,~] = pointwise_offline(final_part,th_sample,explicit_sol,options); %inwhich_region; %iter_offline1 = pointwize_iter_offline(final_partition,th_sample,false); %(:,k)

        %% Online MIQP
        options.detail_online = true;     % get all info of online solver
        options.onlineMiqp_trad = false;  % true: apply tradition online_miqp code as well
        t_s = cputime;
        [iter_online,info_online,iter_online2,info_online2] = compute_online_mip(prob,th_sample,options,default_opts,cert_opts);
        detail_MIP.samples.time = cputime - t_s;
        
        %%Difference between online and offline complexity certification
        % difference with written online miqp iteration
        [diff_iter,max_diff,min_diff,mean_diff,relat_err] = find_vec_diffs(iter_offline,iter_online);
        diff_iter_ind = find(diff_iter); err= length(diff_iter_ind)/ length(diff_iter)

        if options.onlineMiqp_trad
            [diff_iter2,max_diff2,min_diff2,mean_diff2] = find_vec_diffs(iter_offline,iter_online2);
            % difference between two online miqp iteration
            [diff_iter_on,max_diff_on,~,mean_diff_on] = find_vec_diffs(iter_online,iter_online2); %min_diff_on
        end
    end
    %__________________________________________________________________________
%% Chebychev centers as samples 
    if options.gen_chebyCenters
        %Find chebychev center of all regions in final partition
        [th_chebyCenter,iter_off_chCen,xmin_off_chCen,Jmin_off_chCen,~] = find_chebycenter(final_part,P_theta,explicit_sol,options); %cheby_center
        N_chebyCenter = size(th_chebyCenter,2);

        %%Online Certified MIP in chebychev centers
        %options.detail_online = true;     % get all info of online solver
        t_s2 = cputime;
        [iter_online_chCen,info_online_chC,iter_online2_chCen,info_online2_chC] = compute_online_mip(prob,th_chebyCenter,options,default_opts,cert_opts);
        detail_MIP.chebyCen.time = cputime - t_s2;
        
        %%Difference of iteration number in chebycenter online and offline
        [diff_iter_chCen,max_diff_chCen,min_diff_chCen,mean_diff_chCen,relat_err_chCen] = find_vec_diffs(iter_off_chCen,iter_online_chCen);
        diff_iter_chCen_ind = find(diff_iter_chCen); err_ch=length(diff_iter_chCen_ind)/ length(diff_iter_chCen)

        % difference with tradition online miqp 
        if options.onlineMiqp_trad
            [diff_iter2_chCen,max_diff2_chCen,min_diff2_chCen,mean_diff2_chCen] = find_vec_diffs(iter_off_chCen,iter_online2_chCen);
            [diff_iter_on_chCen,max_diff_on_chCen,min_diff_on_chCen,mean_diff_on_chCen] = find_vec_diffs(iter_online_chCen,iter_online2_chCen); % difference of two online miqp 
        end
    end % gen_chebyCenters
     %__________________________________________________________________________
%% Inputted points
    if options.existed_point
        N_chebyCenter = size(th_points,2);
        
        % calculate iteration number of explicit (offline) solution
        options.find_sol_off = false;   %find solution at each sample point
        [iter_off_points,xmin_off_points,Jmin_off_points,~] = pointwise_offline(final_part,th_points,explicit_sol,options); %inwhich_region; %iter_offline1 = pointwize_iter_offline(final_partition,th_sample,false); %(:,k)

        %%Online Certified MIP in chebychev centers
        %options.detail_online = true;     % get all info of online solver
        t_s2 = cputime;
        [iter_online_points,info_online_points,iter_online2_points,info_online2_points] = compute_online_mip(prob,th_points,options,default_opts,cert_opts);
        detail_MIP.chebyCen.time = cputime - t_s2;
        
        %%Difference of iteration number in chebycenter online and offline
        [diff_iter_points,max_diff_points,min_diff_points,mean_diff_points,relat_err_points] = find_vec_diffs(iter_off_points,iter_online_points);
        diff_iter_points_ind = find(diff_iter_points); err_points=length(diff_iter_points_ind)/ length(diff_iter_points)

        % difference with tradition online miqp 
        if options.onlineMiqp_trad
            [diff_iter2_points,max_diff2_points,min_diff2_points,mean_diff2_points] = find_vec_diffs(iter_off_points,iter_online2_points);
            [diff_iter_on_points,max_diff_on_points,min_diff_on_points,mean_diff_on_points] = find_vec_diffs(iter_online_points,iter_online2_points); % difference of two online miqp 
        end
    end % existed_point  
    %% Evaluate result: Solve MIQP problem (not Certify- not parametrically)
    options.online_miqp_gur = false;
    if options.online_miqp_gur
        cons_bin = true;   %solve MIQP problem not QP
        %______________________________________________________________________
        %at random points
        [xmin_gur,Jmin_gur,~,sol_gur]=sol_gurobi(prob,th_sample,cons_bin);
        % difference of optimization variables in certified online miqp and online np-miqp (not-param)
        [diff_x_gur,max_diff_x_gur,min_diff_x_gur,mean_diff_x_gur] = find_vec_diffs(info_online.x,xmin_gur);
        % difference of value function in different online solvers
        [diff_J_gur,max_diff_J_gur,min_diff_J_gur,mean_diff_J_gur] = find_vec_diffs(info_online.J,Jmin_gur');
        %______________________________________________________________________
        % at chebychev centers
        [xmin_chCen_gur,Jmin_chCen_gur,~,sol_chCen_gur]=sol_gurobi(prob,th_chebyCenter,cons_bin);
        % difference of optimization variables in different online solvers
        [diff_x_gur_chCen,max_diff_x_gur_chCen,min_diff_x_gur_chCen,mean_diff_x_gur_chCen] = find_vec_diffs(info_online_chC.x,xmin_chCen_gur);
        % difference of value function in different online solvers
        [diff_J_gur_chCen,max_diff_J_gur_chCen,min_diff_J_gur_chCen,mean_diff_J_gur_chCen] = find_vec_diffs(info_online_chC.J,Jmin_chCen_gur');
    end
    
    %% Analyze difference between online and offline
    options.try_find_diff_on_off = false; %Dig into where there are difference and why
    if options.try_find_diff_on_off
        analyze_diffs(prob,P_theta,th_chebyCenter,diff_iter_chCen,final_part,inds,AS0,x0,options,default_opts,cert_opts);    
    end

    %%Analyze difference between two online miqp solver
    options.find_online_diffs = false;
    if options.onlineMiqp_trad && options.find_online_diffs    
        options.analyze2online = true;  %find diffs in 2 online results
        analyze_diffs(prob,P_theta,th_chebyCenter,diff_iter_on_chCen,final_part,inds,AS0,x0,options,default_opts,cert_opts)
    end 
    
%% Algorithm 3: apply reduce overlap function in QPcomparision
    if options.apply_red_overlap_func
        time_cmpQP = (bnb_info_off.time_compQP)/60;
        time_not_cmpQP = (bnb_info_off.time_bnb - bnb_info_off.time_compQP)/60;
        % in case of reducing overlap in upper bound sol at each iteration
        options.reduceOverlap_eachstep = false;
        if options.reduceOverlap_eachstep 
            [final_partition_u_redeach, explicit_sol_u_redeach, bnb_info_off_u_redeach,cert_info_off_u_redeach] = Bnb_mpMIQP_V2(prob,P_theta,AS0,x0,options,default_opts,cert_opts); 
        end
    end
    
    %% Root node: fully relaxed QP
 options.find_root_node = false;
 if options.find_root_node
     root_node = cert_info_off.part_detail{1};  %bnb_info_cert.seq_all_regions{1};
     %Find lower and upper iteration number of root node
     [iter_root,max_iter_root] = struct2vec(final_part,'iter');    

    %%Plot results
    inds = [1,2];                  % index of parameter to plot based on
    options.plot_allreg = true;    % plot all certified regions
    options.plot_points = false;   % specify parameters in figure
    options.plot3_iter = false;    % plot iterations in parameter space in 3D
    plot_parts(prob,P_theta,root_node,max_iter_root,inds,options,[],[],[],[])
 end
 
    %% Enumerate all nodes and certify MIQP
    %%solve enumerated certified mpMIQP problem
    options.try_enum_miqp = false;
    if options.try_enum_miqp
        options.enum_check_inf = true;   %true: in enumeration case, cut when reached to infeasible node
        [final_part_enum, expl_sol_enum, bnb_info_enum,th_chCen_enum,iter_enum_chCen_off,iter_enum_chCen_on,info_enum_on] = enum_bnb(prob,P_theta,AS0,x0,options,default_opts,cert_opts);
        % Find difference of iteration number in chebycenter online and offline
        [diff_iter_chCen_enum,max_diff_chCen_enum,min_diff_chCen_enum,mean_diff_chCen_enum] = find_vec_diffs(iter_enum_chCen_off,iter_enum_chCen_on);
    end

    %%Analyze difference in enumeration case
    options.find_enum_diffs = false;
    if options.find_enum_diffs
        options.analyze2online = false;  %find diffs in offline and online results
        analyze_diffs(prob,P_theta,th_chCen_enum,diff_iter_chCen_enum,final_part,inds,AS0,x0,options,default_opts,cert_opts);    
    end 
    %%
    options.lowmemory = true;
    if options.lowmemory
       subprob = cell(0); %erase to have less memory required
       cert_info_off = cell(0);
       %bnb_info_off = cell(0);
    end
    
    %% Plot results
    options.plot_result = false;
    if options.plot_result
        %close all
        inds = [1,2];                  % index of parameter to plot based on
        options.plot_allreg = true;    % plot all certified regions
        options.plot_points = options.plot_allreg;    % specify parameters in figure; %th_chebyCenter %th_sample
        options.plot3_iter = ~options.plot_allreg;    % plot iterations in parameter space in 3D
        if options.gen_chebyCenters
            plot_parts(prob,P_theta,final_part,max(iter_off_chCen),inds,options,th_chebyCenter,iter_off_chCen,iter_online_chCen,diff_iter_chCen)
        elseif options.gen_rand_samples
            plot_parts(prob,P_theta,final_part,max_iter_off,inds,options,th_sample,iter_offline,iter_online,diff_iter)
        end
    end
    %% Save results
    % result of sample point
    if options.gen_rand_samples
        detail_MIP.samples.th_sample = th_sample;
        detail_MIP.samples.iter_off = iter_offline;
        detail_MIP.samples.xmin_off = xmin_off;
        detail_MIP.samples.Jmin_off = Jmin_off;

        detail_MIP.samples.iter_on = iter_online;
        detail_MIP.samples.xmin_on = info_online.x;
        detail_MIP.samples.Jmin_on = info_online.J;
        detail_MIP.samples.bnb_info_on = info_online.bnb_info;

        detail_MIP.samples.diff_iter = diff_iter;
        detail_MIP.samples.Nr_diffs_iter = length(find(diff_iter));
        detail_MIP.samples.max_diff_iter = max_diff;
        detail_MIP.samples.min_diff_iter = min_diff;
        detail_MIP.samples.mean_diff_iter = mean_diff;
        detail_MIP.samples.err_diff_iter = err;
    end
    % result in chebychevcenters
    if options.gen_chebyCenters
        detail_MIP.chebyCen.th_chebyCenter = th_chebyCenter;
        detail_MIP.chebyCen.iter_off = iter_off_chCen;
        detail_MIP.chebyCen.xmin_off = xmin_off_chCen;
        detail_MIP.chebyCen.Jmin_off = Jmin_off_chCen;

        detail_MIP.chebyCen.iter_on = iter_online_chCen;
        detail_MIP.chebyCen.xmin_on = info_online_chC.x;
        detail_MIP.chebyCen.Jmin_on = info_online_chC.J;
        detail_MIP.chebyCen.bnb_info_on = info_online_chC.bnb_info;

        detail_MIP.chebyCen.diff_iter = diff_iter_chCen;
        detail_MIP.chebyCen.max_diff_iter = max_diff_chCen;
        detail_MIP.chebyCen.min_diff_iter = min_diff_chCen;
        detail_MIP.chebyCen.mean_diff_iter = mean_diff_chCen;
        detail_MIP.chebyCen.err_diff_iter = err_ch;
    end
    
   
    % result of existed_point
    if options.existed_point
        detail_MIP.points.th_sample = th_points;
        detail_MIP.points.iter_off = iter_off_points;
        detail_MIP.points.xmin_off = xmin_off_points;
        detail_MIP.points.Jmin_off = Jmin_off_points;

        detail_MIP.points.iter_on = iter_online_points;
        detail_MIP.points.xmin_on = info_online_points.x;
        detail_MIP.points.Jmin_on = info_online_points.J;
        detail_MIP.points.bnb_info_on = info_online_points.bnb_info;

        detail_MIP.points.diff_iter = diff_iter_points;
        detail_MIP.points.Nr_diffs_iter = length(find(diff_iter_points));
        detail_MIP.points.max_diff_iter = max_diff_points;
        detail_MIP.points.min_diff_iter = min_diff_points;
        detail_MIP.points.mean_diff_iter = mean_diff_points;
        detail_MIP.points.err_diff_iter = err_points;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % root node 
    if options.find_root_node
        detail_MIP.rootnode.root_node = root_node;
        detail_MIP.rootnode.iter_root = iter_root;
    end
    
    % enumeration case
    if options.try_enum_miqp
        detail_MIP.enum.final_part_enum = final_part_enum;
        detail_MIP.enum.expl_sol_enum = expl_sol_enum;
        detail_MIP.enum.bnb_info_enum = bnb_info_enum;
        detail_MIP.enum.th_chCen_enum = th_chCen_enum;
        detail_MIP.enum.iter_enum_off_chCen = iter_enum_chCen_off;
        detail_MIP.enum.iter_enum_on_chCen = iter_enum_chCen_on;
        detail_MIP.enum.info_enum_on = info_enum_on;
        detail_MIP.enum.diff_iter = diff_iter_chCen_enum;
        detail_MIP.enum.max_diff_iter = max_diff_chCen_enum;
        detail_MIP.enum.min_diff_iter = min_diff_chCen_enum;
        detail_MIP.enum.mean_diff_iter = mean_diff_chCen_enum;
    end
    
end %function