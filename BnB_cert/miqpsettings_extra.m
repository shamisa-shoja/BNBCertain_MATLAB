% Do some extra settings for miqp solver

function opts = miqpsettings_extra(options,default_opts,prob)
    
    opts = options;
    
    if isfield(options,'solver')
        switch options.solver
            case {'lp','linprog','lpnag','qp','quadprog','qpnag','qp_dantz'}
                opts.solver = options.solver;
            otherwise
                %warning('solver is not implemented: will take default solver')
                opts.solver = default_opts.solver;
        end
    else
        opts.solver = default_opts.solver;
    end

    if isfield(options,'mp_solver')
        switch options.mp_solver
            case {'MPT2','MPT3','MPT3_dense','MPT3_sparse'}
                opts.mp_solver = options.mp_solver;
            otherwise
                warning('mp-solver is not implemented: will take default solver')
                opts.mp_solver = default_opts.mp_solver;
        end   
    else
        opts.mp_solver = default_opts.mp_solver;   
    end

    if isfield(options,'method')
        switch options.method
            case {'depth','breadth','best','bestdepth'}
                opts.method = options.method;
            otherwise
                warning('method is not implemented: will take default method')
                opts.method = default_opts.method;
        end
    else
        opts.method = default_opts.method;
    end

    if isfield(options,'branchrule')
        switch options.branchrule
            case {'first','max','min'}
                opts.branchrule = options.branchrule;
                opts.branching_priority = [];
            case {'predefined'}
                opts.branchrule = options.branchrule;
                opts.branching_priority = options.branching_priority;
            otherwise
                warning('branchrule is not implemented: will take default branchrule')
                opts.branchrule = default_opts.branchrule;
                opts.branching_priority = [];
        end
    else
        opts.branchrule = default_opts.branchrule;
        opts.branching_priority = [];
    end

    if isfield(options,'order')
        switch options.order
            case {0,1}
                opts.order = options.order;
            case {'0','1'}
                opts.order = str2num(options.order);
            otherwise
                warning('order is not implemented: will take default order')
                opts.order = default_opts.order;
        end
    else
        opts.order = default_opts.order;
    end

    if isfield(options,'verbose')
        switch options.verbose
            case {0,1,2}
                opts.verbose = options.verbose;
            case {'0','1','2'}
                opts.verbose = str2num(options.verbose);
            otherwise
                warning('verbose is not implemented: will take default verbose')
                opts.verbose = default_opts.verbose;
        end
    else
        opts.verbose = default_opts.verbose;
    end

    if isfield(options,'maxqp')
        if options.maxqp >= 1
            opts.maxqp = floor(options.maxqp);
        else
            warning('maxqp is negative or 0: will take default maxqp')
            opts.maxqp = default_opts.maxqp;
        end
    else
        opts.maxqp = default_opts.maxqp;
    end


    if isfield(options,'inftol')
        if options.inftol > 0
            opts.inftol = options.inftol;
        else
            warning('inftol is negative or 0: will take default inftol')
            opts.inftol = default_opts.inftol;
        end
    else
        opts.inftol = default_opts.inftol;
    end

    if isfield(options,'matrixtol')
        if options.matrixtol >= 0
            opts.matrixtol = options.matrixtol;
        else
            warning('matrixtol is negative: will take default matrixtol')
            opts.matrixtol = default_opts.matrixtol;
        end
    else
        opts.matrixtol = default_opts.matrixtol;
    end

    if isfield(options,'postol')
        if options.postol >= 0
            opts.postol = options.postol;
        else
            warning('postol is negative: will take default postol')
            opts.postol = default_opts.postol;
        end
    else
        opts.postol = default_opts.postol;
    end

    if isfield(options,'integtol')
        if options.integtol >= 0
            opts.integtol = options.integtol;
        else
            warning('integtol is negative: will take default integtol')
            opts.integtol = default_opts.integtol;
        end
    else
        opts.integtol = default_opts.integtol;
    end

    if isfield(options,'maxQPiter')
        if options.maxQPiter >= 1
            opts.maxQPiter = floor(options.maxQPiter);
        else
            warning('maxQPiter is negative or 0: will take default maxQPiter')
            opts.maxQPiter = default_opts.maxQPiter;
        end
    else
        opts.maxQPiter = default_opts.maxQPiter;
    end

    if isfield(options,'deletefile')
        if options.deletefile == 1
            opts.deletefile = 1;
        elseif options.deletefile == 0
            opts.deletefile = 0;
        else
            warning('unallowed value for deletefile: will take default')
            opts.deletefile = default_opts.deletefile;
        end
    else
        opts.deletefile = default_opts.deletefile;
    end

    if isfield(options,'round')
        if options.round == 1
            opts.rounding = 1;
        elseif options.round == 0
            opts.rounding = 0;
        else
            warning('unallowed value for round: will take default')
            opts.rounding = default_opts.round;
        end
    else
        opts.rounding = default_opts.round;
    end

    if isfield(options,'zstar')
        opts.zstar_cut = options.zstar;
    else
        opts.zstar_cut = inf;
    end

    if isfield(options,'zstar_fixed')
        opts.zstar_cut_fixed = options.zstar_fixed;
    else
        opts.zstar_cut_fixed = true;
    end

    if isfield(options,'diff_QP_solver')
        opts.diff_QP_solver = options.diff_QP_solver;
    else
        opts.diff_QP_solver = 'bmibnb';
    end

    if isfield(options,'diff_QP_presolve')
        opts.diff_QP_presolve = options.diff_QP_presolve;
    else
        opts.diff_QP_presolve = 'bmibnb';
    end

    if isfield(options,'store_expl_sols')
        opts.store_expl_sols = options.store_expl_sols;
    else
        opts.store_expl_sols = false;
    end

    if isfield(options,'explicit_eval') % -1 feasibility, 0 MI, 1 explicit, 2 adaptive MI/explicit
        opts.explicit_eval = options.explicit_eval;
    else
        opts.explicit_eval = 1;
    end

    if isfield(options,'max_explore_depth')
        opts.max_explore_depth = options.max_explore_depth;
    else
        opts.max_explore_depth = inf;
    end

    if isfield(options,'max_expl_explore_depth')
        opts.max_expl_explore_depth = options.max_expl_explore_depth;
    else
        opts.max_expl_explore_depth = inf;
    end

    if isfield(options,'integer_times')
        opts.integer_times = options.integer_times(:);
    else
        if strcmp(opts.method,'keep_order')
            warning('No binary variable order is given. Binary variables are assigned an enumerated number.')
        end

        opts.integer_times = (1:length(vartype))';  
    end

    if isfield(options,'reltol')
        opts.reltol = options.reltol;
    else
        opts.reltol = 0;
    end

    if isfield(options,'abstol')
        opts.abstol = options.abstol;
    else
        opts.abstol = 0;
    end

    if isfield(options,'quadtol')
        opts.quadtol = options.quadtol;
    elseif isfield(prob,'W')
        opts.quadtol = zeros(size(prob.W,2));
    end

    if isfield(options,'maxtime')
        opts.maxtime = options.maxtime;
        opts.starttime = cputime;
    else
        opts.maxtime = inf;
        opts.starttime = 0;
    end

    if isfield(options,'add_complete_int_seq_sols')
        opts.add_complete_int_seq_sols = options.add_complete_int_seq_sols;
    else
        opts.add_complete_int_seq_sols = false; %true
    end

    if isfield(options,'tighten_feasible_set_limit')
        opts.tighten_feasible_set_limit = options.tighten_feasible_set_limit;
    else
        opts.tighten_feasible_set_limit = 0; % Disabled
    end

    if isfield(options,'spat_decomp_hole_ratio')
        opts.spat_decomp_hole_ratio = options.spat_decomp_hole_ratio;
    else
        opts.spat_decomp_hole_ratio = inf;
    end

    if isfield(options,'postpone_overlap_rem_to_end')
        opts.postpone_overlap_rem_to_end = options.postpone_overlap_rem_to_end;
    else
        opts.postpone_overlap_rem_to_end = true;
    end

    if isfield(options,'add_glob_min_sol')
        opts.add_glob_min_sol = options.add_glob_min_sol;
    else
        opts.add_glob_min_sol = true;
    end


    if ~isfield(options,'solver')
        if max(svd(prob.H)) < opts.matrixtol
            if options.verbose >= 1
                warning('This is a MILP')
            end
            opts.solver = default_opts.solverlp;  % default solver for LPs
        end
    elseif strcmp(options.solver,'qp_dantz')
        if max(svd(prob.H)) < opts.matrixtol
            error('MILPs cannot be solved with qp_dantz')
        end
    end
    
    if((opts.explicit_eval > 0) && opts.store_expl_sols)
        error('store_expl_sols not supported')
    else
        opts.low_bound_explicit = [];
        opts.up_bound_explicit = [];
    end

    if isfield(options,'acc_int_problems')
        opts.acc_int_problems = options.acc_int_problems;
        opts.acc_int_problems_precomputed = true;
    else
        opts.acc_int_problems = [];
        opts.acc_int_problems_precomputed = false;
    end

    if isfield(options,'acc_int_regions')
        opts.nr_of_only_int_sol_track = length(options.acc_int_problems);
        opts.acc_int_regions = options.acc_int_regions;
    else
        opts.acc_int_regions = {};
    end

    if isfield(options,'acc_int_seqs')
        opts.acc_int_seqs = options.acc_int_seqs;
    else
        opts.acc_int_seqs = [];
    end

    if isfield(options,'check_subtree_int_feas')
        opts.check_subtree_int_feas = options.check_subtree_int_feas;
    else
        opts.check_subtree_int_feas = true;
    end

    if ~isfield(options,'search_type')
        opts.search_type= 'fifo'; %'lifo';
    end
    
    opts.search_type = 'fifo';
    opts.cons_bin = true; % consider binary variables!
    opts.problem_type = [];
    
end  % main function