
function [explicit_sol_remOverlap,time_overlap_rem,dominance_cut] = reduce_overlap_sol(explicit_int_sol,relaxed_sol,compare_upp_low)
     
    % some initialization
    explicit_sol_remOverlap = cell(0);
    dominance_cut = false;
    
    t_start = cputime;
    
    %______________________________________________________________________
        % reduce overlap in upper_bound (integer solution)
    if ~compare_upp_low
        
        % specify options
        redoverlapopts.chk_part = 0;
        redoverlapopts.verbose = 0; % 2 gives lots of information
        redoverlapopts.mark_only = false;
        redoverlapopts.complex = 0;
        
        if(~isempty(explicit_int_sol)) 
            disp('Overlap removal...')
            %explicit_sol_remOverlap = reduceOverlaps_intsol(explicit_int_sol,redoverlapopts);
            
            %%use Gurobi solver in solving convex and non convex QPs
            %solve_convex_qp: CPLEX -> Gurobi, more reliable than CPLEX
            %solve_nonconvex_qp: bmibnb -> Gurobi, much faster than bmibnb
            explicit_sol_remOverlap = reduceOverlaps_gur(explicit_int_sol,redoverlapopts);
            
            if(~iscell(explicit_sol_remOverlap))
                explicit_sol_remOverlap = {explicit_sol_remOverlap};
            end
            disp('Done!')
        end        
        %______________________________________________________________________
        % reduce overlap in upper_bound (integer solution)
    else
        redoverlapopts.chk_part = 1;
        redoverlapopts.verbose = 0;
        redoverlapopts.mark_only = true;
        redoverlapopts.complex = 0; % 0 direct facet enumeration, 1 BnB
        redoverlapopts.rmv_regs = true;
        redoverlapopts.newP_only = true;
        if (~isempty(explicit_int_sol) && ~isempty(relaxed_sol))
        [~,explicit_int_sol_tmp,~,~,partitionskept,~]=reduceOverlaps_against_gur(explicit_int_sol,relaxed_sol,redoverlapopts);%x_wc_tmp
        end
        %__________________________________________________________________
        % try it once more! there are some error here, when it runs twise,
        % it differs from first time!
        check_dominance_again = false;
        if check_dominance_again
            partitionskept2 = partitionskept;
            if ~partitionskept(end)
               [~,explicit_int_sol_tmp,~,~,partitionskept,~]=reduceOverlaps_against_Sh(explicit_int_sol,relaxed_sol,redoverlapopts);%x_wc_tmp
            end
        end
        dominance_cut = ~partitionskept(end);
        explicit_sol_remOverlap = explicit_int_sol_tmp;  % to compatible with above result
        if(~iscell(explicit_sol_remOverlap))
            explicit_sol_remOverlap = {explicit_sol_remOverlap};
        end
    end
    %______________________________________________________________________
    time_overlap_rem = cputime - t_start;
    %______________________________________________________________________
end %function
        
