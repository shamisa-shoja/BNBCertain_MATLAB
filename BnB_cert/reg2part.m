% convert regions to partitions  

function part_reg =reg2part(Part)
        
    part_reg = cell(0);
    npart = length(Part.Fi);
    
    for ireg = 1: npart
        
        if isfield(Part,'Region')
            part_reg{ireg}.Region = Part.Region(ireg);
        else
            part_reg{ireg}.Region = Part.Pn(ireg);
        end
        
        if isfield(Part,'Pn')
            part_reg{ireg}.Pn = Part.Pn(ireg);
        else
            part_reg{ireg}.Pn = Part.Region(ireg);
        end
        
        if isfield(Part,'Fx')
            part_reg{ireg}.Fi = Part.Fx{ireg};
        else
            part_reg{ireg}.Fi = Part.Fi{ireg};
        end
        
        if isfield(Part,'Gx')
            part_reg{ireg}.Gi = Part.Gx{ireg};
        else
            part_reg{ireg}.Gi = Part.Gi{ireg};
        end
        
        if isfield(Part,'QJ')
            part_reg{ireg}.Ai = Part.QJ{ireg};
        else
            part_reg{ireg}.Ai = Part.Ai{ireg};
        end
        
        if isfield(Part,'RJ')
            part_reg{ireg}.Bi = Part.RJ{ireg};
        else
            part_reg{ireg}.Bi = Part.Bi{ireg};
        end
        
        if isfield(Part,'SJ')
            part_reg{ireg}.Ci = Part.SJ{ireg};
        else
            part_reg{ireg}.Ci = Part.Ci{ireg};
        end
        
        if isfield(Part,'AS')
            part_reg{ireg}.activeConstraints = Part.AS{ireg};
            part_reg{ireg}.AS = Part.AS{ireg};
        else
            part_reg{ireg}.activeConstraints = Part.activeConstraints{ireg};
        end
        
        if isfield(Part,'iter')
            if iscell(Part.iter{ireg})
                part_reg{ireg}.iter = cell2mat(Part.iter{ireg});
            else
                part_reg{ireg}.iter = Part.iter{ireg};
            end
        end
        
        if isfield(Part,'parent')
            part_reg{ireg}.parent = Part.parent{ireg}; %(ireg)
        end
        
        if isfield(Part,'state')
            part_reg{ireg}.state = Part.state{ireg};
        end
        if isfield(Part,'FLOPs')
            part_reg{ireg}.FLOPs = Part.FLOPs{ireg}; 
        end
        if isfield(Part,'feas_iter')
            if ireg <=length(Part.feas_iter) %otherwise subprob is infeasible
                part_reg{ireg}.feas_iter = Part.feas_iter{ireg};
            end
        end
        if isfield(Part,'GI_FLOPs')
            part_reg{ireg}.GI_FLOPs = Part.GI_FLOPs{ireg};
        end
        if isfield(Part,'GI_sqrts')
            part_reg{ireg}.GI_sqrts = Part.GI_sqrts{ireg};
        end
        if isfield(Part,'flags')
            part_reg{ireg}.flags = Part.flags{ireg};
        end
        if isfield(Part,'L')
            part_reg{ireg}.L = Part.L{ireg};
        end
        if isfield(Part,'D')
            part_reg{ireg}.D = Part.D{ireg};
        end
        
         if isfield(Part,'dynamics')
            part_reg{ireg}.dynamics = Part.dynamics(ireg);
         end
        
        if isfield(Part,'overlaps')
            part_reg{ireg}.overlaps = Part.overlaps; %{ireg}
        end
        
        if isfield(Part,'iter_stack')
            if iscell(Part.iter_stack{ireg})
                part_reg{ireg}.iter_stack = cell2mat(Part.iter_stack{ireg});
            else
                part_reg{ireg}.iter_stack = Part.iter_stack{ireg};
            end
        elseif isfield(Part,'iter')
            part_reg{ireg}.iter_stack = Part.iter{ireg};
        end
        
        if isfield(Part,'parent_stack')
            %if iscell(Part{ireg}.parent_stack)
            %    reg_part.parent_stack{ireg} = cell2mat(Part{ireg}.parent_stack);
            %else
                part_reg{ireg}.parent_stack = Part.parent_stack{ireg};
            %end
        elseif isfield(Part,'parent')
            part_reg{ireg}.parent_stack = Part.parent{ireg};
        end               
      
      if isfield(Part,'activeConstraints_stack')
         part_reg{ireg}.activeConstraints_stack = Part.activeConstraints_stack{ireg};
      elseif isfield(Part,'activeConstraints')
         part_reg{ireg}.activeConstraints_stack = Part.activeConstraints{ireg};          
      end
      
      if isfield(Part,'AS_log')
            part_reg{ireg}.AS_log = Part.AS_log{ireg}; 
      end
      
      if isfield(Part,'theta_cheby')
         part_reg{ireg}.theta_cheby = Part.theta_cheby{ireg};
      end
      
      if isfield(Part,'xi')
         part_reg{ireg}.xi = Part.xi{ireg};
      end
      
    end %for
    
end %function
