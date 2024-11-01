% convert partitions to regions 

function reg_part =part2reg(Part)
        
    reg_part = cell(0);
    npart = length(Part);
    
    for ireg = 1: npart
      if iscell(Part) 
        if isfield(Part{ireg},'Region')
            reg_part.Region(ireg) = Part{ireg}.Region;
        else
            reg_part.Region(ireg) = Part{ireg}.Pn;
        end
        
        if isfield(Part{ireg},'Pn')
            reg_part.Pn(ireg) = Part{ireg}.Pn;
        else
            reg_part.Pn(ireg) = Part{ireg}.Region;
        end
        
        if isfield(Part{ireg},'Fx')
            reg_part.Fi{ireg} = Part{ireg}.Fx;
        elseif isfield(Part{ireg},'Fi')
            reg_part.Fi{ireg} = Part{ireg}.Fi;
        end
        
        if isfield(Part{ireg},'Gx')
            reg_part.Gi{ireg} = Part{ireg}.Gx;
        elseif isfield(Part{ireg},'Gi')
            reg_part.Gi{ireg} = Part{ireg}.Gi;
        end
        
        if isfield(Part{ireg},'QJ')
            reg_part.Ai{ireg} = Part{ireg}.QJ;
        elseif isfield(Part{ireg},'Ai')
            reg_part.Ai{ireg} = Part{ireg}.Ai;
        end
        
        if isfield(Part{ireg},'RJ')
            reg_part.Bi{ireg} = Part{ireg}.RJ;
        else
            reg_part.Bi{ireg} = Part{ireg}.Bi;
        end
        
        if isfield(Part{ireg},'SJ')
            reg_part.Ci{ireg} = Part{ireg}.SJ;
        else
            reg_part.Ci{ireg} = Part{ireg}.Ci;
        end
        
        if isfield(Part{ireg},'AS')
            reg_part.activeConstraints{ireg} = Part{ireg}.AS;
            reg_part.AS{ireg} = Part{ireg}.AS;
        else
            reg_part.activeConstraints{ireg} = Part{ireg}.activeConstraints;
        end
        
        if isfield(Part{ireg},'iter')
            if iscell(Part{ireg}.iter)
                reg_part.iter{ireg} = cell2mat(Part{ireg}.iter);
            else
                reg_part.iter{ireg} = Part{ireg}.iter;
            end
        end
        
        if isfield(Part{ireg},'parent')
            reg_part.parent{ireg} = Part{ireg}.parent; %(ireg)
        end
        
        if isfield(Part{ireg},'state')
            reg_part.state{ireg} = Part{ireg}.state;
        end
        if isfield(Part{ireg},'FLOPs')
            reg_part.FLOPs{ireg} = Part{ireg}.FLOPs; 
        end
        if isfield(Part{ireg},'feas_iter')
            reg_part.feas_iter{ireg} = Part{ireg}.feas_iter;
        end
        if isfield(Part{ireg},'GI_FLOPs')
            reg_part.GI_FLOPs{ireg} = Part{ireg}.GI_FLOPs;
        end
        if isfield(Part{ireg},'GI_sqrts')
            reg_part.GI_sqrts{ireg} = Part{ireg}.GI_sqrts;
        end
        if isfield(Part{ireg},'flags')
            reg_part.flags{ireg} = Part{ireg}.flags;
        end
        if isfield(Part{ireg},'L')
            reg_part.L{ireg} = Part{ireg}.L;
        end
        if isfield(Part{ireg},'D')
            reg_part.D{ireg} = Part{ireg}.D;
        end
        
         if isfield(Part{ireg},'dynamics')
            reg_part.dynamics{ireg} = Part{ireg}.dynamics;
         end
        
        if isfield(Part{ireg},'overlaps')
            reg_part{ireg}.overlaps = Part.overlaps{ireg};
        end
        
        if isfield(Part{ireg},'iter_stack')
            if iscell(Part{ireg}.iter_stack)
                reg_part.iter_stack{ireg} = cell2mat(Part{ireg}.iter_stack);
            else
                reg_part.iter_stack{ireg} = Part{ireg}.iter_stack;
            end
        elseif isfield(Part{ireg},'iter')
            reg_part.iter_stack{ireg} = Part{ireg}.iter;
        end
        
        if isfield(Part{ireg},'parent_stack')
            %if iscell(Part{ireg}.parent_stack)
            %    reg_part.parent_stack{ireg} = cell2mat(Part{ireg}.parent_stack);
            %else
                reg_part.parent_stack{ireg} = Part{ireg}.parent_stack;
            %end
        elseif isfield(Part{ireg},'parent')
            reg_part.parent_stack{ireg} = Part{ireg}.parent;
        end               
      
      if isfield(Part{ireg},'activeConstraints_stack')
         reg_part.activeConstraints_stack{ireg} = Part{ireg}.activeConstraints_stack;
      elseif isfield(Part{ireg},'activeConstraints')
         reg_part.activeConstraints_stack{ireg} = Part{ireg}.activeConstraints;          
      end
      
      if isfield(Part{ireg},'AS_log')
            reg_part.AS_log{ireg} = Part{ireg}.AS_log; 
      end
      
      if isfield(Part{ireg},'theta_cheby')
         reg_part.theta_cheby{ireg} = Part{ireg}.theta_cheby;
      end
      
      if isfield(Part{ireg},'xi')
         reg_part.xi{ireg} = Part{ireg}.xi;
      end
    %______________________________________________________________________
        %if Part is not cell one part and one reg
      else
            if isfield(Part,'Region')
                reg_part.Region(ireg) = Part.Region;
            else
                reg_part.Region(ireg) = Part.Pn;
            end

            if isfield(Part,'Pn')
                reg_part.Pn(ireg) = Part.Pn;
            else
                reg_part.Pn(ireg) = Part.Region;
            end

            if isfield(Part,'Fx')
                reg_part.Fi{ireg} = Part.Fx;
            else
                reg_part.Fi{ireg} = Part.Fi;
            end

            if isfield(Part,'Gx')
                reg_part.Gi{ireg} = Part.Gx;
            else
                reg_part.Gi{ireg} = Part.Gi;
            end

            if isfield(Part,'QJ')
                reg_part.Ai{ireg} = Part.QJ;
            else
                reg_part.Ai{ireg} = Part.Ai;
            end

            if isfield(Part,'RJ')
                reg_part.Bi{ireg} = Part.RJ;
            else
                reg_part.Bi{ireg} = Part.Bi;
            end

            if isfield(Part,'SJ')
                reg_part.Ci{ireg} = Part.SJ;
            else
                reg_part.Ci{ireg} = Part.Ci;
            end

            if isfield(Part,'AS')
                reg_part.activeConstraints{ireg} = Part.AS;
                reg_part.AS{ireg} = Part.AS;
            else
                reg_part.activeConstraints{ireg} = Part.activeConstraints;
            end

            if isfield(Part,'iter')
                if iscell(Part.iter)
                    reg_part.iter{ireg} = cell2mat(Part.iter);
                else
                    reg_part.iter{ireg} = Part.iter;
                end
            end

            if isfield(Part,'parent')
                reg_part.parent{ireg} = Part.parent; %(ireg)
            end

            if isfield(Part,'state')
                reg_part.state{ireg} = Part.state;
            end
            if isfield(Part,'FLOPs')
                reg_part.FLOPs{ireg} = Part.FLOPs; 
            end
            if isfield(Part,'feas_iter')
                reg_part.feas_iter{ireg} = Part.feas_iter;
            end
            if isfield(Part,'GI_FLOPs')
                reg_part.GI_FLOPs{ireg} = Part.GI_FLOPs;
            end
            if isfield(Part,'GI_sqrts')
                reg_part.GI_sqrts{ireg} = Part.GI_sqrts;
            end
            if isfield(Part,'flags')
                reg_part.flags{ireg} = Part.flags;
            end
            if isfield(Part,'L')
                reg_part.L{ireg} = Part.L;
            end
            if isfield(Part,'D')
                reg_part.D{ireg} = Part.D;
            end

             if isfield(Part,'dynamics')
                reg_part.dynamics{ireg} = Part.dynamics;
             end

            if isfield(Part,'overlaps')
                reg_part{ireg}.overlaps = Part.overlaps{ireg};
            end

            if isfield(Part,'iter_stack')
                if iscell(Part.iter_stack)
                    reg_part.iter_stack{ireg} = cell2mat(Part.iter_stack);
                else
                    reg_part.iter_stack{ireg} = Part.iter_stack;
                end
            elseif isfield(Part,'iter')
                reg_part.iter_stack{ireg} = Part.iter;
            end

            if isfield(Part,'parent_stack')
                %if iscell(Part.parent_stack)
                %    reg_part.parent_stack{ireg} = cell2mat(Part.parent_stack);
                %else
                    reg_part.parent_stack{ireg} = Part.parent_stack;
                %end
            elseif isfield(Part,'parent')
                reg_part.parent_stack{ireg} = Part.parent;
            end               

          if isfield(Part,'activeConstraints_stack')
             reg_part.activeConstraints_stack{ireg} = Part.activeConstraints_stack;
          elseif isfield(Part,'activeConstraints')
             reg_part.activeConstraints_stack{ireg} = Part.activeConstraints;          
          end

          if isfield(Part,'AS_log')
                reg_part.AS_log{ireg} = Part.AS_log; 
          end

          if isfield(Part,'theta_cheby')
             reg_part.theta_cheby{ireg} = Part.theta_cheby;
          end

          if isfield(Part,'xi')
             reg_part.xi{ireg} = Part.xi;
          end
        
      end %if Part is not cell
      %____________________________________________________________________
      
    end %for
    
    reg_part.Pfinal = PolyUnion(reg_part.Pn(ireg));
    
end %function
