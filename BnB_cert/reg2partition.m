% convert regions to partitions

function reg_part =reg2partition(Part)
        
    reg_part = cell(0);
    nreg = length(Part.Fi);
    for ireg = 1: nreg
        if isfield(Part,'Region')
            reg_part{ireg}.Region = Part.Region(ireg);
        else
            reg_part{ireg}.Region = Part.Pn(ireg);
        end
        reg_part{ireg}.Pn = Part.Pn(ireg);
        reg_part{ireg}.Fi = Part.Fi{ireg};
        reg_part{ireg}.Gi = Part.Gi{ireg};
        reg_part{ireg}.Ai = Part.Ai{ireg};
        reg_part{ireg}.Bi = Part.Bi{ireg};
        reg_part{ireg}.Ci = Part.Ci{ireg};
        
        reg_part{ireg}.activeConstraints = Part.activeConstraints{ireg};
        
        if isfield(Part,'dynamics')
            if ~isempty(Part.dynamics)
                reg_part{ireg}.dynamics = Part.dynamics{ireg};
            else
                reg_part{ireg}.dynamics = [];
            end
        end
        if isfield(Part,'overlaps')
            if ~isempty(Part.overlaps)
                reg_part{ireg}.overlaps = Part.overlaps{ireg};
            else
                reg_part{ireg}.overlaps = [];
            end
        end
        if isfield(Part,'parent')
            reg_part{ireg}.parent = Part.parent(ireg);
        end
%         if isfield(Part,'iter')
%             reg_part{ireg}.iter = Part.iter(ireg);
%         end
        if isfield(Part,'iter')
            if iscell(Part.iter{ireg})
                reg_part{ireg}.iter = cell2mat(Part.iter{ireg});
            else
                reg_part{ireg}.iter = Part.iter{ireg};
            end
        end
        if isfield(Part,'state')
            reg_part{ireg}.state = Part.state(ireg);
        end
        if isfield(Part,'FLOPs')
            reg_part{ireg}.FLOPs = Part.FLOPs{ireg}; 
        end
        if isfield(Part,'feas_iter')
            reg_part{ireg}.feas_iter = Part.feas_iter{ireg};
        end
        if isfield(Part,'GI_FLOPs')
            reg_part{ireg}.GI_FLOPs = Part.GI_FLOPs{ireg};
        end
        if isfield(Part,'GI_sqrts')
            reg_part{ireg}.GI_sqrts = Part.GI_sqrts{ireg};
        end
        if isfield(Part,'flags')
            reg_part{ireg}.flags = Part.flags{ireg};
        end
        if isfield(Part,'L')
            reg_part{ireg}.L = Part.L{ireg};
        end
        if isfield(Part,'D')
            reg_part{ireg}.D = Part.D{ireg};
        end
        if isfield(Part,'Pfinal')
            reg_part{ireg}.Pfinal = Part.Pn(ireg);
        end
        
        if isfield(Part,'iter_intsct')
            if iscell(Part.iter_intsct{ireg})
                reg_part{ireg}.iter_intsct = cell2mat(Part.iter_intsct{ireg});
            else
                reg_part{ireg}.iter_intsct = Part.iter_intsct{ireg};
            end
        end
        
%         if isfield(Part,'iter_intsct')
%             if iscell(Part.iter_intsct{ireg})
%                 reg_part{ireg}.iter_intsct = cell2mat(Part.iter_intsct{ireg});
%             else
%                 reg_part{ireg}.iter_intsct = Part.iter_intsct{ireg};
%             end
%         end
        
       if isfield(Part,'iter_intsct_stack2')
            if iscell(Part.iter_intsct_stack2{ireg})
                reg_part{ireg}.iter_intsct_stack2 = cell2mat(Part.iter_intsct_stack2{ireg});
            else
                reg_part{ireg}.iter_intsct_stack2 = Part.iter_intsct_stack2{ireg};
            end
       end
       
       if isfield(Part,'iter_intsct_stack')
            if iscell(Part.iter_intsct_stack{ireg})
                reg_part{ireg}.iter_intsct_stack = cell2mat(Part.iter_intsct_stack{ireg});
            else
                reg_part{ireg}.iter_intsct_stack = Part.iter_intsct_stack{ireg};
            end
       end
              
       if isfield(Part,'parent_intsct')
            if iscell(Part.parent_intsct{ireg})
                reg_part{ireg}.parent_intsct = cell2mat(Part.parent_intsct{ireg});
            else
                reg_part{ireg}.parent_intsct = Part.parent_intsct{ireg};
            end
       end
       
       if isfield(Part,'iter_intsct1')
          reg_part{ireg}.iter_intsct1 = Part.iter_intsct1{ireg};       
       end
       if isfield(Part,'parent_intsct1')
          reg_part{ireg}.parent_intsct1 = Part.parent_intsct1{ireg};       
       end
       
      if isfield(Part,'iter_int')
        if iscell(Part.iter_int{ireg})
            reg_part{ireg}.iter_int = cell2mat(Part.iter_int{ireg});
        else
            reg_part{ireg}.iter_int = Part.iter_int{ireg};
        end
      end
      
      if isfield(Part,'iter_int2')
        if iscell(Part.iter_int2{ireg})
            reg_part{ireg}.iter_int2 = cell2mat(Part.iter_int2{ireg});
        else
            reg_part{ireg}.iter_int2 = Part.iter_int2{ireg};
        end
      end
      if isfield(Part,'iter_int3')
        if iscell(Part.iter_int3{ireg})
            reg_part{ireg}.iter_int3 = cell2mat(Part.iter_int3{ireg});
        else
            reg_part{ireg}.iter_int3 = Part.iter_int3{ireg};
        end
      end
      if isfield(Part,'iter_int4')
        if iscell(Part.iter_int4{ireg})
            reg_part{ireg}.iter_int4 = cell2mat(Part.iter_int4{ireg});
        else
            reg_part{ireg}.iter_int4 = Part.iter_int4{ireg};
        end
      end
      if isfield(Part,'iter_int5')
        if iscell(Part.iter_int5{ireg})
            reg_part{ireg}.iter_int5 = cell2mat(Part.iter_int5{ireg});
        else
            reg_part{ireg}.iter_int5 = Part.iter_int5{ireg};
        end
      end
      if isfield(Part,'iter_int6')
        if iscell(Part.iter_int6{ireg})
            reg_part{ireg}.iter_int6 = cell2mat(Part.iter_int6{ireg});
        else
            reg_part{ireg}.iter_int6 = Part.iter_int6{ireg};
        end
      end
      
      if isfield(Part,'activeConstraints_acc')
         reg_part{ireg}.activeConstraints_acc = Part.activeConstraints_acc{ireg};
      end
      
      if isfield(Part,'parent_acc')
         reg_part{ireg}.parent_acc = Part.parent_acc{ireg};
      end
      
      if isfield(Part,'iter_prev')
         reg_part{ireg}.iter_prev = Part.iter_prev{ireg};
      end
      
      if isfield(Part,'theta_cheby')
         reg_part{ireg}.theta_cheby = Part.theta_cheby{ireg};
      end
      
      if isfield(Part,'xi')
         reg_part{ireg}.xi = Part.xi{ireg};
      end
      
    end %for
    
end %function
