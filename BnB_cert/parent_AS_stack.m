% stack parent and active set info

function [parent_acc,AS_acc,AS_log_acc] = parent_AS_stack(old_parent,new_parent,old_AS,new_AS,old_AS_log,new_AS_log)

        if ~isempty(new_AS)
%             [nr_activ_old,mii] = size(old_AS);
%             [nr_activ_new,mjj] = size(new_AS);
%             diff_activ = nr_activ_old - nr_activ_new;
% 
%             if (diff_activ >= 0)
%                 AS_acc = [old_AS -5*ones(nr_activ_old,1) [new_AS; -ones(diff_activ,mjj)] ];
%             else
%                 AS_acc = [ [old_AS ; -ones(-diff_activ,mii)] -5*ones(nr_activ_new,1) new_AS ];                    
%             end
            AS_acc = [old_AS, new_AS];
        else
            AS_acc = old_AS;
        end
        %__________________________________________________________________
        
            if (~isempty(old_parent) && iscell(old_parent))
                old_parent =cell2mat(old_parent);
            end

            if  ~isempty(old_parent) %cell2mat(old_parent)) 
                if ~isempty(new_parent) 

                    [ni,mi] = size(old_parent);
                    [nj,mj] = size(new_parent);
                    di = ni-nj;
                    if (di >= 0 )
                        % TODO: change logical to number 
                        parent_acc = [old_parent -5*ones(ni,1) [new_parent; -ones(di,mj)] ];
                    else
                        di=abs(di);
                        parent_acc = [[old_parent; -ones(di,mi)] -5*ones(nj,1) new_parent];
                    end 
                else 
                   parent_acc = old_parent;
                end
            else 
                parent_acc = new_parent;
            end
        %__________________________________________________________________
        if ~isempty(new_AS_log)
            [nr_activ_old_log,mii] = size(old_AS_log);
            [nr_activ_new_log,mjj] = size(new_AS_log);
            diff_activ_log = nr_activ_old_log - nr_activ_new_log;

            if (diff_activ_log >= 0)
                AS_log_acc = [old_AS_log -5*ones(nr_activ_old_log,1) [new_AS_log; -ones(diff_activ_log,mjj)] ];
            else
                AS_log_acc = [ [old_AS_log ; -ones(-diff_activ_log,mii)] -5*ones(nr_activ_new_log,1) new_AS_log ];                    
            end
%            AS_acc = [old_AS, new_AS];
        else
            AS_log_acc = old_AS_log;
        end

end %function