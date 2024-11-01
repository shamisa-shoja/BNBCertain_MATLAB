function entiresol = subsol2entiresol_sh(H,f,f_theta,subsol,vartype,ivalues)

if(~isempty(subsol))
 if(~isempty(subsol.Fi))
    Fi = subsol.Fi;
    Gi = subsol.Gi;
    
    entiresol = subsol;
    
    nparam = size(Fi{1},2);
    
    integers = sort(vartype(ivalues ~= -1));
    Fi_entire = cell(1,length(Fi));
    Gi_entire = cell(1,length(Gi));
    
    for i = 1:length(Fi)
        Fi_entire{i} = Fi{i};
        Gi_entire{i} = Gi{i};
        for j = 1:length(integers)
            Fi_entire{i} = [Fi_entire{i}(1:integers(j)-1,:);zeros(1,nparam);Fi_entire{i}(integers(j):end,:)];
            Gi_entire{i} = [Gi_entire{i}(1:integers(j)-1,:);ivalues(j);Gi_entire{i}(integers(j):end,:)];

%             % update Active set if there is: Warm start
%             if isfield(subsol,'AS')    %SH
%                 AS_entire{i} = subsol.AS{i};
%                 %AS_entire{i} = [AS_entire{i}(1:integers(j)-1,:);1;AS_entire{i}(integers(j):end,:)];
%                 m = length(subprob.b); %number of inequality constraints
%                 AS_entire{i} = [AS_entire{i}; m + integers(j)];
%                 entiresol.AS = AS_entire;
%             end
        end
        
        if(entiresol.Ci{i} < inf)
            
            if isempty(f_theta)
                entiresol.Ai{i} = 0.5*Fi_entire{i}'*H*Fi_entire{i};
                entiresol.Bi{i} = (Gi_entire{i}'*H + f')*Fi_entire{i};
                entiresol.Ci{i} = 0.5*Gi_entire{i}'*H*Gi_entire{i} + f'*Gi_entire{i};
            else %Ss
                entiresol.Ai{i} = 0.5*Fi_entire{i}'*H*Fi_entire{i} + f_theta'* Fi_entire{i};   
                entiresol.Bi{i} = (Gi_entire{i}'*H + f')*Fi_entire{i} + Gi_entire{i}'* f_theta;
                entiresol.Ci{i} = 0.5*Gi_entire{i}'*H*Gi_entire{i} + f'*Gi_entire{i};
            end
            
        end
    end
    
    entiresol.Fi = Fi_entire;
    entiresol.Gi = Gi_entire;
else
    entiresol = subsol;
end
end %SH