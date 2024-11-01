% do some settings on mp-miqp problem

function [prob,x0] = probsettings(prob,x0,options)
    
    nx = size(prob.H,1);
    
    if (~isfield(prob,'lb') || isempty(prob.lb) )
        prob.lb=-inf(nx,1);
        prob.ub=inf(nx,1);
    end

    prob.lb(prob.vartype)=zeros(size(prob.vartype));
    prob.ub(prob.vartype)=ones(size(prob.vartype)); 

    if ~issymmetric(prob.H)
            prob.H=(prob.H+prob.H')/2;
    end
    if ~isfield(prob,'f_theta')
       prob.f_theta=[];
    end
    if ~isfield(prob,'Aeq')
        prob.Aeq=[];
    end
    
    if ~isfield(prob,'beq')
        prob.beq=[];
    end

    if size(prob.H,1) ~= size(prob.H,2)
        error('H is not square')
    end
    
    if ~(max(max(abs(prob.H-prob.H'))) <= eps^(2/3)*max(max(abs(prob.H))))
        prob.H=0.5*(prob.H+prob.H');
        if options.verbose >= 1
            warning('H is not symmetric: replaced by its symmetric part')
        end
    end

    if (size(prob.f,1) ~= 1) & (size(prob.f,2) ~= 1)
        %error('f must be a vector')
        warning('f must be a vector')
    end

    if ~isempty(prob.b) & (size(prob.b,1) ~= 1) & (size(prob.b,2) ~= 1)
        %error('b must be a vector')
        warning('b must be a vector')
    end

    prob.f = prob.f(:);   % f and b are column vectors
    prob.b = prob.b(:);

    if isempty(prob.b)
        prob.b=zeros(0,1);
    end
    if isempty(prob.A)
        prob.A=zeros(0,nx);
    end

    if size(prob.A,1) ~= size(prob.b,1)
        error('A and b have incompatible dimensions')
    end

    if size(prob.A,2) ~= nx
        error('A and H have incompatible dimensions')
    end

    if size(prob.Aeq,1) ~= size(prob.beq,1)
        error('Aeq and beq have incompatible dimensions')
    end

    if (size(prob.Aeq,2) ~= nx)&(~isempty(prob.Aeq))
        error('Aeq and H have incompatible dimensions')
    end

    prob.lb      = prob.lb(:);
    prob.ub      = prob.ub(:);
    prob.vartype = prob.vartype(:);
    
    
    if ischar(prob.vartype)
        if length(prob.vartype) ~= nx
            error('wrong dimension of vartype as character array')
        end
        Ccons = find(prob.vartype=='C');
        Bcons = find(prob.vartype=='B');
        Icons = find(prob.vartype=='I');
        if length([Ccons(:); Bcons(:); Icons(:)]) ~=  nx;
            error('wrong entries in vartype')
        end
        if ~isempty(Icons)
            error('specifications for integer variables are not supported')
        end
        % deletes the variable in character syntax and uses the syntax with vector
        % of indices instead
        prob.vartype = Bcons(:);
    end

    if size(prob.lb,1)~=nx & ~isempty(prob.lb)
        error('lb has wrong dimensions')
    end
    if size(prob.ub,1)~=nx & ~isempty(prob.ub)
        error('ub has wrong dimensions')
    end

    if max(prob.vartype) > nx
        error('largest index in vartype is out of range')
    end
    if min(prob.vartype) < 1
        error('smallest index in vartype is out of range')
    end
    if find(diff(sort(prob.vartype)) == 0)
        error('binary variables are multiply defined')
    end
    if length(prob.vartype) > nx
        error('too many entries in vartype')
    end
    if floor(prob.vartype) ~= prob.vartype
        error('fractional number in vartype not allowed')
    end

    % Define default values for lb,ub,x0
    % ------------------------------------
    if isempty(prob.lb)
        prob.lb =-inf*ones(nx,1);
    end
    if isempty(prob.ub)
        prob.ub = inf*ones(nx,1);
    end
    
    if find(prob.ub-prob.lb < 0)
        error('infeasible constraints specified in ub and lb')
    end
    
    % settings on x0
    if isempty(x0)
        x0  = zeros(nx,1);
    end
    if ~isstruct(x0)
        if size(x0,1)~=nx & ~isempty(x0)
            error('x0 has wrong dimensions')
        end
        x0 = x0(:);
    end
    
    % settings on vartype
    % checking whether the bounds 0,1 on the binary variables are already present in
    % the problem constraints, if not, add them
    aux1         = prob.lb(prob.vartype);
    index1       = find(aux1<0);
    aux1(index1) = 0;
    prob.lb(prob.vartype)  = aux1;

    aux2         = prob.ub(prob.vartype);
    index2       = find(aux2>1);
    aux2(index2) = 1;
    prob.ub(prob.vartype)  = aux2;
    
    % define continuous and binary constraints
    prob.cont = (1:nx)';
    prob.cont(prob.vartype) = [];       % Indices of continuous variables
    prob.bin  = (1:nx)';
    prob.bin(~prob.vartype) = [];       % Indices of continuous variables
    
end %main function