% Examples of multi-parametric Mixed Integere Linear Prigramming
% Reference: An Algorithm for the Solution of Multiparametric Mixed Integer Linear Programming Problems 
%   J =  min(x) 1/2 * f'*x
%   s.t. Ax <= b+ W*theta
%        x_i = {0,1}  i \in B

function [prob, P_theta, AS0, x0, vartype] = example_MILP(v,n, nb, nth, m, bnd)
    
    switch(v)

        case(1)
            % Cost fubction:
            f = [-3; -8; 4; 2];

            % Constraints
            A = [1 1 0 0; 5 -4 0 0; -8 22 0 0; -4 -1 0 0; 1 0 -10 0; ...
                 0 1 0 -15; -1 0 0 0; 0 -1 0 0];
            b = [13; 20; 121; -8; 0; 0; 0; 0];
            W = [1 0; 0 0; 0 1;0 0 ; 0 0; 0 0; 0 0; 0 0];

            %equality constraints
            A_eq = []; b_eq = []; 
            
            N_x = size(f,1);
            [N_cons,N_theta] = size(W);

            % Parameter set
            lb_theta = 0*ones(N_theta,1);
            ub_theta = 10*ones(N_theta,1);

            vartype = [3; 4]; %index of binary variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       case(2)
           % Cost fubction:
            f = [-3; -2; 10; 5];

            % Constraints
            A = [1 0 0 0; 0 1 0 0; 1 1 0 0; 1 2 0 0;1 0 -20 0; 0 1 0 -20; ...
                 1 -1 0 0; 0 0 -1 -1;-1 0 0 0;0 -1 0 0];
            b = [10; 10; 20; 12; 0; 0; -4; -1; 0; 0];
            W = [1 2 0; -1 1 0; 0 -1 0;1 0 -1; 0 0 0; 0 0 0; 0 0 1; 0 0 0;  0 0 0; 0 0 0];

            %equality constraints
            A_eq = []; b_eq = []; 
            
            N_x = size(f,1);
            [N_cons,N_theta] = size(W);

            % Parameter set
            lb_theta = 0*ones(N_theta,1);
            ub_theta = 5*ones(N_theta,1);

            P_theta = Polyhedron('lb', lb_theta, 'ub', ub_theta); 
            P_theta.minHRep();

            vartype = [3; 4]; %index of binary variables 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case(3)
            % Third example from the reference
            % Cost fubction:
            Nc = 23; Nb = 8; N_x = Nc+ Nb;
            f = zeros(N_x,1);
            f(19:23) = 1;

            % Constraints
            N_cons = 10;
            N_theta = 4;
            A = zeros(N_cons,N_x); b = zeros(N_cons,1);  W = zeros(N_cons,N_theta); 
            A(1:Nc,:) = - eye(Nc);
            %A(1,:) = - eye(Nc);
            
            
            %equality constraints
            N_con_eq = 14;
            A_eq = zeros(N_con_eq,N_x); b_eq = zeros(N_con_eq,1); 
            A_eq(1,1) = 1;   A_eq(1,3:7) = -1;
            A_eq(2,11) = 1;  A_eq(2,6:10) = -1;
            A_eq(3,12) = 1;  A_eq(3,2:5) = -1;  A_eq(12,8:10) = 1;
            A_eq(4,13) = 1;  A_eq(4,3) = -42.6; A_eq(13,6) = -121;
            A_eq(5,14) = 1;  A_eq(5,4) = -42.6; A_eq(14,7) = -121;
            A_eq(6,15) = 1;  A_eq(6,8) = -78.4; 
            A_eq(7,16) = 1;  A_eq(7,9) = -78.4;
            A_eq(8,17) = 1;  A_eq(8,13:14) = -1;
            A_eq(9,18) = 1;  A_eq(9,15:16) = -1;
            A_eq(10,19) = 1;  A_eq(10,1) = -90000; A_eq(19,24) = -9600;
            A_eq(11,20) = 1;  A_eq(11,1) = -40000; A_eq(20,25) = -8500;
            A_eq(12,21) = 1;  A_eq(12,26:27) = -45000; A_eq(21,15) = -25;
            A_eq(13,22) = 1;  A_eq(13,28:29) = -25000; A_eq(22,18) = -14.5;
            A_eq(14,23) = 1;  A_eq(14,20:31) = -20000; 
            

            %N_x = size(f,1);
            %[N_con,N_theta] = size(W);

            % Parameter set
            lb_theta = 0*ones(N_theta,1);
            ub_theta = 1*ones(N_theta,1);

            P_theta = Polyhedron('lb', lb_theta, 'ub', ub_theta); 
            P_theta.minHRep();

            vartype = [Nc+1; Nx]; %index of binary variables 
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case(4)
            % Random example
            N_x = n; Nb = nb; N_c = N_x - Nb;
            N_cons = m; %N_x + 2; %floor(N_x/3); %in simplex method: m>n %
            %N_theta = 4;
            N_theta = nth;
            
            % Cost fubction:
            f = randn(N_x,1); %10*

            % Constraints
            A = randn(N_cons,N_x);%2*
            b = 2+rand(N_cons,1); %5 + 
            W = randn(N_cons,N_theta); %2*

            %equality constraints
            A_eq = []; b_eq = []; 
            
            %N_x = size(f,1);
            %[N_con,N_theta] = size(W);

            % Parameter set
            bnd = 2;
            lb_th = -bnd;  ub_th = bnd;
            lb_theta = lb_th*ones(N_theta,1);
            ub_theta = ub_th*ones(N_theta,1);

            vartype = [N_c+1: N_x]; %index of binary variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case(5)
            % Example 6.3 from reference book: Hybrid MPC changed to be
            % MILP
            % Cost fubction:
            f = [1; 1; 0; 0];

            % Constraints
            A = [-1 0 0 0; -1 0 -1 0; -1 0 0 0; -1 0 1 0; 0 -1 -1 0; ...
                 0 -1 -1 -1; 0 -1 1 0; 0 -1 1 1];
            b = [0; 0; 0; 0; 0; 0; 0; 0];
            W = [1 1; 0 1; -1 1;0 -1 ; 1 2; 0 1; -1 -2; 0 -1];

            %equality constraints
            A_eq = []; b_eq = []; 
            
            N_x = size(f,1);
            [N_cons,N_theta] = size(W);

            % Parameter set
            lb_theta = -2.5 *ones(N_theta,1);
            ub_theta = 2.5 *ones(N_theta,1);

            vartype = [3; 4]; %index of binary variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
        case(6)
            % Random MIQP example
            Nc = 2; Nb = 2; N_x = Nc+ Nb;
            N_cons = N_x; 
            N_theta = 2;
            
            % Cost fubction i QP case
            H1 = 2*randn(N_x); 
            H = (H1*H1')/2; 
            f_theta= randn(N_x,N_theta);
            prob.H = H;
            prob.f_theta = f_theta;
            %
            f = randn(N_x,1); %10*

            % Constraints
            A = rand(N_cons,N_x); 
            b = rand(N_cons,1) + ones(N_cons,1);
            W = 2* rand(N_cons,N_theta); 

            %equality constraints
            A_eq = []; b_eq = []; 
            
            %N_x = size(f,1);
            %[N_con,N_theta] = size(W);

            % Parameter set
            bnd = .5;
            lb_th = -bnd;  ub_th = bnd;
            lb_theta = lb_th*ones(N_theta,1);
            ub_theta = ub_th*ones(N_theta,1);

            vartype = [Nc+1: N_x]; %index of binary variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end     
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%%%%%%%%%%%%%%%%%%%
    % Outputs:
    prob.f = f;
    prob.A = A;
    prob.b = b;
    prob.W = W;
    prob.A_eq = A_eq;
    prob.b_eq = b_eq;
    prob.vartype = vartype;
    prob.N_x = N_x; prob.N_theta = N_theta;
    prob.N_cons = N_cons; 
    
    P_theta = Polyhedron('lb', lb_theta, 'ub', ub_theta); 
    P_theta.minHRep();
        
    AS0 = []; %cold-start; % false *ones(N_x,1);
    x0.F = zeros(N_x,N_theta);
    x0.G = zeros(N_x,1);
    
    %---------------------------------------------
    if (~isfield(prob,'lb') || isempty(prob.lb) )
        prob.lb=-inf(N_x,1);
        prob.ub=inf(N_x,1);
    end
    %---------------------------------------------     
    if ~isfield(prob,'upper_ind')
        prob.box_ind = [1:N_x]; %or nth?
        prob.upper_ind = [1:N_x/2];
        prob.lower_ind = [N_x/2 +1 : N_x];
    end
end