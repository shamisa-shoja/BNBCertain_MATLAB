
function th_sample = sampling(N_sample,lb_th,ub_th,sample_kind)
    
    N_theta = size(lb_th,1);    
    th_sample = []; 
    
    switch sample_kind
        
        % random sapmles in parameter set --> Sparse points
    case 'rand'         
        delta_th=ub_th - lb_th;
        th_sample = zeros(N_theta, N_sample);
        
        for i=1: N_sample
            %th_sample(:,i)=lb+rand*delta_th;   % generate 3 similar points
            for j =1 : N_theta
                th_sample(j,i)=lb_th(j)+ rand*delta_th(j);
            end
        end
        
        % deterministic sapmles in parameter set --> Dense points
    case 'deter'

        vect = zeros(N_theta, N_sample);
        for j =1 : N_theta
            vect(j,:) = linspace(lb_th(j),ub_th(j),N_sample);
        end
        %th_sample2 = meshgrid(vec(1,:),vec(2,:));
        %th_sample = zeros(N_theta, N_sample^N_theta); --> dimension of theta        
%        block = cell(0);
%         if N_theta == 2 % 2dimension
%             for j = 1: N_sample
%                 block{j} = [vec(1,:); vec(2,j)*ones(1,N_sample)];
%                 th_sample = [th_sample, block{j}];
%     %                 k = j + N_sample;
%     %                 th_sample(1,j:k)= vec(1,:);
%     %                 th_sample(2,j:k)= vec(2,j)*ones(1,N_sample);
%             end           
%         elseif N_theta == 3 % 3dimension
%             th_sample3 = []; block3 = cell(0);
%             for j = 1: N_sample
%                 block3{j} = [th_sample, vec(3,j)*ones(1,N_sample)];
%                 th_sample3 = [th_sample3, block3{j}];
%             end
%             th_sample = th_sample3;              
%         end

%         if N_theta == 2 % 2D
%              th_sample = meshgrid(vec(1,:),vec(2,:));  
%         elseif N_theta == 3 % 3D
%              th_sample = meshgrid(vec(1,:),vec(2,:),vec(3,:)); 
%         elseif N_theta == 4 % 4D
%              th_sample = meshgrid(vec(1,:),vec(2,:),vec(3,:),vec(4,:));     
%         end
%         %th_sample2 = meshgrid(vec);
%         %param = 1: N_theta;
%         %%th_sample2 = meshgrid(vec(param(1),:),vec(param(2),:),vec(param(N_theta),:));
%         %th_sample2 = meshgrid(vec(param,:));
%         
%         [X,Y]= meshgrid(vec(1,:),vec(2,:));
%         figure; plot(X,Y,'+')
%         theta = zeros(N_theta-1,N_sample^2);

        theta = [];
        if N_theta == 2 % 2D
             %th_sample = meshgrid(vec(1,:),vec(2,:));
             for ii = 1: N_sample
                for jj = 1: N_sample
                    %theta_col = [X(ii),Y(jj)];
                    theta_col = [vect(1,ii),vect(2,jj)];
                    theta = [theta, theta_col'];
                 end % for jj
             end % for ii
             
        elseif N_theta == 3 % 3D
             %th_sample = meshgrid(vec(1,:),vec(2,:),vec(3,:)); 
             for ii = 1: N_sample
                for jj = 1: N_sample
                    for kk = 1: N_sample
                        theta_col = [vect(1,ii),vect(2,jj),vect(3,kk)];
                        theta = [theta, theta_col'];
                    end % for kk
                 end % for jj
             end % for ii
             
        elseif N_theta == 4 % 4D
             %th_sample = meshgrid(vec(1,:),vec(2,:),vec(3,:),vec(4,:));  
             for ii = 1: N_sample
                for jj = 1: N_sample
                    for kk = 1: N_sample
                        for ll = 1: N_sample
                            theta_col = [vect(1,ii),vect(2,jj),vect(3,kk),vect(4,ll)];
                            theta = [theta, theta_col'];
                        end % for ll
                    end % for kk
                 end % for jj
             end % for ii             
        end
        
        th_sample = theta;
        
    end %switch
    
end % function