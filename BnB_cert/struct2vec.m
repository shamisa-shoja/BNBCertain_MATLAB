   % find key in a struct and pass it as a vector or matrix than struct
   
function [vec_out,max_vec_out] = struct2vec(part,in)
   
    vec_out = zeros(length(part),1);

    for k =1 : length(part)
        vec_out(k) = part{k}.(in); 
    end

    [max_vec_out,max_vec_out_ind] = max(vec_out); %upper bound of in
    [min_vec_out,min_vec_out_ind] = min(vec_out); %lower bound of in
end