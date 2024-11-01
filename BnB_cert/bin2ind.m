function ind_vec = bin2ind(bin_vec)

[n,m] = size(bin_vec);

if(m > 1 && n > 1)
    error('Matrices not supported');
elseif(n >= m)
    ind_vec = find(bin_vec);
else % n < m
    ind_vec = find(bin_vec')';
end