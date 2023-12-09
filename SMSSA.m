function [x_c,J] = SMSSA(X, n_comp, W,L, delta ,eps) 
if ~exist('delta', 'var')
 delta = 1; 
end

if ~exist('eps', 'var')
    eps = 0.03;
end


[N,n_c] = size(X);
n_co = max(n_comp); % using max of all comps as the n_comp
x_hat = zeros(N,n_co,n_c);

for p = 1:delta: N-W+1
    
    nc = 1 + (W-1)/2;
    m = p+W -1;
    
    x_tilde = mssa_decomp( X(p:m,:), L, n_comp, eps);

    if p ==1
        x_hat(1:nc,:,:) = x_tilde(1:nc,:,:);
        
        
    else
        
        for f = 1:n_c
            
            J = [];

            for j = 1:n_co

                x_tilde_j = x_tilde(1:nc,j,f);
              
                x_hat_k = x_hat(p-delta:p-delta+nc -1,:,f);
                
                for k = 1:n_co
                    dists(k) = norm(x_tilde_j - x_hat_k(:,k));
                end
                [~, id] = sort(dists);
                j_ = id(1);
                

                e = 2;
                while ismember(j_,J) && e<=n_co
    %                 disp('in')
                    j_ = id(e);
                    e = e+1;
                end
    %             disp(j_)
                J = [J,j_];
                
    if p < N - W+ 1
                    x_hat(p-1+nc:p-1+nc+delta-1, j_,f) = x_tilde(nc:nc+delta-1,j,f);


                else

                    x_hat( p-1+nc: p-1+W, j_,f) = x_tilde(nc:W,j,f);


    end
                



            end
        
        end
        
        
    end
    
end
x_c = x_hat;
end