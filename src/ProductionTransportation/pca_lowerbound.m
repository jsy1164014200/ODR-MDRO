function [ f_opt,x,s,z_k,q,Q,lambda_k,CPUTime ] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1)

[l, data_size] = size(A);
[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
A1Complete = A1Complete*B;

tic 
cvx_begin sdp

    variable x(m) nonnegative;
    variable z_k(K,data_size);
    variables q(m1) s;
    variable Q(m1,m1) symmetric;
    variable lambda_k(K, l);

    %%
    minimize(    s + sum(sum(  (gamma2*eye(m1)).*Q  )) + sqrt(gamma1)*norm(q,2)    );
    %%
    subject to

    vec(z_k) >= 0;
    vec(lambda_k) >= 0;

    for k = 1:K
        [ s-c'*x-beta_k(k)-lambda_k(k,:)*(b-A*mu) - alpha_k(k)*z_k(k,:)*mu     1/2*( q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') )' ; 1/2*( q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') ) Q] >= 0 ;

        z_tep = reshape(z_k(k,:),n,m)';
        for j = 1 : n                
            sum (z_tep(:,j)) == d(j);
        end
        for i = 1 : m          
            sum (z_tep(i,:)) == x(i);
        end

    end

    x <= 1;

cvx_end
CPUTime = toc;

disp(['Problem is ' cvx_status])
if ~strfind(cvx_status,'Solved')
    return
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['optimal value of cvx:',num2str(cvx_optval)]);

f_opt=cvx_optval;
 

end