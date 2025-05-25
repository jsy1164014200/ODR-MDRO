function [ f_opt,x,z_k,q,Q,lambda_k,u_k,ubar_k,CPUTime ] = admm_revisted_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1,theta,rho_const)

[l, data_size] = size(A);
[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);


tic 
cvx_begin sdp

    variable x(m) nonnegative;
    variable z_k(K,data_size);
    variables q(data_size) s;
    variable Q(m1,m1) symmetric;
    variable lambda_k(K, l);
    variable u_k(K, K);
    variable ubar_k(K, data_size);

     
    obj = theta(1,:)*( ubar_k(1,:)' - B*u_k(1, :)' ) + (rho_const/2)*( ubar_k(1,:)' - B*u_k(1, :)' )' * ( ubar_k(1,:)' - B*u_k(1, :)' );
    for k = 2:K
        obj = obj + theta(k,:)*( ubar_k(k,:)' - B*u_k(k, :)' ) + (rho_const/2)*( ubar_k(k,:)' - B*u_k(k, :)' )' * ( ubar_k(k,:)' - B*u_k(k, :)' );
    end 

    %%
    minimize(    s + sum(sum(  (gamma2*eye(m1)).*Q  )) + sqrt(gamma1)*norm(q,2)  + obj  );    
    %%
    subject to

    vec(z_k) >= 0;
    vec(lambda_k) >= 0;

    for k = 1:K
        [ s-c'*x-beta_k(k)-lambda_k(k,:)*(b-A*mu) - alpha_k(k)*z_k(k,:)*mu     1/2*u_k(k, 1:m1) ; 1/2*u_k(k, 1:m1)'  Q] >= 0 ;
        q + A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') == ubar_k(k,:)';

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
x_opt=x;
z_opt=z_k;

end