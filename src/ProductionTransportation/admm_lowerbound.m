function [ f_opt,v,w_k,u_k,t_k,p_k,P_k,tp_k,CPUTime ] = admm_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1,theta,rho_const)

[l, data_size] = size(A);
[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

tic 
cvx_begin sdp

    variable v(m) nonnegative;
    variable w_k(K,n);
    variable u_k(K,m);
    variables p_k(K,m1) t_k(K);
    variable P_k(m1,m1,K) symmetric;

    variable tp_k(K,data_size);
    

    %%
    obj = -sum(v);
    for k = 1:K 
        obj = obj + t_k(k)*beta_k(k);
    end 
    for k = 1:K 
        for j = 1:n
            obj = obj + w_k(k,j)*d(j); 
        end 
    end 

    for k = 1:K 
        obj = obj + theta(k,:)*(tp_k(k,:)'-B*p_k(k,:)') - (rho_const/2)*(tp_k(k,:)'-B*p_k(k,:)')'*(tp_k(k,:)'-B*p_k(k,:)') ;
    end 

    maximize(    obj   ); 
    %%
    subject to

    1 - sum(t_k) == 0;
    
    tep_p = p_k(1,:);
    for k = 2:K 
        tep_p = tep_p + p_k(k,:);
    end 
    sqrt(gamma1) - norm(tep_p) >= 0;

    for k = 1:K 
        t_k(k)*(A*mu-b)' + tp_k(k,:)*A1Complete'*A' <= 0;
    end 

    tep_P = P_k(:,:,1);
    for k = 2:K 
        tep_P = tep_P + P_k(:,:,k);
    end
    eye(m1) - tep_P >=0;

    tep_x = t_k(1)*c + u_k(1,:)';
    for k = 2:K
        tep_x = tep_x + t_k(k)*c + u_k(k,:)';
    end
    tep_x + v >= 0;

    for k = 1:K
        [t_k(k) p_k(k,:); p_k(k,:)' P_k(:,:,k)] >=0;

        expression temp_w(m*n);
        expression temp_u(n*m);
        for i = 1:m 
            temp_w( 1+(i-1)*n:i*n ) = w_k(k,:);
        end 
        for i = 1:m
            temp_u( 1+(i-1)*n:i*n ) = u_k(k,i);
        end 
        alpha_k(k)*t_k(k)*mu' + alpha_k(k)*tp_k(k,:)*A1Complete' - temp_w' - temp_u' >= 0;
    end

 

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