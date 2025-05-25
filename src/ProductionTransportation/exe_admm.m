clear;
clc;

record = zeros(10,100);    
Number = 3;  

for iterationNumber = 1 : Number
    load("./Data/15-8-25-"+iterationNumber) % TODO
    m1=17; % TODO

    B = zeros(data_size,m1);
    for i = 1 : m1
        B(i,i) = 1;  
    end
    
    theta = zeros(K, data_size);
    rho_const = 100;

    totalTime = 0;
    lastOptimalValue = 9999;
    for iterNum = 1:500
        [ f_opt2,v,w_k,u_k,t_k,p_k,P_k,tp_k,CPUTime2 ] = admm_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1,theta,rho_const);
        % record(iterationNumber,10 + 6*(iterNum-1) + 1) = f_opt2;
        % record(iterationNumber,10 + 6*(iterNum-1) + 2) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        
        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
            % record(iterationNumber,10 + 6*(iterNum-1) + 5) = iterNum;
            % record(iterationNumber,10 + 6*(iterNum-1) + 6) = totalTime;
            record(iterationNumber,1) = iterNum;
            record(iterationNumber,2) = f_opt2;
            record(iterationNumber,3) = totalTime;
            for k = 1:K
                disp( (tp_k(k,:)'-B*p_k(k,:)')'*(tp_k(k,:)'-B*p_k(k,:)') );
                disp(theta(k,:));
            end
            break;
        else
            lastOptimalValue = f_opt2;
        end
        SofB = -theta(1,:)'*p_k(1,:) + rho_const*tp_k(1,:)'*p_k(1,:);
        for k = 2:K 
            SofB = SofB - theta(k,:)'*p_k(k,:) + rho_const*tp_k(k,:)'*p_k(k,:); 
        end  
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';

        for k = 1:K 
            theta(k,:) = theta(k,:) - rho_const*(tp_k(k,:) - p_k(k,:)*B');
        end  
    end

    [lowerbound,~,~,~,~,~,~,lb_time] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
    record(iterationNumber,4) = lowerbound;

    [ f_opt,x,s,z_k,q,Q,lambda_k,CPUTime ] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
    [ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
    S = zeros(1, K);
    M = zeros(data_size,data_size,K);
    for k = 1:K 
        S(k) = s-c'*x-beta_k(k)-lambda_k(k,:)*(b-A*mu) - alpha_k(k)*z_k(k,:)*mu; 
        M(:,:,k) = ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') ) * ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') )'; 
    end 
    
    
    real_S = max(S); 
    P = 0;
    for k = 1:K 
        P = P + (gamma2/4) * trace( eye(data_size)*M(:,:,k) );
    end 
    
    if sqrt(P) < real_S
        record(iterationNumber,5) = f_opt + P/real_S;
    else
        record(iterationNumber,5) = f_opt + 2*sqrt(P) - real_S;
    end 



    %%% Upper Bound
    B = zeros(data_size,m1);
    for i = 1 : m1
        B(i,i) = 1;  
    end
    
    theta = zeros(K, data_size);
    rho_const = 100;

    totalTime = 0;
    lastOptimalValue = 9999;
    for iterNum = 1:500
        [ f_opt2,x,z_k,q,Q,lambda_k,u_k,ubar_k,CPUTime2 ] = admm_upperbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1,theta,rho_const);
        % record(iterationNumber,2*(iterNum-1) + 41) = f_opt2;
        % record(iterationNumber,2*(iterNum-1) + 42) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
            record(iterationNumber,6) = iterNum;
            record(iterationNumber,7) = f_opt2;
            record(iterationNumber,8) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = theta(1,:)'*u_k(1,:) + rho_const*ubar_k(1,:)'*u_k(1,:);
        for k = 2:K 
            SofB = SofB + theta(k,:)'*u_k(k,:) + rho_const*ubar_k(k,:)'*u_k(k,:); 
        end  
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';

        for k = 1:K 
            theta(k,:) = theta(k,:) + rho_const*(ubar_k(k,:) - u_k(k,:)*B');
        end  
    end
    [upperbound,~,~,~,~,~,~,lb_time] = pca_upperbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
    record(iterationNumber,9) = upperbound; 

    [ f_opt,x,s,z_k,q,Q,lambda_k,CPUTime ] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
    [ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
    S = zeros(1, K);
    M = zeros(data_size,data_size,K);
    for k = 1:K 
        S(k) = s-c'*x-beta_k(k)-lambda_k(k,:)*(b-A*mu) - alpha_k(k)*z_k(k,:)*mu; 
        M(:,:,k) = ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') ) * ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') )'; 
    end 
    
    
    real_S = max(S); 
    P = 0;
    for k = 1:K 
        P = P + (gamma2/4) * trace( eye(data_size)*M(:,:,k) );
    end 
    
    if sqrt(P) < real_S
        record(iterationNumber,10) = f_opt + P/real_S;
    else
        record(iterationNumber,10) = f_opt + 2*sqrt(P) - real_S;
    end


    %%% Revisted Lower Bound
    if m1 > K
        continue;
    end
    B = zeros(data_size,K); % B = [Bpri, Bbar]
    for i = 1 : K
        B(i,i) = 1;  
    end
    
    theta = zeros(K, data_size);
    rho_const = 100;

    totalTime = 0;
    lastOptimalValue = 9999;
    for iterNum = 1:500
        [ f_opt2,x,z_k,q,Q,lambda_k,u_k,ubar_k,CPUTime2 ] = admm_revisted_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1,theta,rho_const);
        % record(iterationNumber,2*(iterNum-1) + 41) = f_opt2;
        % record(iterationNumber,2*(iterNum-1) + 42) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
            record(iterationNumber,11) = iterNum;
            record(iterationNumber,12) = f_opt2;
            record(iterationNumber,13) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = theta(1,:)'*u_k(1,:) + rho_const*ubar_k(1,:)'*u_k(1,:);
        for k = 2:K 
            SofB = SofB + theta(k,:)'*u_k(k,:) + rho_const*ubar_k(k,:)'*u_k(k,:); 
        end  
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';

        for k = 1:K 
            theta(k,:) = theta(k,:) + rho_const*(ubar_k(k,:) - u_k(k,:)*B');
        end  
    end

    [lowerbound,~,~,~,~,~,~,lb_time] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,K); % consider the entire B (data_size x K)
    record(iterationNumber,14) = lowerbound;

    [lowerbound,~,~,~,~,~,~,lb_time] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B(:,1:m1),m1); % consider the partial B (data_size x m1)
    record(iterationNumber,15) = lowerbound;
    

    [ f_opt,x,s,z_k,q,Q,lambda_k,CPUTime ] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,K);
    [ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
    S = zeros(1, K);
    M = zeros(data_size,data_size,K);
    for k = 1:K 
        S(k) = s-c'*x-beta_k(k)-lambda_k(k,:)*(b-A*mu) - alpha_k(k)*z_k(k,:)*mu; 
        M(:,:,k) = ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') ) * ( B*q+ A1Complete'*(A'*lambda_k(k,:)'-alpha_k(k)*z_k(k,:)') )'; 
    end 
    
    
    real_S = max(S); 
    P = 0;
    for k = 1:K 
        P = P + (gamma2/4) * trace( eye(data_size)*M(:,:,k) );
    end 
    
    if sqrt(P) < real_S
        record(iterationNumber,16) = f_opt + P/real_S;
    else
        record(iterationNumber,16) = f_opt + 2*sqrt(P) - real_S;
    end
 
 
end