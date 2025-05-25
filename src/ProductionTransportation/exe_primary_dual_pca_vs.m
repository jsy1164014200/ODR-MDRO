clear;
clc;

record = zeros (10,100);    
Number = 3;  

for iterationNumber = 1 : Number 
    load("./Data/15-4-25-"+iterationNumber) % TODO

    m1_list = [1*data_size, 0.8*data_size, 0.6*data_size, 0.4*data_size, 0.2*data_size, K];
    P_list = [2, 4, 5];
    
   
    % Mosek 9.3 solver
    if data_size <= 250
        [ f_opt1,~,~,~,~,~,CPUTime1 ] = primary_sdp( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n);
        record(iterationNumber,1) = f_opt1;
        record(iterationNumber,2) = CPUTime1;
%         [ dual_opt,~,~,~,~,~,~,dual_CPUTime ] = dual_sdp( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n);
%         record(iterationNumber,3) = dual_opt;
%         record(iterationNumber,4) = dual_CPUTime;
    end
    f_opt1 = 1;

     


    % Burer Algorithm
    % [ f_opt,X_opt,CPUTime ] = BurerAlgorithm( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g);
    % record(iterationNumber,3) = f_opt;
    % record(iterationNumber,4) = CPUTime;
    

    % [ f_opt,X_opt,CPUTime ] = BurerAlgorithmForcerank( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g);
    % record(iterationNumber,5) = f_opt;
    % record(iterationNumber,6) = CPUTime;
    
    % record(iterationNumber, 7) = abs(f_opt1 - f_opt) / abs(f_opt1) * 100;

%     for P_index = 1:3
%         P = P_list(P_index); 
%         if data_size/P > 320
%             continue;
%         end 
%         [ VS_opt,~,~,VS_Time ] = vs_upperbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,P);
%         record(iterationNumber,2*(P_index-1)+5) = VS_opt;
%         record(iterationNumber,2*(P_index-1)+6) = VS_Time;
%     end 

    for m1_index = 2 : 6
        m1 = m1_list(m1_index); 
        m1 = round(m1);
        if m1 > 250
            continue;
        end

        B = zeros(data_size,m1);
        for i = 1 : m1
            B(i,i) = 1;  
        end

        [PCA_opt,~,~,~,~,~,~,PCA_time] = pca_lowerbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
        record(iterationNumber,4*(m1_index-1)+11) = PCA_opt;
        record(iterationNumber,4*(m1_index-1)+12) = PCA_time;

        [PCA_opt,~,~,~,~,~,~,PCA_time] = pca_upperbound( ZIGMA,mu,A,b,gamma1,gamma2,c,d,alpha_k,beta_k,K,m,n,B,m1);
        record(iterationNumber,4*(m1_index-1)+13) = PCA_opt;
        record(iterationNumber,4*(m1_index-1)+14) = PCA_time;
        
    end

    
 


end