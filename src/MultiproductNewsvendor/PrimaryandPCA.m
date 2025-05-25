clear;
clc;



record = zeros (10,100);    
Number = 5;  
ID = 1;

for iterationNumber = 1 : Number 
    load("./Data/2000-"+iterationNumber) % TODO
    m1_list = [1*m, 0.8*m, 0.6*m, 0.4*m, 0.2*m, 2];
    %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
    
    % Mosek 9.3 solver
    if m <= 320
        [optimalB, f_opt1,X_opt1,CPUTime1 ] = PrimaryNewsvendor( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g);
        record(ID,1) = f_opt1;
        record(ID,2) = CPUTime1;
    end
    f_opt1 = 1;
%     [PCA_opt,~,PCA_time] = PrimaryNewsvendorUpperGivenB(ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,2, optimalB);
%     [~,~,~,~,~,~,PCA_opt,PCA_time] = PrimaryNewsvendorGivenB(ZIGMA,mu,2,A,b,gamma1,gamma2,c,v,g, optimalB);
%     record(ID,3) = PCA_opt;
%     record(ID,4) = PCA_time;
%     return;
    
    % Burer Algorithm
    [ f_opt,X_opt,CPUTime ] = BurerAlgorithm( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g);
    record(ID,3) = f_opt;
    record(ID,4) = CPUTime;
    
%    [ f_opt,X_opt,CPUTime ] = BurerAlgorithmForcerank( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g);
%    record(ID,5) = f_opt;
%    record(ID,6) = CPUTime;
    
    record(ID, 7) = abs(f_opt1 - f_opt) / abs(f_opt1) * 100;

    

    
    

    for m1_index = 1 : 6
        m1 = m1_list(m1_index); 
        m1 = round(m1);
        if m1 > 320
            continue;
        end

        B = zeros(m,m1);
        for i = 1 : m1
            B(i,i) = 1;  
        end

        [~,~,~,~,~,~,PCA_opt,PCA_time] = PrimaryNewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B);
        record(ID,4*(m1_index-1)+11) = PCA_opt;
        record(ID,4*(m1_index-1)+12) = PCA_time;

        [PCA_opt,~,PCA_time] = PrimaryNewsvendorUpperGivenB(ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,m1, B);
        record(ID,4*(m1_index-1)+13) = PCA_opt;
        record(ID,4*(m1_index-1)+14) = PCA_time;
        
    end

            
    ID = ID + 1;
end




