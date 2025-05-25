clear;
clc;


record = zeros (10,100);    
Number = 5; 
ID = 1;

for iterationNumber = 1 : Number
    load("./Data/2000-"+iterationNumber) % TODO 
    %%%%%%%%% Solving original and LB problems and saving their results %%%%%%
    
    f_opt1 = 1;
    
    %LowerBound plan A
    m1=2;
    % init B 
    B = zeros(m,m1);
    % randomRows = randperm(m);
    for i = 1 : m1
        B(i,i) = 1; %i  or m-m1+i or randomRows(i)
    end
    
    theta = {};
    theta{1} = zeros(m, 1);
    theta{2} = zeros(m, 1);

    rho_const = 100;

    totalTime = 0;
    lastOptimalValue = f_opt1;
    for iterNum = 1:500
        [ p,tp,t,P,f_opt2,CPUTime2] = ADMMnewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B, theta, rho_const);
%                 record(ID,10 + 6*(iterNum-1) + 1) = f_opt2;
%                 record(ID,10 + 6*(iterNum-1) + 2) = CPUTime2;
        totalTime = totalTime + CPUTime2; 


       
        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
%                     record(ID,10 + 6*(iterNum-1) + 5) = iterNum;
%                     record(ID,10 + 6*(iterNum-1) + 6) = totalTime;
            record(ID,1) = iterNum;
            record(ID,2) = f_opt2;
            record(ID,3) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = -( theta{1}*p{1}' + theta{2}*p{2}' - rho_const*tp{1}*p{1}' - rho_const*tp{2}*p{2}' )  ;
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';


        theta{1} = theta{1} - rho_const*(tp{1} - B*p{1});
        theta{2} = theta{2} - rho_const*(tp{2} - B*p{2});
    end

    [X,s,Landa1, Landa2, qrr, Qr,real_opt,~] = PrimaryNewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B);
    record(ID,4) = real_opt;





    % Upper Bound
    m1=2;
    B = zeros(m,m1);
    for i = 1 : m1
        B(i,i) = 1; 
    end
    
    beta1 = zeros(m, 1);
    beta2 = zeros(m, 1);
    rho_const = 200;

    totalTime = 0;
    lastOptimalValue = f_opt1;
    for iterNum = 1:500
        [X,s,Landa1,Landa2,q,Qr,u1,u2,ubar1,ubar2, f_opt2, CPUTime2 ] = upperboundADMMGivenB(ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,m1, B, beta1, beta2, rho_const);
%         record(ID,2*(iterNum-1) + 41) = f_opt2;
%         record(ID,2*(iterNum-1) + 42) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
            record(ID,6) = iterNum;
            record(ID,7) = f_opt2;
            record(ID,8) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = ( beta1*u1' + beta2*u2' + rho_const*ubar1*u1' + rho_const*ubar2*u2' );
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';

        beta1 = beta1 + rho_const*(ubar1 - B*u1);
        beta2 = beta2 + rho_const*(ubar2 - B*u2);

        %rho_const = rho_const + 0.1*rho_const; 
    

    
    end
    [real_opt,~,~] = PrimaryNewsvendorUpperGivenB( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,m1, B);
    record(ID,9) = real_opt;

    [X,s,Landa1, Landa2, qrr, Qr,real_opt,~] = PrimaryNewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B);
    record(ID,10) = real_opt;


    ID = ID + 1;
    continue; % TODO

    m1=1;
    
    % init B 
    B = zeros(m,m1);
    % randomRows = randperm(m);
    for i = 1 : m1
        B(i,i) = 1; %i  or m-m1+i or randomRows(i)
    end
    
    theta = {};
    theta{1} = zeros(m, 1);
    theta{2} = zeros(m, 1);

    rho_const = 100;

    totalTime = 0;
    lastOptimalValue = f_opt1;
    for iterNum = 1:500
        [ p,tp,t,P,f_opt2,CPUTime2] = ADMMnewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B, theta, rho_const);
%                 record(ID,10 + 6*(iterNum-1) + 1) = f_opt2;
%                 record(ID,10 + 6*(iterNum-1) + 2) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        %return;

       
        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
%                     record(ID,10 + 6*(iterNum-1) + 5) = iterNum;
%                     record(ID,10 + 6*(iterNum-1) + 6) = totalTime;
            record(ID,11) = iterNum;
            record(ID,12) = f_opt2;
            record(ID,13) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = -( theta{1}*p{1}' + theta{2}*p{2}' - rho_const*tp{1}*p{1}' - rho_const*tp{2}*p{2}' )  ;
        [Umat,~,Vmat] = svd(SofB, "econ");
        B = Umat * Vmat';


        theta{1} = theta{1} - rho_const*(tp{1} - B*p{1});
        theta{2} = theta{2} - rho_const*(tp{2} - B*p{2});
    end

    [X,s,Landa1, Landa2, qrr, Qr,real_opt,~] = PrimaryNewsvendorGivenB(ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B);
    record(ID,14) = real_opt;






    m1=2;
    
    % init B 
    newB = zeros(m,m1);
    % randomRows = randperm(m);
    for i = 1 : m1
        newB(i,i) = 1; %i  or m-m1+i or randomRows(i)
    end
    
    beta1 = zeros(m, 1);
    beta2 = zeros(m, 1);

    rho_const = 500;

    totalTime = 0;
    lastOptimalValue = f_opt1;
    for iterNum = 1:500
        [X,s,Landa1,Landa2,q,Qr,u1,u2,ubar1,ubar2, f_opt2, CPUTime2 ] = PBlowerboundGivenB( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,1, newB, beta1, beta2, rho_const);
%                 record(ID,10 + 6*(iterNum-1) + 1) = f_opt2;
%                 record(ID,10 + 6*(iterNum-1) + 2) = CPUTime2;
        totalTime = totalTime + CPUTime2; 

        %return;

       
        if abs(f_opt2 - lastOptimalValue) / abs(lastOptimalValue) <= 0.001
%                     record(ID,10 + 6*(iterNum-1) + 5) = iterNum;
%                     record(ID,10 + 6*(iterNum-1) + 6) = totalTime;
            record(ID,16) = iterNum;
            record(ID,17) = f_opt2;
            record(ID,18) = totalTime;
            break;
        else
            lastOptimalValue = f_opt2;
        end

        SofB = ( beta1*u1' + beta2*u2' + rho_const*ubar1*u1' + rho_const*ubar2*u2'   );
        %SofBbar =  beta3*h1' + beta4*h2' + rho_const*hbar1*h1' + rho_const*hbar2*h2';
        [Umat,~,Vmat] = svd(SofB, "econ");
        %[Umat2,~,Vmat2] = svd(SofBbar, "econ");
        newB = Umat * Vmat';
        %Bbar = Umat2 * Vmat2';


        beta1 = beta1 + rho_const*(ubar1 - newB*u1);
        beta2 = beta2 + rho_const*(ubar2 - newB*u2);

    end

    [X,s,Landa1, Landa2, qrr, Qr,real_opt,~] = PrimaryNewsvendorGivenB(ZIGMA,mu,1,A,b,gamma1,gamma2,c,v,g, newB(:,1));
    record(ID,19) = real_opt;





            
    ID = ID + 1;
        
   
    
end




