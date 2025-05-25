
clear;
clc;
for iterateNumber = 1:5
    
    
    
    n = xxx;  %TODO
    m = n; 
    
    gamma1 = 1;  
    gamma2 = 2;  
    
    mu = 10 * rand(m,1);  
    SDcosi = 1 * rand(m,1) + 1; 
    MULTI = (SDcosi)*(SDcosi');
    ZIGMA1 = gallery('randcorr',m); 
    ZIGMA = ZIGMA1.*MULTI;  
    
    ZIGMA = eye(m); 
    for i = 1:m
        ZIGMA(i,i) = 1+rand(1)*0.2; %TODO
    end
        
    I= eye(m);
    A= [I ; -1*I]; %   ones(1,m);-ones(1,m)];
    b = [2*SDcosi+mu; 2*SDcosi-mu];% sum(mu)+3.5; -(sum(mu)-3.5)];
     
    c = zeros(n,1); % to save purchase (wholesale) prices
    v = zeros(n,1); % to save selling (retail) prices
    g = zeros(n,1); % to save salvage prices
        
    for i = 1 : n
        c(i,1) = 0.1*(5+i-1);
        v(i,1) = 0.15*(5+i-1);
        g(i,1) = 0.05*(5+i-1);
    end
    
    save("./Data/" + n + "-" + iterateNumber);
end

