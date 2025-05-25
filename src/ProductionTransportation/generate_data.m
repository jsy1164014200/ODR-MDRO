clear;
clc;
for iterateNumber = 1:3
    
    m=20; %%TODO The number of suppliers 
    n=30; %%TODO The number of customers 
    K = 5; %% TODO

    data_size = m*n;
    gamma1 = 1;
    gamma2 = 2;

    Faci=rand(2,m); %It has two rows because it shows the cordinations of supplier locations (X,Y)
    Dema=rand(2,n); %It has two rows because it shows the cordinations of demand locations(X,Y)
    for i=1:m
        for j=1:n
            dd((i-1)*n+j)=norm(Faci(:,i)-Dema(:,j)); 
            % It calculates the distance between any pair of demand and facility locations
            % the components of dd are \bar{\xi_{ij}} for all i and j
        end
    end

    Samp=zeros(10000,n*m);
    for t=1:10000
        Samp(t,:)=rand(1,n*m).*dd+0.5*dd; % Generating 10000 samples (Scenarios) of random variable \xi
    end

    mu = mean(Samp)';    % \mu of randome vectore \xi (It creates the mean vector (\mu) from 10000 sample \xi )
    ZIGMA = cov(Samp);   % It creates covariance matrix (\Sigma) from 10000 sample \xi
    SDcosi = std(Samp)'; % (MATLAB) S = std(A): If A is a matrix whose columns are random variables and whose rows are observations, then S is a row vector containing the standard deviations corresponding to each column.

    I= eye(data_size);
    A= [I ; -1*I];
    b_list={};
    b_list{1} = [2*SDcosi+mu; 2*SDcosi-mu];
    b_list{2} = [3*SDcosi+mu; 3*SDcosi-mu];
    b_list{3} = [4*SDcosi+mu; 4*SDcosi-mu];
    b = b_list{2};

    c=rand(m,1)*mean(mu)+0.5*mean(mu); % Production cost
    d=(rand(n,1)*0.5+0.5)*m/n;         % the amount of demand

    
    %%% To approximate the disutility function as a piece-wise linear function%%  
    xx=0:(1/K):1; % X = [0,.2,.4,.6,.8,1]
    yy=0.25*(exp(2*xx)-1); % Disutility function
    alpha_k=(yy(2:2+K-1)-yy(1:K))/0.2; % Slopes of disutility functiuon (a_1, a_2, ... , a_5)
    beta_k=yy(1:K)-alpha_k.*xx(1:K);  % Intercepts of disutility functiuon (b_1, b_2, ... , b_5)
    %%% end 
    
    save("./Data/" + K + "-" + m + "-" + n + "-" + iterateNumber);
end