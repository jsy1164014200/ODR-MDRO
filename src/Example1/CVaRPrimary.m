function [ f_opt,X,t,CPUTime ] = CVaRPrimary(ZIGMA,mu,A,b,ALPHA, m, B)
 
 
 [l,~] = size(A);

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
A1Complete = A1Complete*B;

tic 
cvx_begin sdp
 
                variable X(m) nonnegative;
                variable t;
                variables q(m) s;
                variable Q(m,m) symmetric;
                variable Landa1(l) nonnegative;
                variable Landa2(l) nonnegative;
                
                                
                %%
                minimize( s + sum(sum( eye(m).*Q )) );
                %%
                subject to
                
                [s-t-(Landa1')*(b-A*mu)     1/2*(q + A1Complete'*(A'*Landa1))' ; 1/2*(q + A1Complete'*(A'*Landa1)) Q]>=0 ;
                [s-(1-1/ALPHA)*t-Landa2'*(b-A*mu)-1/ALPHA*X'*mu    1/2*(q +A1Complete'*(A'*Landa2-1/ALPHA*X))' ; 1/2*(q +A1Complete'*(A'*Landa2-1/ALPHA*X)) Q]>=0;
                sum(X) == 1; 
%                 X(1) == 1;
                
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

