function [optimalB, f_opt,X_opt,CPUTime ] = PrimaryNewsvendor( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g)
%This function finds the original OPT value (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

m=length(mu);
n = m;

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables q(m) s;
                variable Q(m,m) symmetric;
                variable Landa1(2*m) nonnegative;
                variable Landa2(2*m) nonnegative;
                
                                
                %%
                minimize( s + sum(sum((gamma2*eye(m)).*Q)) + sqrt(gamma1)*norm(q,2) );
                %%
                subject to
                
                [s-(c-v)'*X-(Landa1')*(b-A*mu)     1/2*(q+ A1Complete'*(A'*Landa1))' ; 1/2*(q + A1Complete'*(A'*Landa1)) Q]>=0 ;
                [s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu    1/2*(q +A1Complete'*(A'*Landa2+(v-g)))' ; 1/2*(q +A1Complete'*(A'*Landa2+(v-g))) Q]>=0;
                
 cvx_end
    
CPUTime = toc;  
 
 
    disp(['Problem is ' cvx_status])
    if ~strfind(cvx_status,'Solved')
      return
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['optimal value of cvx:',num2str(cvx_optval)]);
    f_opt=cvx_optval;
    X_opt=X;
    
    % [A1Complete'*(A'*Landa1)
    addOne = zeros(m,1);
    addOne(1) = 1;
    beforeSchmidt = [addOne, A1Complete'*(A'*Landa2+(v-g))];
    afterSchmidt = schmidt(beforeSchmidt);
    
    optimalB = afterSchmidt;
    
    


end

