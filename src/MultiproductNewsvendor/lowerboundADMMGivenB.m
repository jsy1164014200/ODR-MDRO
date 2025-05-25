function [X,s,Landa1,Landa2,q,Qr,u1,u2,ubar1,ubar2, f_opt, CPUTime ] = lowerboundADMMGivenB( ZIGMA,mu,A,b,gamma1,gamma2,c,v,g,m1, newB, beta1, beta2, rho_const)
%This function finds the original OPT value (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

[l,m] = size(A);
m=length(mu);
n = m;

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);



tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables q(m) s;
                variable Qr;
                variable Landa1(l) nonnegative;
                variable Landa2(l) nonnegative;
                variables u1(2) u2(2) ubar1(m) ubar2(m);

%                 variables h1 h2 hbar1(m) hbar2(m);
                 
                
                                
                %% 
                minimize( s +  gamma2*Qr + sqrt(gamma1)*norm(q,2) +beta1'*(ubar1-newB*u1) +beta2'*(ubar2-newB*u2)  + (rho_const/2)*(ubar1-newB*u1)'*(ubar1-newB*u1) + (rho_const/2)*(ubar2-newB*u2)'*(ubar2-newB*u2)    );
                %%
                subject to
                [s-(c-v)'*X-(Landa1')*(b-A*mu)  0.5* u1(1)'; 0.5* u1(1)  Qr] >= 0;
                [s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu  0.5* u2(1)'; 0.5* u2(1)  Qr] >= 0;
                q + A1Complete'*(A'*Landa1) == ubar1 ;
                q + A1Complete'*(A'*Landa2+(v-g)) == ubar2;
                

                
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

