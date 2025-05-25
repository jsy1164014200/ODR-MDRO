function [X,s,Landa1, Landa2, qrr, Qr, f_opt,CPUTime] = PrimaryNewsvendorGivenB( ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B)

[l,m] = size(A);

m=length(mu);
n = m;

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

A1 = A1Complete*B;



tic 

cvx_begin sdp
 
                variable X(n) nonnegative;
                variables qrr(m1) s;
                variable Qr(m1,m1) symmetric;
                variable Landa1(l) nonnegative;
                variable Landa2(l) nonnegative;

                                
                %%
                minimize( s + sum(sum((gamma2*eye(m1)).*Qr)) + sqrt(gamma1)*norm(qrr,2) );
                %%
                subject to

                [s-(c-v)'*X-(Landa1')*(b-A*mu)     1/2*(qrr+ A1'*(A'*Landa1))' ; 1/2*(qrr + A1'*(A'*Landa1)) Qr]>=0 ;
                [s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu    1/2*(qrr +A1'*(A'*Landa2+(v-g)))' ; 1/2*(qrr +A1'*(A'*Landa2+(v-g))) Qr]>=0;
                
               
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

