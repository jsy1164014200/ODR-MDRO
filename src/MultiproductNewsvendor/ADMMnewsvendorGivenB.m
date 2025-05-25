function [ p,tp,t,P,f_opt,CPUTime ] = ADMMnewsvendorGivenB( ZIGMA,mu,m1,A,b,gamma1,gamma2,c,v,g, B, theta, rho_const)
%This function finds the LB OPT value (General moment-based ambiguity
%set)for Multiproduct Newsvendor problem.

m=length(mu);
n = m;

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

tic 

cvx_begin sdp
 
                variables p1(m1) p2(m1);
                variables tp1(m) tp2(m);
                variable P1(m1,m1) symmetric;
                variable P2(m1,m1) symmetric;
                variables t1 t2;

                

                                
                %%- (beta_const/2)*(p1-curp{1})'*(p1-curp{1})- (beta_const/2)*(p2-curp{2})'*(p2-curp{2})
                maximize( (t2*mu' + tp2'*A1Complete')*(g-v) + theta{1}'*(tp1 - B*p1)+ theta{2}'*(tp2 - B*p2)- (rho_const/2)*(tp1 - B*p1)'*(tp1 - B*p1) - (rho_const/2)*(tp2 - B*p2)'*(tp2 - B*p2)   );
                %%
                subject to

                t1 + t2 == 1;
                t1*(A*mu-b)' + tp1'*A1Complete'*A' <= 0;
                t2*(A*mu-b)' + tp2'*A1Complete'*A' <= 0;
                sqrt(gamma1) - norm(p1+p2,2) >= 0;
                gamma2*eye(m1) - (P1+P2) >= 0;
                t1*(c-v)' + t2*(c-g)' >= 0;
                [t1   p1'; p1  P1] >= 0;
                [t2   p2'; p2  P2] >= 0;
                

                               

 cvx_end
    
CPUTime = toc;  
 
 
    disp(['Problem is ' cvx_status])
    if ~strfind(cvx_status,'Solved')
      return
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['optimal value of cvx:',num2str(cvx_optval)]);
    p = {};
    p{1} = p1;
    p{2} = p2;
    tp = {};
    tp{1} = tp1;
    tp{2} = tp2;
    t = {};
    t{1} = t1;
    t{2} = t2;
    P = {};
    P{1} = P1;
    P{2} = P2;

    f_opt=cvx_optval;

end

