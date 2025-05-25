function [B] = Insight(ZIGMA,v,g,m)
%INSIGHT Summary of this function goes here
%   Detailed explanation goes here
[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

para = (g-v)'*A1Complete;


cvx_begin sdp
 
                variable B(m,2);

                p = para*B;

                                
                %%
                maximize(  p*p'  );
                %%
                subject to

                
                [eye(m)  B; B' eye(2)] >= 0;
               
 cvx_end
    
 
    disp(['Problem is ' cvx_status])
    if ~strfind(cvx_status,'Solved')
      return
    end
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['optimal value of cvx:',num2str(cvx_optval)]);
    f_opt=cvx_optval;





end

