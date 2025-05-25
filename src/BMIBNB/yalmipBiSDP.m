clear;
clc;


 
load("./Data/200-1.mat") 
[l,m] = size(A);
m=length(mu);
n = m;

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);

m1 = 2;

tic

yalmip('clear')
X = sdpvar(n, 1); 
q = sdpvar(m, 1);
s = sdpvar(1,1); 
Qr = sdpvar(m1, m1);
Landa1 = sdpvar(l, 1);
Landa2 = sdpvar(l, 1);
u1 = sdpvar(m1, 1);
u2 = sdpvar(m1, 1);
ubar1 = sdpvar(m, 1);
ubar2 = sdpvar(m ,1); 

B = sdpvar(m, m1, 'full');


F = [X>=0, Landa1>=0, Landa2>=0]; 
F = [F, q + A1Complete'*(A'*Landa1) == ubar1];
F = [F, q + A1Complete'*(A'*Landa2+(v-g)) == ubar2]; 
SDPcons1 = [s-(c-v)'*X-(Landa1')*(b-A*mu)  0.5* u1'; 0.5* u1  Qr]; 
F = [F, SDPcons1>=0];
SDPcons2 = [s-((c-v)'+(v-g)')*X-Landa2'*(b-A*mu)+(v-g)'*mu  0.5* u2'; 0.5* u2  Qr]; 
F = [F, SDPcons2>=0];
F = [F, ubar1-B*u1 == 0]; 
F = [F, ubar2-B*u2 == 0]; 
F = [F, B'*B == eye(m1)];
 
options = sdpsettings('solver','bmibnb', 'bmibnb.lowersolver', 'mosek', 'bmibnb.maxtime', 7200, 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 7200);
% s + sum(sum((gamma2*eye(m1)).*Qr)) + sqrt(gamma1)*norm(q,2) 
optimize(F,s + sum(sum((gamma2*eye(m1)).*Qr)) + sqrt(gamma1)*norm(q,2),options);


CPUTime = toc;

disp(CPUTime)
 




