clear;
clc;


 

m = 3; 

% mu = [1;2;3];
mu = [0.5;1.5;2.5];
% ZIGMA = [1 0 0;0 3 0; 0 0 2];
ZIGMA = [1 0.2 0.1;0.2 3 0.15; 0.1 0.15 2];

A = [-1 0 0;1 0 0;0 -1 0;0 1 0;0 0 -1;0 0 1];

b = [0;
    8;
    -1;
    12;
    -2;
    16];
% b = [-1;7;-2;11;-3;15];
ALPHA = 0.05;

B1 = [1 0 0; 0 0 0; 0 0 0];
B2 = [0 0 0; 0 1 0; 0 0 0];
B3 = [0 0 0; 0 0 0; 0 0 1];

B4 = [1 0 0; 0 1 0; 0 0 1];


[ f_opt4,X4,t4,CPUTime4 ]  = xxxCVaRPrimary(ZIGMA,mu,A,b,ALPHA, m, B4);
[ f_opt1,X1,t1,CPUTime1 ]  = xxxCVaRPrimary(ZIGMA,mu,A,b,ALPHA, m, B1);
[ f_opt2,X2,t2,CPUTime2 ]  = xxxCVaRPrimary(ZIGMA,mu,A,b,ALPHA, m, B2);
[ f_opt3,X3,t3,CPUTime3 ]  = xxxCVaRPrimary(ZIGMA,mu,A,b,ALPHA, m, B3);


disp([f_opt4,f_opt1,f_opt2,f_opt3]);
disp([X4,X1,X2,X3]);
disp([t4,t1,t2,t3]);

[ A1Complete,U,delta ] = covTransformerDecomposer(ZIGMA);
disp(delta);

 

         



