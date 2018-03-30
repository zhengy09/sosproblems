function [x,y,cost,info] = solveWithSCSdirect(At,b,C,K,params)

%% SOLVEWITHSCS.m Solve with scs
% 
% Solve conic program in sedumi format with SCS.


% Solve full problem with SCS, default settings
[data,K] = sedumi2scs(At,b,C,K);
[y,x,~,info] = scs_direct(data,K,params);
cost = -data.c'*y;
X = {};%cell(K.f+K.l+K.q+K.s,1);
X = blockify(X,x,K);
x = flatten(x,X,0);


