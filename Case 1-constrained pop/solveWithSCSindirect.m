function [x,y,cost,info] = solveWithSCSindirect(At,b,C,K,params)

%% SOLVEWITHSCS.m Solve with scs
% 
% Solve conic program in sedumi format with SCS.


% Solve full problem with SCS, default settings
[data,K] = sedumi2scs(At,b,C,K);
[x,y,~,info] = scs_indirect(data,K,params);
cost = -data.c'*x;

