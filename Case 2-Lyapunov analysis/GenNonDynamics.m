function [f,x] = GenNonDynamics(R,rho,d)
% create a polynomial system with banded correlative sparsity R and asymptotically stable
% linearisation - to give a locally asymptotically system.
% R: correlative sparsity pattern of the dynamics
% rho: controls the number of monomials in the dynamical system 
% d: degree of the nonlinear systems

n = size(R,1);
A = ceil(10*(rand(n,n)-0.5));
A = full(A.*R);

%% create local stable linear equations
lambda = max(real(eig(A)));
A = A -(ceil(lambda)+10)*eye(n);

%% variables x = [x1, x2, ..., xn];
for i=1:n
    x(i)=pvar(strcat('x',num2str(i)));
end

%% dynamics \dot x_i = f_i(x1, ..., xn)
for i=1:n
    f(i)=pvar(strcat('f',num2str(i)));
end

%% create dynamics for each node
for i=1:n
    R = logical(R);
    CandidateVariables = x(R(i,:));
    if isempty(CandidateVariables)
        CandidateVariables = x(i);
    end
    Z = monomials(CandidateVariables,2:d);  %% candidate monomails in fi(x)
    m = length(Z);
    tmp = ceil(rho*m);
    u = randperm(m,tmp);
    if isempty(u)
        f(i) = 0;  %% linear dyanmics
    else
        f(i) = ceil(10*(rand(1,tmp)-0.5))*Z(u);
    end
    f(i) = f(i) + A(i,:)*x.';
end
end
