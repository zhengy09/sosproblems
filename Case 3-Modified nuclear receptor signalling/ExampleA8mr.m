%% A.8 SBToolbox file for Matlab containing the model of nuclear receptor dignalling
%  A modified version of the Reduced model in A.8 SBToolbox
%  for the original reduced model, see the following thesis
%
%  Khoshnaw, Sarbaz Hamza Abdullah. Model Reductions in Biochemical Reaction Networks. 
%  Diss. Department of Mathematics, 2015, Appendix A.8
%

clc;clear

% parameters
paraA8r

% variables
pvar HNr Lr L2Nr ABr Nr HNABr L2NABr GRr GRBr HNAr GRACr Gr GRR1r RR1r R1r DRR1r Fr
pvar NFr ENr NTr DRFr DRr GRBNr NCr GRBP1r L2NRr L4NRr H2NRr H3NRr GREr HNEr L2NTr GRTr
pvar HNTr GRHNEr GRL2NTr GRHNTr

x = [HNr; Lr; L2Nr; ABr; Nr; HNABr; L2NABr; GRr; GRBr; HNAr; GRACr; Gr; GRR1r; RR1r; R1r; DRR1r; Fr; ...
    NFr; ENr; NTr; DRFr; DRr; GRBNr; NCr; GRBP1r; L2NRr; L4NRr; H2NRr; H3NRr; GREr; HNEr; L2NTr; GRTr;...
    HNTr; GRHNEr; GRL2NTr; GRHNTr];

u1  = k1*HNr*Lr;
u3  = k3*HNr*ABr;
u4  = k4*L2Nr*ABr;
u8  = k8*L2Nr;
u10 = k10*HNABr*GRr;
u11 = k11*GRr*HNAr*(0.6-GRACr-NCr);
u13 = k13*GRr*(0.6-GRACr-NCr);
u15 = k15*GRr*(0.6-GRACr-NCr)*0.4;
u16 = k16*HNr*Lr;
u18 = k18*L2Nr*(0.6-L2NRr-L4NRr-H2NRr-H3NRr);
u19 = k19*L2Nr*L2NRr;
u20 = k20*HNr*(0.6-L2NRr-L4NRr-H2NRr-H3NRr);
u22 = k22*HNr*H2NRr;
u24 = k24*HNr*0.9;
u25 = k25*GRr*HNTr;
u26 = k26*HNr*GRTr;
u27 = k27*GRr*0.9;
u29 = k29*L2Nr*0.9;
u30 = k30*GRr*L2NTr;
u31 = k31*HNr*0.9;
u32 = k32*GRr*HNEr;
u33 = k33*GRr*0.9;
u34 = k34*HNr*GREr;
u41 = k41*GRHNEr;
u42 = k42*GRL2NTr;
u45 = k45*HNEr;
u46 = k46*L2NTr;
u48 = k48*GRHNTr;
u49 = k49*HNTr;
u52 = k52*GRR1r;
u53 = k53*Gr*RR1r;
u54 = k54*DRR1r;
u55 = k55*R1r*DRFr;
u56 = k56*Nr*Fr;
u59 = k59*Nr*(0.6-GRACr-NCr);
u62 = k62*NCr;
u63 = k63*Nr*GRBr;
u65 = k65*NFr;
u66 = k66*Nr*0.9;
u68 = k68*NTr;
u70 = k70*GRBr*0.8;
u71 = k71*0.4*0.3;
u72 = k72*0.4*GRBP1r;
u73 = k73*GRBP1r;
u74 = k74*GRr;
u75 = k75*GRr*(0.6-GRACr-NCr)*0.4;
u77 = k77*GRACr*0.8;
u79 = k79*Fr*DRr;

% dynamics 
f(1)  = -u1-u3+u11-u16-u20-u22-u24-u26-u31-u34+u45+u49;
f(2)  = -u1-u16 - k1*x(2);
f(3)  = u1-u4-u8+u13+u15+u16-u18-u19-u29+u46;
f(4)  = -u3-u4+u72 -k3*x(4);
f(5)  = -u56-u59+u62-u63+u65-u66+u68;
f(6)  = u3-u10 -k3*x(6);
f(7)  = u4+u8 -k4*x(7);
f(8)  = -u10-u11-u13-u15-u25-u27-u30-u32-u33+u41+u42+u48+u52+u73-u74-u75;
f(9)  = u10-u63-u70;
f(10) = u10-u11 -k10*x(10);
f(11) = u11+u13+u15+u75-u77;
f(12) = -u53 - k53*x(12);
f(13) = -u52+u53;
f(14) = -u53+u54 - k53*x(14);
f(15) = u52-u55 - k52*x(15);
f(16) = -u54+u55;
f(17) = u55-u56+u65-u79 - k55*x(17);
f(18) = u56-u65;
f(19) = u66 - k66*x(19);
f(20) = -u68;
f(21) = -u55+u79 - k55*x(21);
f(22) = u74+u77-u79 - k74*x(22);
f(23) = u63 - k63*x(23);
f(24) = u59-u62;
f(25) = u70-u72-u73;
f(26) = u18-u19 - k18*x(26);
f(27) = u19 - k19*x(27);
f(28) = u20-u22 - k20*x(28);
f(29) = u22 - k22*x(29);
f(30) = u33-u34 - k20*x(30);
f(31) = u31-u32+u41-u45;
f(32) = u29-u30+u42-u46;
f(33) = -u26+u27 - k26*x(33);
f(34) = u24-u25+u48-u49;
f(35) = u32+u34-u41;
f(36) = u30-u42;
f(37) = u25+u26-u48;

f0 = double(subs(f,x,zeros(37,1)));               %% origin is an equilibrium point


%% build a quadratic Lyapunov equation
Maxiter = 2000;
pro = tic;
% local region; g(x) <= 0
gamma = 0.01;        % a ball with a radius of 0.1
g = sum(x.^2)-gamma;

%% find Lyapunov function
Degree = 2;
% =============================================
% First, initialize the sum of squares program
prog = sosprogram(x);

% =============================================
% The Lyapunov function V(x) : 
CandidateMonomails = monomials(x,[0:Degree]);         % full candidate
[prog,V] = sospolyvar(prog,CandidateMonomails);     

% =============================================
% Next, define SOSP constraints
kappa1 = 0.1;
PolyIneq1 = V - kappa1*sum(x.^2);
prog = sosineq(prog,PolyIneq1);

% Constraint 2: -dV/dx*f + r*g >= 0
[prog,r] = sospolyvar(prog,CandidateMonomails); 
prog = sosineq(prog,r);

PolyIneq2 = r*g;
for i = 1:length(x)
    PolyIneq2 = PolyIneq2 - diff(V,x(i))*f(i);
end
prog = sosineq(prog,PolyIneq2);

Time = toc(pro);

% =============================================
% And call solver
options.solver          = 'CDCS';
options.params.solver   = 'sos';
options.params.relTol   = 1e-4;
options.params.dispIter = 200;
options.params.maxIter  = 2000;
prog = sossolve(prog,options);
SOLV = sosgetsol(prog,V);
SOLR = sosgetsol(prog,r);

%% Check the Lyapunov function 
Flag = CheckLyapunov(f,x,SOLV,SOLR,g,kappa1);