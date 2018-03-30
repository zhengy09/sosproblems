
%% Minimize the quartic polynomial over the hypercube
% f(x) = \sum_{1<i<j<n} xixj + xi^2xj-xj^3-xi^2xj^2
% gi(x) = \sum xi^2<1
clc;clear

N = round(logspace(log10(10),log10(50),10));
d = 2;

TimePro   = zeros(length(N),1);    % time for problem generation
TimeTotal = zeros(length(N),7);    % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
TimeSetup = zeros(length(N),3);    % scs-direct, scs-indirect, cdcs-sos
TimeADMM  = zeros(length(N),3);    % scs-direct, scs-indirect, cdcs-sos
TimeAver  = zeros(length(N),3);    % scs-direct, scs-indirect, cdcs-sos
Cost      = zeros(length(N),7);    % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
Iter      = zeros(length(N),6);    % sedumi, sdpt3, sdpa, scs-direct, scs-indirect, cdcs-sos
Density   = zeros(length(N),4);

data = cell(length(N),1);
%%
Maxiter = 2e3;
Tol     = 1e-3;
TolIPM  = 1e-8;

for Index = 1:length(N)
    n = N(Index);
    fprintf('testing : N = %i', N(Index))
    tdata = tic;
    %% generating POP via GloptiPloy 3
    mpol('x',n,1)
    f = 0;
    for i = 1:n
        for j = i:n
            f = f + x(i)*x(j) + x(i)^2*x(j) - x(j)^3 - x(i)^2*x(j)^2;
        end
    end
    g = 0;
    for i = 1:n
        g = g+x(i)^2;
    end
    K = [g <=1];
    P = msdp(min(f), K, d);
    [A,b,c,K] = msedumi(P);

    %% Problem data in Sedumi form
    [m,n] = size(A);
    Density(Index,:) = [m,n,sum(sum(spones(A)))/m/n,(min(K.s)*(min(K.s)+1))/2];
    
    TimePro(Index) = toc(tdata);
    At = A';
    data{Index}.At = At;
    data{Index}.b  = b;
    data{Index}.c  = c;
    data{Index}.K  = K;
    
    %% solutions using different method
    % 1 by sedumi
    x1 = zeros(length(c),1);
    x2 = zeros(length(c),1);
    x3 = zeros(length(c),1);
    try
        opts.eps = TolIPM;
        [x1,y1,infoSedumi] = sedumi(At,b,c,K,opts);
    catch
        warning('out of memory')
    end
    
    % 2 by sdpt3
    try
        [blk,At2,C2,b2]      = read_sedumi(At,b,c,K); 
        tsdpt = tic;
        opts.gaptol = TolIPM;
        [obj,X,y,Z,infoSDPT] = sqlp(blk,At2,C2,b2,opts); 
        timeSdpt = toc(tsdpt);
    catch
        warning('out of memory')
    end
    
    % 3 by sdpa
    if N(Index) < 29  % N = 29, matlab on my desktop crashed
    try 
        if isfield(K,'q') && ~isempty(K.q) && K.q == 0
            K.q = [];  
            tSdpa = tic;
            opts.epsilonStar = TolIPM;
            [x3,y3,infoSdpa] = sedumiwrap(A,b,c,K,[],opts);
            timeSdpa = toc(tSdpa);
        else
            
        end
    catch
        warning('out of memory')
    end
    end
    
    % 4 by csdp
    try
        opts.axtol  = TolIPM;
        opts.aytol  = TolIPM;
        opts.objtol = TolIPM;
        tCsdp       = tic;
        [x4,y4,z4,info_csdp] = csdp(At,b,c,K,opts); 
        timeCsdp = toc(tCsdp);
    catch
        warning('out of memory')
    end
    
    % 5 by SCS-direct
    params.max_iters = Maxiter;
    params.eps       = Tol;
    [x5,y5,cscs5,infoSCSdirect] = solveWithSCSdirect(At,full(b),full(c),K,params);
    
    % 6 by SCS-indirect
    params.max_iters = Maxiter;
    params.eps = Tol;
    [x6,y6,cscs6,infoSCSindirect] = solveWithSCSindirect(At,full(b),full(c),K,params);
    

    % 7 by cdcs - sos
    opts.relTol = Tol;
    opts.solver = 'sos';
    opts.maxIter = Maxiter;
    [x7,y7,z7,infoCDCSsos] = cdcs(At,b,c,K,opts);
    
    
    %% statistics
    TimeTotal(Index,:) = [infoSedumi.wallsec,timeSdpt,timeSdpa,timeCsdp, ...
                          (infoSCSdirect.solveTime+infoSCSdirect.setupTime)/1e3, ...
                          (infoSCSindirect.solveTime+infoSCSindirect.setupTime)/1e3, ...
                          infoCDCSsos.time.total];   
    TimeSetup(Index,:) = [infoSCSdirect.setupTime/1e3,infoSCSindirect.setupTime/1e3,infoCDCSsos.time.setup]; 
    TimeADMM(Index,:) = [infoSCSdirect.solveTime/1e3,infoSCSindirect.solveTime/1e3,infoCDCSsos.time.admm]; 
    Cost(Index,:) = [c'*x1,obj(1),c'*x3,c'*x4,cscs5,cscs6,c'*x7];
    Iter(Index,:) = [infoSedumi.iter,infoSDPT.iter,infoSdpa.iteration,infoSCSdirect.iter,infoSCSindirect.iter,infoCDCSsos.iter];
    TimeAver(Index,:)  = TimeADMM(Index,:)./Iter(Index,4:end);
    
    save QuarticPoly.mat
end


