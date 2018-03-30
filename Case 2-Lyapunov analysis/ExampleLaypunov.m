
%% Testing Lyapunov functions

%% create a nonlinear dynamic system
clc;clear   
n = round(logspace(log10(10),log10(50),10)); %% number of variable in x, x1,x2,...,xn.
d = 3;              %% degree of fi(x)
density = 0.2;
rho     = 0.2;

warning on;
Maxiter = 2e3;
Tol     = 1e-3;

TimePro   = zeros(length(n),1);     % time for problem generation
TimeSOS   = zeros(length(n),1);     % time for problem generation
TimeTotal = zeros(length(n),7);     % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
TimeSetup = zeros(length(n),3);     % scs-direct, scs-indirect, cdcs-sos
TimeADMM = zeros(length(n),3);      % scs-direct, scs-indirect, cdcs-sos
TimeAver = zeros(length(n),3);      % scs-direct, scs-indirect, cdcs-sos
Iter = zeros(length(n),6);
Flag = zeros(length(n),4);          % whether the resulting Lyapunov function is valid?

InfoTerm = zeros(length(n),4)+NaN;

DynamicsModel = cell(length(n),1);  % Dynamics
SOSProg = cell(length(n),1);        % SOSTOOLS prog

SeDuMiData = cell(length(n),1);  
Dimension  = zeros(length(n),3);    % n,m, K.f,


for i = 1:length(n)

    fprintf('Testing n(i) = %i \n', n(i));
    R = sprand(n(i),n(i),density);  %% generate random sparse matrix, determine the pattern of f(x).
    R = full(spones(R));
    tic
    [f,x] = GenNonDynamics(R,rho,d); 
    TimePro(i) = toc;
    DynamicsModel{1}.f = f;
    DynamicsModel{1}.x = x;

    % local region; g(x) <= 0
    gamma = 0.01;     % a ball with a radius of 0.1
    g = sum(x.^2)-gamma;
    
    fprintf('SOSTOOLS processing ... \n');
    tic
    %% find Lyapunov function, a quardractic candidate
    Degree = 2;
    % =============================================
    % First, initialize the sum of squares program
    prog = sosprogram(x);

    % =============================================
    % The Lyapunov function V(x): 
    CandidateMonomails = monomials(x,[2:Degree]);         % candidate monomials
    [prog,V] = sospolyvar(prog,CandidateMonomails); 

    % =============================================
    % Next, define SOSP constraints
    % Constraint 1 : V(x) - 0.01*(x1^2 + x2^2 + x3^2 + ... + xn^2)>= 0
    kappa1 = 0.01;
    PolyIneq1 = V - kappa1*sum(x.^2);
    prog = sosineq(prog,PolyIneq1);

    % Constraint 2: -dV/dx*f + r*g >= 0;
    [prog,r] = sospolyvar(prog,CandidateMonomails); 
    prog = sosineq(prog,r);
    PolyIneq2 = r*g;
    for j = 1:length(x)
        PolyIneq2 = PolyIneq2 - diff(V,x(j))*f(j);
    end
    prog = sosineq(prog,PolyIneq2);
    
    TimeSOS(i) = toc;
    fprintf('SOSTOOLS processed in %6.2f\n', TimeSOS(i));
    SOSProg{i} = prog;
    % =============================================     
    

    %% by SeDuMi
    %if i == 1
    %if n(i) < 24
    try
        prog.solinfo.x = [];
        opts.solver = 'sedumi';
        prog1 = sossolve(prog,opts);
        SOLV = sosgetsol(prog1,V);
        SOLR = sosgetsol(prog1,r);
        Flag1 = CheckLyapunov(f,x, SOLV,SOLR,g,kappa1);
    catch
        prog1.solinfo.info.wallsec = NaN;
        prog1.solinfo.info.iter    = NaN;
        prog1.solinfo.info.feasratio    = NaN;
        warning('out of memory');
    end
    %end
    
   
    %% call SCS solver
    prog.solinfo.x = [];
    options.solver = 'scs-indirect';
    options.params.max_iters = Maxiter;
    options.params.eps = Tol;
    prog5 = sossolve(prog,options);
    SOLV = sosgetsol(prog5,V);
    SOLR = sosgetsol(prog5,r);
    Flag5 = CheckLyapunov(f,x, SOLV,SOLR,g,kappa1);

    %% call SCS solver
    prog.solinfo.x = [];
    options.solver = 'scs-direct';
    options.params.max_iters = Maxiter;
    options.params.eps = Tol;
    prog6 = sossolve(prog,options);
    SOLV = sosgetsol(prog6,V);
    SOLR = sosgetsol(prog6,r);
    Flag6 = CheckLyapunov(f,x, SOLV,SOLR,g,kappa1);

    %% by CDCS-SOS
    prog.solinfo.x = [];
    options.solver = 'cdcs';
    options.params.solver = 'sos';
    options.params.relTol = Tol;
    options.params.maxIter = Maxiter;
    prog7 = sossolve(prog,options);
    SOLV = sosgetsol(prog7,V);
    SOLR = sosgetsol(prog7,r);
    Flag7 = CheckLyapunov(f,x, SOLV,SOLR,g,kappa1);
    
    SeDuMiData{i} = prog7.solinfo.info.SeDuMiData;
    At = SeDuMiData{i}.At;
    b = SeDuMiData{i}.b;
    c = SeDuMiData{i}.c;
    K = SeDuMiData{i}.K;
    [nn,mm] = size(SeDuMiData{i}.At);
    Dimension(i,:) = [nn,mm,SeDuMiData{i}.K.f];
    
    % 2 by sdpt3
    timeSdpt = NaN;
    obj = NaN;
    try
        [blk,At2,C2,b2]      = read_sedumi(At,b,c,K); 
        tsdpt = tic;
        [obj,X,y,Z,infoSDPT] = sqlp(blk,At2,C2,b2,opts); 
        timeSdpt = toc(tsdpt);
    catch
        warning('out of memory')
    end

    x3 = ones(length(c),1);
    x  = ones(length(c),1);
    timeSdpa = NaN;
    timeCsdp = NaN;
    if n(i) < 29
        % 3 by sdpa
        %if n(i) < 29  % N = 29, matlab on my computer crashed
        try  % sedumiwrap cannot handle k.q
            if ~isfield(K,'q') || isempty(K.q) ||(~isempty(K.q) && K.q == 0)
                K.q = [];  
                tSdpa = tic;
                [x3,y3,infoSdpa] = sedumiwrap(At',b,c,K,[],opts);
                timeSdpa = toc(tSdpa);
            end
        catch
            warning('out of memory')
        end


        % 4 by csdp  N = 29, storage allocation failed 
        try
            if (isfield(K,'f')) %Convert free vars to non-negative LP vars
                n_free = K.f; 
                [A1,b1,c1,K1] = convertf(At,b,c,K); %K.f set to zero
                At1 = A1';
            end
            c1 = full(c1);
            tCsdp    = tic;
            [x,y,z,info_csdp] = csdp(At1,b1,c1,K1);  %JA updated handling of info flag
            timeCsdp = toc(tCsdp);
        catch
            warning('out of memory')
        end
     end

    %% record time
    TimeTotal(i,:) = [prog1.solinfo.info.wallsec, timeSdpt,timeSdpa,timeCsdp,...
        (prog5.solinfo.info.setupTime + prog5.solinfo.info.solveTime)/1e3,(prog6.solinfo.info.setupTime + prog6.solinfo.info.solveTime)/1e3,prog7.solinfo.info.time.total];   
    TimeSetup(i,:) = [(prog5.solinfo.info.setupTime)/1e3,(prog6.solinfo.info.setupTime)/1e3,prog7.solinfo.info.time.setup];  
    TimeADMM(i,:) = [(prog5.solinfo.info.solveTime)/1e3,(prog6.solinfo.info.solveTime)/1e3,prog7.solinfo.info.time.admm]; 
    Iter(i,:) = [prog1.solinfo.info.iter,infoSdpa.iteration, infoSDPT.iter, ...
                prog5.solinfo.info.iter,prog6.solinfo.info.iter,prog7.solinfo.info.iter]; 
    TimeAver(i,:) = TimeADMM(i,:)./Iter(i,4:end);
    Flag(i,:) = [Flag1,Flag5,Flag6,Flag7]; 
    
    InfoTerm(i,:) = [prog1.solinfo.info.feasratio,infoSDPT.termcode, infoSdpa.dualityGap,prog7.solinfo.info.problem] ; 

    save ResultRandom
end


