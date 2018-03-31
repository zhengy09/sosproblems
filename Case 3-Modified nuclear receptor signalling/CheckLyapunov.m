function Flag = CheckLyapunov(f,x,V,r,g,kappa)
% verify the Lyapunov Candidate

    PolyIneq1 = V - kappa*sum(x.^2);
    PolyIneq2 = r*g;
    for i = 1:length(x)
        PolyIneq2 = PolyIneq2 - diff(V,x(i))*f(i);
    end
    

    options.solver = 'cdcs';
    options.params.solver = 'sos';
    options.params.relTol = 1e-4;
    options.params.maxIter = 2e3;
    
    prog1 = sosprogram(x); 
    prog1 = sosineq(prog1,PolyIneq1);
    prog1 = sossolve(prog1,options);
    
    prog2 = sosprogram(x); 
    prog2 = sosineq(prog2,PolyIneq2);
    prog2 = sossolve(prog2,options);
    
    if prog1.solinfo.info.problem == 0 && prog2.solinfo.info.problem == 0  %% Ture Lyapunov function
        Flag = 1;
    elseif prog1.solinfo.info.problem == 0 && prog2.solinfo.info.problem ~= 0
        Flag = 2;
    elseif prog1.solinfo.info.problem ~= 0 && prog2.solinfo.info.problem == 0
        Flag = 3;
    else
        Flag = 4;
    end

end

