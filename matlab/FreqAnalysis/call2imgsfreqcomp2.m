function outmat = call2imgsfreqcomp2(n_times)
    
    outmat = {};
    for ii = 1:n_times
    tic
    [baseIm, noiseData] = noiseBaseImDescFile(1);
    
    optParsStart = noiseData.pars2opt_start;
    lb = noiseData.pars2opt_lb;
    ub = noiseData.pars2opt_ub;
    stpSiz = noiseData.pars2opt_stpSiz;
    
    

    % fminsearch call **
    % #####################################################################
    %[x,fval] = fminsearch(@ImgSFreqComp_v5,optParsStart); %fminsearch..
    % #####################################################################

    % fmincon only call: **
    % #####################################################################
%     options = optimoptions(@fmincon,'FiniteDifferenceStepSize',stpSiz,'MaxIterations',150);
%     options
%     [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@ImgSFreqComp_v5,optParsStart,[],[],[],[],lb,ub,[],options); % defines a set of lower and upper bounds on the design variables in x, 
%                                           % so that the solution is always in the range lb ≤ x ≤ ub. If no equalities 
%                                           % exist, set Aeq = [] and beq = []. If x(i) is unbounded below, set 
%                                           % lb(i) = -Inf, and if x(i) is unbounded above, set ub(i) = Inf.
    % #####################################################################

    % globalsearch call **
    % #####################################################################
%    %Start Parallel Pool
%     poolobj = gcp;
%     files2add = listFilesInDir({"/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen","/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/functions","/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/FreqAnalysis"});
%     % Add Attached Files to Parallel Pool
%     addAttachedFiles(poolobj,files2add);

% 
%     %gs = GlobalSearch('XTolerance',0.0001,'NumTrialPoints',20,'NumStageOnePoints',20,'MaxTime', 1200);
%     %gs = GlobalSearch('XTolerance', 0.000000000000001,'NumTrialPoints',1000,'NumStageOnePoints',1000,'MaxTime', 600);
%     gs = GlobalSearch('NumTrialPoints',1000,'NumStageOnePoints',1000);
%     
%     %opts = optimoptions(@fmincon,'FiniteDifferenceStepSize',stpSiz,'MaxIterations',150);
% 
%     %opts = optimoptions(@fmincon,'OptimalityTolerance', 0.000000000000001,'StepTolerance',0.00000000000000000001,'ConstraintTolerance',0.000000000000001,'MaxFunctionEvaluations', 10000, 'MaxIterations',10000,'FiniteDifferenceType','central','UseParallel',false);
%     opts = optimoptions(@fmincon,'FiniteDifferenceType','central','UseParallel',false);
%     
%     %opts = optimoptions(@fmincon,'UseParallel',false);
%     %opts
%     problem = createOptimProblem('fmincon','x0',optParsStart,...
%     'objective',@ImgSFreqComp_v6,'lb',lb,'ub',ub,options=opts);
%     [x,fval,exitflag,output,solutions] = run(gs,problem);
    % #####################################################################

    % multistart call
    % #####################################################################
%     ms = MultiStart;
%     opts = optimoptions(@fmincon,'FiniteDifferenceStepSize',stpSiz,'MaxIterations',150);
%     opts
%     problem = createOptimProblem('fmincon','x0',optParsStart,...
%     'objective',@ImgSFreqComp_v5,'lb',lb,'ub',ub,options=opts);
%     [x,fval,exitflag,outpt,solutions] = run(ms,problem,10);
    % #####################################################################
    
    % patternsearch call
    % #####################################################################
    %opts = optimoptions(@patternsearch,'MaxFunEvals',150);
    opts = optimoptions(@patternsearch,"TolMesh", 1e-14,'FunctionTolerance',1e-14,'StepTolerance',1e-14,'MaxFunEvals',150);
    %opts = optimoptions(@patternsearch,'MaxFunEvals',150);

    opts
    [x,fval,exitflag,output] = patternsearch(@ImgSFreqComp_v7,optParsStart,[],[],[],[],lb,ub,[],opts);
    
    % #####################################################################
    
    % surrogateopt call
    % #####################################################################
%     opts = optimoptions(@surrogateopt,'MaxFunctionEvaluations',10,'PlotFcn',[],'OutputFcn',[]);
%     %opts
%     %[x,fval,exitflag,output] = surrogateopt(@ImgSFreqComp_v6,lb,ub,[1],[],[],[],[],opts);
%     [x,fval,exitflag,output] = surrogateopt(@ImgSFreqComp_v7,lb,ub,[],[],[],[],[],opts);
    
    
    %x = surrogateopt(objconstr,lb,ub,intcon,A,b,Aeq,beq)
    % #####################################################################
    this_output = {x,fval,exitflag,output,optParsStart};
    outmat{end+1} = this_output;
    
    end
    toc
    solutions_mat = zeros(2,size(outmat,2),length(outmat{1,1}{1,1}));
    
    for jj = 1:size(outmat,2)
        solutions_mat(1,jj,:) = outmat{1,jj}{1,1};
        solutions_mat(2,jj,:) = outmat{1,jj}{1,2};
    end
    clear jj;
    
    for jj=1:size(solutions_mat,3)
    figure
    hold
    titlestr=strcat("sig",num2str(jj));
    title(titlestr);
    scatter(solutions_mat(1,:,jj),solutions_mat(2,:,jj))
    end
    clear jj

    for jj=1:size(solutions_mat,3)
    figure
    histogram(solutions_mat(1,:,jj))
%     mean_sol = mean(solutions_mat(1));
%     min_sol = min(solutions_mat(1));
%     max_sol = max(solutions_mat(1));
%     pause = "";
    end
end