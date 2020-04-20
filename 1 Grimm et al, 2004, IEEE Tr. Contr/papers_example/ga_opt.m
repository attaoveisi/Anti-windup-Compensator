function [x,fval,exitflag,output,population,score] = ga_opt(nvars,lb,ub,PopInitRange_Data,PopulationSize_Data,CrossoverFraction_Data,MigrationInterval_Data,MigrationFraction_Data,Generations_Data,TolFun_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'CrossoverFraction', CrossoverFraction_Data);
options = gaoptimset(options,'MigrationInterval', MigrationInterval_Data);
options = gaoptimset(options,'MigrationFraction', MigrationFraction_Data);
options = gaoptimset(options,'Generations', Generations_Data);
options = gaoptimset(options,'TolFun', TolFun_Data);
options = gaoptimset(options,'SelectionFcn', @selectionuniform);
options = gaoptimset(options,'MutationFcn', @mutationadaptfeasible);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv });
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'UseParallel', 1 );
[x,fval,exitflag,output,population,score] = ...
ga(@find_N,nvars,[],[],[],[],lb,ub,[],[],options);
