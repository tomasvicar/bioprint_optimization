clc;clear all; close all;

table_name = 'test.xlsx';
rng(42)

copyfile(table_name,['backup/' replace(table_name,'.xlsx','') replace(datestr(now),':','_')  '.xlsx'])




data = readtable(table_name);


featrues = removevars(data,{'Exp_No','FibreFormation','LayerStacking','score_abs','prediction','prediction_sigma'});


score_abs = data{:,'score_abs'};




C_ = optimizableVariable('C_',[0,90],'Type','real');
A_ = optimizableVariable('A_',[0,55],'Type','real');
G_ = optimizableVariable('G_',[0,86],'Type','real');
Temperature_ink_ = optimizableVariable('Temperature_ink_',[24,130],'Type','real');
pressure = optimizableVariable('pressure',[5,300],'Type','real');
Speed = optimizableVariable('Speed',[5,20],'Type','real');

vars = [C_,A_,G_,Temperature_ink_,pressure,Speed];

InitialX = featrues;
InitialObjective = score_abs;

if any(isnan(InitialObjective))
    error('please add score')
end

if isempty(InitialX)
    
    for kk = 1:3
        rand_vars_value = {};
        for k = 1:length(vars)
            rand_var_value = vars(k).Range(1) + rand*vars(k).Range(2)- vars(k).Range(1);
            rand_vars_value = [rand_vars_value,rand_var_value];
        end
        rand_vars_value = [1,rand_vars_value,nan,nan,nan,nan,nan];
    
        data = [data;rand_vars_value];
        
    end
    
else
    bayesianOptimization = bayesopt_custom(@sum,vars,...
        'InitialX',InitialX,'InitialObjective',...
        InitialObjective,'MaxObjectiveEvaluations',0,'NumSeedPoints',0);
    
    nextPoint = bayesianOptimization.NextPoint;
    [objective,sigma] = predictObjective(bayesianOptimization,nextPoint);
    
    Exp_No = size(data,1)+1;
    FibreFormation = nan;
    LayerStacking = nan;
    score_abs = nan;
    prediction = objective;
    prediction_sigma = sigma;
    nextPoint = addvars(nextPoint,Exp_No,'Before','C_');
    nextPoint = addvars(nextPoint,FibreFormation,LayerStacking,score_abs,prediction,prediction_sigma);
    
    data = [data;nextPoint];
    
    disp(nextPoint)
end



writetable(data,table_name)

