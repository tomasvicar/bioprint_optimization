clc; clear all; close all force;


data = readtable('../Amir-data2.xlsx');



data_clean = removevars(data,{'Exp_No','Temperature_platform_','FibreFormation','LayerStacking'});


featrues = removevars(data_clean,{'score'});
% features_names = featrues.Properties.VariableNames;
% featrues = featrues{:,:};


label = data_clean{:,'score'};

score_abs = abs(label);



C_ = optimizableVariable('C_',[0,90],'Type','real');
A_ = optimizableVariable('A_',[0,55],'Type','real');
G_ = optimizableVariable('G_',[0,86],'Type','real');
Temperature_ink_ = optimizableVariable('Temperature_ink_',[24,130],'Type','real');
pressure = optimizableVariable('pressure',[5,300],'Type','real');
Speed = optimizableVariable('Speed',[5,20],'Type','real');



InitialX = featrues;
InitialObjective = score_abs;

% InitialX(31,:) = [];
% InitialObjective(31) = [];

    

bayesianOptimization = bayesopt_custom(@sum,[C_,A_,G_,Temperature_ink_,pressure,Speed],...
    'InitialX',InitialX,'InitialObjective',...
    InitialObjective,'MaxObjectiveEvaluations',0,'NumSeedPoints',0);



[x,CriterionValue] = bestPoint(bayesianOptimization,'Criterion','min-upper-confidence-interval');

disp(x)
% disp(CriterionValue)

[objective,sigma] = predictObjective(bayesianOptimization,x);



disp(objective)
disp(sigma)
