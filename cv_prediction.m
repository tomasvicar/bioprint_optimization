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




prediction = zeros(size(label));
prediction_sigma = zeros(size(label));
for cv_it = 1:size(featrues,1)

InitialX = featrues;
InitialObjective = score_abs;

InitialX(cv_it,:) = [];
InitialObjective(cv_it) = [];

    
bayesianOptimization = bayesopt_custom(@sum,[C_,A_,G_,Temperature_ink_,pressure,Speed],...
    'InitialX',InitialX,'InitialObjective',...
    InitialObjective,'MaxObjectiveEvaluations',0,'NumSeedPoints',0);


[objective,sigma] = predictObjective(bayesianOptimization,featrues(cv_it,:));
drawnow;

prediction(cv_it) = objective;
prediction_sigma(cv_it) = sigma;


close all force;
end



final_table = data;


dif_scoreabs_prediction = score_abs-prediction;

final_table = addvars(final_table,score_abs, prediction,dif_scoreabs_prediction, prediction_sigma);

writetable(final_table,'predictions_crossval.xlsx')

