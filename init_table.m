clc;clear all;close all;

table_name = 'test.xlsx';

var_names = {'Exp_No','C_','A_','G_','Temperature_ink_','pressure','Speed','FibreFormation','LayerStacking','score_abs', 'prediction', 'prediction_sigma'};


T = array2table(zeros(0,length(var_names)),'VariableNames',var_names);



writetable(T,table_name,'WriteMode','overwritesheet')
