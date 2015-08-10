%Greg Koytiger
%Aug 2013
%Runs all MATLAB code

load('Input/SH2_Domains_wNC.mat');

%% Imports all training data for modeling
%2=strong bind 1=weak bind 0=not bind -1=not test

[peptides_cesarini, binds_cesarini, quant_cesarini] = import_Cesarini2013_HTP(SH2_Domains);
[peptides_ltp, binds_ltp] = import_Cesarini2013_LTP(SH2_Domains);
[peptides_jones, binds_jones, quant_jones] = import_Jones2012(SH2_Domains);
[peptides_macbeath, binds_macbeath, quant_macbeath] = import_MacBeath2013(SH2_Domains);
[peptides_nash, binds_nash] = import_Nash2012(SH2_Domains);

save('Output/TrainingData.mat'); clear all;

