load('Output/TrainingData.mat', 'peptides_cesarini', 'peptides_nash', 'peptides_macbeath', 'peptides_ltp', 'peptides_jones');
load('Input/SH2_Domains_wNC.mat');
path='Input/othermodels/Netphorest/';

peptides = cat(1, peptides_cesarini, peptides_jones, peptides_ltp, peptides_macbeath, peptides_nash);

%%
%Generates input file for model
forprediction = cell(length(peptides),1);
header = cell(length(peptides),1);

for i = 1:length(peptides)
    forprediction{i} = upper(peptides{i});
    header{i}=num2str(i);
end

fastawrite([path 'peptides.fasta'], header, forprediction);
