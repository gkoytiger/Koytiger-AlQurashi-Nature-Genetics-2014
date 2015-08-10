load('Output/TrainingData.mat', 'peptides_cesarini', 'peptides_nash', 'peptides_macbeath', 'peptides_ltp', 'peptides_jones');
load('Input/SH2_Domains_wNC.mat');
path='Input/othermodels/SH2PepInt/';

peptides = cat(1, peptides_cesarini, peptides_jones, peptides_ltp, peptides_macbeath, peptides_nash);
%%
%Generates input file for model
forprediction = cell(length(peptides),1);
header = cell(length(peptides),1);

for i = 1:length(peptides)
    forprediction{i} = upper(peptides{i}(6:12));
    header{i}=num2str(i);
end

fastawrite([path 'peptides.fasta'], header, forprediction);
%%

models = importdata([path 'models.txt']);

map = {'APS','SH2B2';
    'BRDG1', 'STAP1';
    'CTEN', 'TNS4';
    'E105251', 'SHD';
    'E109111', 'SUPT6H';
    'E185634', 'SHC4';
    'EAT2', 'SH2D1B';
    'MIST', 'CLNK';
    'SH2B', 'SH2B1';
    'TENS1', 'TNS3';
    'TNS', 'TNS1'};



%%
predictions = importdata([path 'predictions.txt']); %Load predictions after running model

for i = 1:length(models)
    [tf,loc]=ismember(models{i}, map);
    if(tf)
        models{i}=map{loc,2};
    end
end   

[~,sh2idx] = ismember(models, SH2_Domains(:,2));
i = 1; k = 1;
probs= -1* ones(length(SH2_Domains), length(peptides)); %Set default to unanalyzed

while i < length(predictions.textdata)

    probs(sh2idx(1),k) = 2^(predictions.data(i)); %Raise to power of 2 to make everything positive, wont change AUC
    
    for j = 2:length(models)
        probs(sh2idx(j),k) = 2^str2double((predictions.textdata{i+j-1,2}));
    end
    
    k= k + 1;
    i = i + j ;
end

save([path 'sh2pepint'], 'probs')
xlswrite([path 'sh2pepint.xlsx'], probs)