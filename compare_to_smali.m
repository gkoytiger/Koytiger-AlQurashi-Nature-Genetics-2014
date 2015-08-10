load('Output/TrainingData.mat', 'peptides_cesarini', 'peptides_nash', 'peptides_macbeath', 'peptides_ltp', 'peptides_jones');
load('Input/SH2_Domains_wNC.mat');
path='Input/othermodels/SMALI/';
L = dir([path '*sm.txt']);
filenames = {L(:).name};
peptides = cat(1, peptides_cesarini, peptides_jones, peptides_ltp, peptides_macbeath, peptides_nash);
probs= -1* ones(length(SH2_Domains), length(peptides)); %Set default to unanalyzed
positions=[6,7,9,10,11,12]; %Positions evaluated by PSSM
thresholds = importdata([path 'domains_thresholds.txt']); 

map= {'APS','SH2B2';
    'BRDG1','STAP1';
    'BRK','PTK6';
    'SH2B','SH2B1';
    'SHIP1','INPP5D';
    'SHIP2','INPPL1';
    'SLNK','SH2D6';
    'YES', 'YES1'};
    
for i = 1:length(filenames)
    pssm = importdata([path filenames{i}]);
    sh2= filenames{i}(1:strfind(filenames{i}, ' ')-1); %extract sh2 name
    thresh = thresholds.data(strcmp([sh2 ' SH2'], thresholds.textdata));
    sh2=strrep(sh2,'_','-');
    sh2_idx= strcmp(sh2,SH2_Domains(:,2));
    
    if sum(sh2_idx)==0
        sh2=map(strcmp(sh2,map(:,1)),2);
        sh2_idx= strcmp(sh2,SH2_Domains(:,2));
    end

    for k = 1:length(peptides)
        
            pep_prob= 0;
            for j = 1:length(positions)
                idx = strcmp(peptides{k}(positions(j)),pssm.colheaders);
                if ~sum(idx)==0
                    pep_prob=pep_prob + pssm.data(j,idx); %Add position scores
                end
            end
            
            probs(sh2_idx,k) = pep_prob/thresh;
        
        
    end
end

save('Output/SMALI_predictions.mat', 'SH2_Domains','peptides', 'probs');
        