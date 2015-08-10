%Greg Koytiger Sept 2013
%Process Bead Images of GCSAM (GCET2) Binding

im_path = '..\RawData\Images\GCSAM\Beads\';
L = dir([im_path '*.tiff']);
filenames = {L(:).name};

rows = zeros(numel(filenames),1);
columns = rows; fields = rows; channels  = rows;

for i = 1:numel(filenames)
    rows(i) = sscanf(filenames{i}(2:3),'%u');
    columns(i) = sscanf(filenames{i}(5:6),'%u');
    fields(i) = sscanf(filenames{i}(8:9),'%u');
    channels(i) = sscanf(filenames{i}(15),'%u');
end

data = zeros(length(unique(rows)), length(unique(columns)), length(unique(fields)), length(unique(channels)));
r0 = min(rows); c0 = min(columns);

for i = 1:numel(filenames)
        I = imread([im_path filenames{i}]);
        data(rows(i)+1-r0, columns(i)+1-c0, fields(i), channels(i)) = median(I(:));
end

unique_chan = unique(channels);
for i = 1:numel(unique_chan)
    assignin('base', ['channel' num2str(unique_chan(i))], mean(data(:,:,:,i),3));
    %assignin('base', ['channel' num2str(unique_chan(i)) '_e'], std(data(:,:,:,i),0,3));
end

clear im_path L filenames rows columns channels fields r0 c0 I i unique_chan

%Analyze stimulated cells
sh2_probed = {'GRB2',	'LCK',	'HCK',	'YES1',	'VAV1',	'BLK',	'FGR',	'TENC1',	'CRK',	'LYN',	'ABL1',	'CRKL'};
pY_protein = {'GCSAM'};

%Normalized signal is [BeadRFP]/([BeadGFP] * [SupRFP])

%Replicate1
replicate1 = zeros(12,1);

for i = 1:12
    j = i*2-1;
    replicate1(i) = channel2(9,j)/(channel1(9,j) * channel2(3,j));
    
end

%Replicate2
replicate2 = zeros(12,1);


for i = 1:12
    j = i*2;
    replicate2(i) = channel2(9,j)/(channel1(9,j) * channel2(3,j));
end

[cf1,gof1] = fit(replicate1,replicate2,'poly1');

minsignal  = min([replicate1;replicate2]);
f1 = replicate1./min(replicate1); f1 = f1./(1+f1);
f2 = replicate2./min(replicate2); f2 = f2./(1+f2);
error = std([f1,f2],1,2);
fmean = mean([f1, f2], 2);
%%
SH2=importdata('Output/SH2s.dat');
proteins = importdata('Output/Phosphoproteins.dat');
prob = importdata('Output/WT Predictions.dat');
[~,sh2idx] = ismember(sh2_probed,SH2); modelp = prob(sh2idx, strcmp(proteins, pY_protein));


%%
[cf,gof] = fit(modelp,fmean,'poly1', 'Robust', 'Bisquare');
figure;errorbar(modelp,fmean,error,'.', 'LineWidth',2, 'MarkerSize',15);  axis([.2 1 .5 .91]); axis square; hold on;  p= plot(cf); set(p, 'LineWidth', 2); xlabel('Model Probability'); ylabel('F/(1+F)'); legend off;
text(modelp,fmean,sh2_probed);hold off; 
print('-depsc2', 'Output/figures/ModelP_vs_Exp');
%xlswrite('Output/figures/ModelP_vs_Exp');

figure;hold on; plot(replicate1,replicate2,'.', 'LineWidth',2, 'MarkerSize',15); axis([0 9*10^-5 0 8*10^-5]);  p= plot(cf1); set(p, 'LineWidth', 2); legend off; axis square;
text(replicate1,replicate2,sh2_probed);hold off; 
print('-depsc2', 'Output/figures/Exp_vs_Exp');

%%
%No stim
nostim = zeros(12,1);
for i = 1:12
    j = i*2-1;
    nostim(i) = channel2(10,j)/(channel1(10,j) * channel2(4,j));
end

meansignal = mean([replicate1,replicate2],2); [~, idx] = sort(meansignal);
bar([meansignal(idx),nostim(idx)],1.25); set(gca,'YLim',[min(nostim)*.9 max(meansignal)*1.1],'XTickLabel',sh2_probed(idx));
print('-depsc2', 'Output/figures/stim_vs_nostim');


blank_gfp = channel1([2,4,6,8], [2,4,6,8,10,12,14,16,18,20,22,24]);
blank_rfp = channel2([2,4,6,8], [2,4,6,8,10,12,14,16,18,20,22,24]);