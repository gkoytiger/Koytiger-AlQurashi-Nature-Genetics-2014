%function [ymut, phospho_ymut, unique_phospho_ymut] = import_cosmic(release, Phosphosites)
%Greg Koytiger August 2012
%Imports COSMIC mutations and keeps only those that are mutations that modify
%the sequences around tyrosine
%Returns file in format:
%GENE_NAME pY ORIGINAL_AA MUTATED_RES NEW_AA CANCER_TYPE TUMOR_ID WT_SEQ MUT_SEQ

load('Output/PhosphositeData.mat');
release = 'CosmicMutantExport';
data_path = 'Input/';
fid = fopen([data_path release '.tsv']);

% if fid==-1 %%If COSMIC release not local, then download & untar it
%     disp('Downloading COSMIC data');
%     gunzip(['ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/' release '.tsv.gz'], data_path);
%     fid = fopen([data_path release '.tsv']);
%     %disp('Downloading gene data');
%     %untar('ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta.tgz', data_path);
% end
%     
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'Delimiter', '\t'); %loads all of the fileds in COSMIC data fileC = C{:};   %Converts to cell of strings
C=[C{:}];
headers=C(1,:);
C(1,:)=[];
fclose(fid);

%%
GENE_NAME=strcmp(headers,'Gene name');
SAMPLE_NAME=strcmp(headers,'Sample name');
PRIMARY_SITE=strcmp(headers,'Primary site');
PRIMARY_HISTOLOGY=strcmp(headers,'Primary histology');
GENOME_WIDE_SCREEN=strcmp(headers,'Genome-wide screen');
MUTATION_ID=strcmp(headers,'Mutation ID');
MUTATION_AA=strcmp(headers,'Mutation AA');
MUTATION_DESCRIPTION=strcmp(headers, 'Mutation Description');
MUTATION_POSITION = strcmp(headers, 'Mutation genome position');
MUTATION_STRAND = strcmp(headers, 'Mutation strand');
eliminate = ~strcmp(C(:,MUTATION_DESCRIPTION), 'Substitution - Missense');
C(eliminate,:)=[];

%%
ymut = cell(size(C,1), 14);
gene = ''; k = 1;

%for i = 1:1000
for i = 1:size(C,1)

    str = C{i,MUTATION_AA}(3:end);

    if(isempty(regexp(str, '[a-z,*,?,>]', 'once')) && ~isempty(str))        %Extracts only data that concerns point mutations
        oAA = str(1);                            %Original AA
        nAA = str(end);
        mutloc = sscanf(str(2:(end-1)), '%i');
        
        
        if ~strcmp(oAA, nAA)                  %If mutation is not silent
            try
                %Reads new fasta file if needed
                
                if ~(strcmp(gene, C{i, GENE_NAME}))
                    gene = C{i, GENE_NAME};
                    [~,seq] = fastaread([data_path 'fasta/',gene(1), '/',gene,'_protein.txt']);
                end
                
                
                %Isolates the region around the mutation and prevents it from
                %referencing improper sequence areas
                
                pep_start = mutloc - 7;        pep_end = mutloc + 7;
                
                if pep_start < 1
                    pep_start = 1;
                end
                
                if pep_end > length(seq)
                    pep_end = length(seq);
                end
                
                yloc = strfind(seq(pep_start:pep_end), 'Y');                 %Finds if any Tyr are in viscinity of mutation
                
                if ~(isempty(yloc))
                    
                    mut_seq = seq;
                    mut_seq(mutloc) = nAA;
                    yloc = mutloc - 8 + yloc;
                    
                    for j = 1:length(yloc)
                        
                        if(yloc(j) - 7 > 0 && yloc(j) + 7 < length(seq))     %Only take "full" peptides
                            
                            wt_pep = seq(yloc(j)-7:yloc(j)+7);
                            mut_pep = mut_seq(yloc(j)-7:yloc(j)+7);
                            
                            
                            if(strcmp(mut_pep(8), 'Y'))                      %"phosphorylate" tyrosine
                                mut_pep(8) = 'y';
                            end
                            
                            wt_pep(8) = 'y';
                            
                            ymut{k,1} = gene;
                            ymut{k,2} = yloc(j);
                            ymut{k,3} = oAA;
                            ymut{k,4} = mutloc-yloc(j)+8;
                            ymut{k,5} = nAA;
                            ymut{k,6} = wt_pep;
                            ymut{k,7} = mut_pep;
                            ymut{k,8} = C{i,GENOME_WIDE_SCREEN};
                            ymut{k,9} = C{i, PRIMARY_SITE};
                            ymut{k,10} = C{i, PRIMARY_HISTOLOGY};
                            ymut{k,11} = C{i, SAMPLE_NAME};
                            ymut{k,12}= C{i, MUTATION_ID};
                            ymut{k,13}= C{i,MUTATION_POSITION};
                            ymut{k,14}= C{i,MUTATION_STRAND};
                            k = k+1;
                        end
                        
                    end
                end
                
            catch err
                
            end
        end
    end
    
end

%%
ymut(k:end,:) = [];
snp_sites = cell(size(ymut,1),1);
p_sites = cell(size(Phosphosites,1),1);


for i = 1:size(ymut,1)
    snp_sites{i} = [ymut{i,1} '_' int2str(ymut{i,2})];
end

for i = 1:size(Phosphosites,1)
    p_sites{i} = [Phosphosites{i,1} '_' int2str(Phosphosites{i,3})];
end

phospho_ymut = ymut(ismember(snp_sites, p_sites),:);

mutations = cell(size(phospho_ymut,1),1);

for i = 1:size(phospho_ymut,1)
    mutations{i} = [phospho_ymut{i,1} '_' int2str(phospho_ymut{i,2}) '_' phospho_ymut{i,3} int2str(phospho_ymut{i,4}) phospho_ymut{i,5}];
end

[~, m] = unique(mutations);
isRedundant = ~ismember(1:size(phospho_ymut,1), m);
unique_phospho_ymut = phospho_ymut(~isRedundant,:);

save('Output/COSMIC_Data.mat', 'ymut', 'unique_phospho_ymut', 'phospho_ymut')

clear ans fid fileToRead gene_name i k mutloc nAA oAA str err homolog humanid orthoseq start tmp unique_loc

%%
%%Import SH2 Domain mutants
load('Input/SH2_Domains.mat');

Domain_idx = zeros(size(SH2_Domains,1),2);

%Converts COSMIC indices to Uniprot domain indices, different offsets due
%to different isoforms being canonical, affects small fraction of SH2s

for i = 1:size(SH2_Domains,1)
    try
        [~,seq] = fastaread([data_path 'fasta/',SH2_Domains{i,1}(1), '/',SH2_Domains{i,1},'_protein.txt']);
        [~, ~, starts] = swalign(SH2_Domains{i,4}, seq);
        Domain_idx(i,1) = starts(2);
        Domain_idx(i,2) = starts(2)+length(SH2_Domains{i,4})-1;
    catch err
        Domain_idx(i,1) = 0;
        Domain_idx(i,2) = 0;
    end
end
%%
isSH2 = ismember(C(:,1), SH2_Domains(:,1));
SH2_mutants = C(isSH2,:);
Domain_mutants= cell(size(SH2_mutants,1),14);

for i = 1:size(SH2_mutants,1)
    str = SH2_mutants{i,MUTATION_AA}(3:end);
    if(isempty(regexp(str, '[a-z,*,?,>]', 'once')) && ~isempty(str)) 
        oAA = str(1);
        nAA = str(end);
        mutloc = sscanf(str(2:(end-1)), '%i');
        
        domain_hit = mutloc >= Domain_idx(:,1) & mutloc < Domain_idx(:,2) & strcmp(SH2_mutants(i,1), SH2_Domains(:,1));
        
        if sum(domain_hit) == 1 && ~strcmp(oAA,nAA)
            Domain_mutants{i,1} = SH2_Domains{domain_hit,1};
            Domain_mutants{i,2} = SH2_Domains{domain_hit,2};
            Domain_mutants{i,3} = oAA;
            Domain_mutants{i,4} = mutloc - Domain_idx(domain_hit,1) + 1; %Converts numbering to domain from protein
            Domain_mutants{i,5} = nAA;
            Domain_mutants{i,8} = SH2_mutants{i,GENOME_WIDE_SCREEN};
            Domain_mutants{i,9} = SH2_mutants{i, PRIMARY_SITE};
            Domain_mutants{i,10} = SH2_mutants{i, PRIMARY_HISTOLOGY};
            Domain_mutants{i,11} = SH2_mutants{i, SAMPLE_NAME};
            Domain_mutants{i,12}= SH2_mutants{i, MUTATION_ID};
            Domain_mutants{i,13}= SH2_mutants{i,MUTATION_POSITION};
            Domain_mutants{i,14}= SH2_mutants{i,MUTATION_STRAND};
        end
    end
        
end
Domain_mutants(cellfun('isempty', Domain_mutants(:,2)),:) = [];

%%
%Maps domain mutations onto SCOP alignment of parent domain

for i = 1:size(Domain_mutants,1)
    idx = strcmp(Domain_mutants{i,2}, SH2_Domains(:,2));   
    hmm_idx   = find(SH2_Domains{idx,8} == Domain_mutants{i,4});
    
    if ~isempty(hmm_idx)
        Domain_mutants{i,6} = SH2_Domains{idx,5};
        Domain_mutants{i,7} = SH2_Domains{idx,5};
        Domain_mutants{i,7}(hmm_idx)= Domain_mutants{i,5}; %Mutate HMM aligned seq
        Domain_mutants{i,4} = hmm_idx;
    end
    
end

Domain_mutants(cellfun('isempty', Domain_mutants(:,6)),:) = [];
save('Output/COSMIC_SH2_mutants', 'Domain_mutants');