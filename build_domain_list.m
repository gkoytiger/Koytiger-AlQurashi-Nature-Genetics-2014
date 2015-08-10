load('Input/SCOP_SH2.mat');
fid = fopen('Input/uniprot_SH2.tab');
D = textscan(fid,'%s %*s %*s %*s %s %*s %*s %*s', 'delimiter', '\t',  'HeaderLines', 1);
fclose(fid);

SH2_Proteins = cell(size(D{1},1),7);

for i = 1:size(D{1,2},1)
    SH2_Proteins{i,1} = textscan(D{1,2}{i},  '%s %*[^\n]');
    SH2_Proteins{i,1} = SH2_Proteins{i,1}{1}{1};
    SH2_Proteins{i,2} = D{1}{i};
    [~, SH2_Proteins{i,7}] = fastaread(['http://www.uniprot.org/uniprot/' D{1}{i} '.fasta']);
    B = textscan(urlread(['http://www.uniprot.org/uniprot/' D{1}{i} '.gff']), '%*s %*s %*s %u %u %*s %*s %*s %s', 'delimiter', '\t');
    loc = find(strcmp(B{1,3}, 'Note=SH2') | strcmp(B{1,3}, 'Note=SH2%3B atypical') | strcmp(B{1,3}, 'Note=SH2-like') | strcmp(B{1,3}, 'Note=SH2 1') |strcmp(B{1,3}, 'Note=SH2 2')  );
    SH2_Proteins{i,3} = B{1}(loc(1));
    SH2_Proteins{i,4} = B{2}(loc(1));
    if length(loc) == 2
        SH2_Proteins{i,5} = B{1}(loc(2));
        SH2_Proteins{i,6} = B{2}(loc(2));
    end
end

%%
[~,idx] = sort(SH2_Proteins(:,1));
SH2_Proteins = SH2_Proteins(idx,:);

SH2_Domains = cell(sum(~cellfun('isempty', [SH2_Proteins(:,4); SH2_Proteins(:,6)])),8);
k = 1;

for i = 1:length(SH2_Proteins)
    SH2_Domains{k,1} = SH2_Proteins{i,1};
    SH2_Domains{k,2} = SH2_Proteins{i,1};
    SH2_Domains{k,3} = SH2_Proteins{i,2};
    
    dstart = SH2_Proteins{i,3} - 10;
    dend = SH2_Proteins{i,4} + 10;
    
    if dstart < 1
        dstart = 1;
    end
    
    if dend > length(SH2_Proteins{i,7})
        dend = length(SH2_Proteins{i,7});
    end
    
    dstart = uint32(dstart); dend=uint32(dend);
    SH2_Domains{k,4} = SH2_Proteins{i,7}(dstart:dend);
    [~,str,ptr] = hmmprofalign(SCOP_SH2, SH2_Domains{k,4});
    s1 = regexp(str, '[a-z]'); str(s1) = [];
    SH2_Domains{k,5} = str;
    SH2_Domains{k,6}=dstart;
    SH2_Domains{k,7}=dend;
    SH2_Domains{k,8}=ptr;
    
    if ~isempty(SH2_Proteins{i,5})
        SH2_Domains{k,2} = [SH2_Proteins{i,1} '-N'];
        k = k + 1;
        SH2_Domains{k,1} = SH2_Proteins{i,1};
        SH2_Domains{k,2} = [SH2_Proteins{i,1} '-C'];
        SH2_Domains{k,3} = SH2_Proteins{i,2};
        
        dstart = SH2_Proteins{i,5} - 10;
        dend = SH2_Proteins{i,6} + 10;
        
        if dstart < 1
            dstart = 1;
        end
        
        if dend > length(SH2_Proteins{i,7})
            dend = length(SH2_Proteins{i,7});
        end
        
        dstart = uint32(dstart); dend=uint32(dend);

        SH2_Domains{k,4} = SH2_Proteins{i,7}(dstart:dend);
        [~,str,ptr] = hmmprofalign(SCOP_SH2, SH2_Domains{k,4});
        s1 = regexp(str, '[a-z]'); str(s1) = [];
        SH2_Domains{k,5} = str;
        SH2_Domains{k,6}=dstart;
        SH2_Domains{k,7}=dend;
        SH2_Domains{k,8}=ptr;
        
    end
    
    k = k + 1;
end
%%
%Fixes pointers so that they point to aligned residues (Random matlab bug I think)        
for i = 1:122
    pointer_seq = '----------------------------------------------------------------------------------------------------------';
    clean_hmm = SH2_Domains{i,5};
    fixed_pointers = SH2_Domains{i,8};
    
    for j = 1:106
        if(~isnan(SH2_Domains{i,8}(j)))
            pointer_seq(j) = SH2_Domains{i,4}(SH2_Domains{i,8}(j));
            if(~(pointer_seq(j) == clean_hmm(j)))
                fixed_pointers(j) = fixed_pointers(j-1) + 1;
                pointer_seq(j) = SH2_Domains{i,4}(fixed_pointers(j));
            end
            
        end
    end
    SH2_Domains{i,8} = fixed_pointers;
    
end
save('Input/SH2_Domains.mat', 'SH2_Domains');

