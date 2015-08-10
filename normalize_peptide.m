function [ out_pep ] = normalize_peptide( in_pep )
%GK Feb 2013
%Outputs standardized peptide sequence

pep_l = length(in_pep);
py_loc = find(in_pep == 'y');

if pep_l > 15
    out_pep = in_pep(py_loc - 7:py_loc+7);

else
    out_pep = '---------------';

    n_len = 9-py_loc;
    c_len = pep_l - py_loc+8;
    
    if c_len > 15
        in_pep(c_len-n_len+1:end) = []; %truncates long peptides
        c_len = 15;
    end

    out_pep(n_len:c_len) = in_pep;
end


end

