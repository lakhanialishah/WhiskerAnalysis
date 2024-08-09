function groups = separating_units_into_conditions(keep4, IDlist)

% L2/3/4 combined, L5/6
L2cwd0_16 = []; L2cwd0_21 = []; L2cwd2 = []; L2cwd7 = [];
L5cwd0_16 = []; L5cwd0_21 = []; L5cwd2 = []; L5cwd7 = [];
L2kwd0_16 = []; L2kwd0_21 = []; L2kwd2 = []; L2kwd7 = [];
L5kwd0_16 = []; L5kwd0_21 = []; L5kwd2 = []; L5kwd7 = [];

L2FScwd0_16 = []; L2FScwd0_21 = []; L2FScwd2 = []; L2FScwd7 = [];
L5FScwd0_16 = []; L5FScwd0_21 = []; L5FScwd2 = []; L5FScwd7 = [];
L2FSkwd0_16 = []; L2FSkwd0_21 = []; L2FSkwd2 = []; L2FSkwd7 = [];
L5FSkwd0_16 = []; L5FSkwd0_21 = []; L5FSkwd2 = []; L5FSkwd7 = [];

groups = {};

for clu = 1:size(keep4,1)
    %code to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep4{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    
    %to combine control KO and WD into their own variables to get indices for RS units only
    %defining what kind of unit to take here
    if strcmp(IDlist{cluindex,2}.clu_info.classify, 'RS') == 1
        if strcmp(IDlist{cluindex,2}.clu_info.layer, 'L4') || strcmp(IDlist{cluindex,2}.clu_info.layer, 'L2/3') == 1
            if strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_16') == 1
                L2cwd0_16 = vertcat(L2cwd0_16, cluindex);
                groups.L2cwd0_16.IDlist_indices = L2cwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_21') == 1
                L2cwd0_21 = vertcat(L2cwd0_21, cluindex);
                groups.L2cwd0_21.IDlist_indices = L2cwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd2') == 1
                L2cwd2 = vertcat(L2cwd2, cluindex);
                groups.L2cwd2.IDlist_indices = L2cwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd7') == 1
                L2cwd7 = vertcat(L2cwd7, cluindex);
                groups.L2cwd7.IDlist_indices = L2cwd7;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_16') == 1
                L2kwd0_16 = vertcat(L2kwd0_16, cluindex);
                groups.L2kwd0_16.IDlist_indices = L2kwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_21') == 1
                L2kwd0_21 = vertcat(L2kwd0_21, cluindex);
                groups.L2kwd0_21.IDlist_indices = L2kwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd2') == 1
                L2kwd2 = vertcat(L2kwd2, cluindex);
                groups.L2kwd2.IDlist_indices = L2kwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd7') == 1
                L2kwd7 = vertcat(L2kwd7, cluindex);
                groups.L2kwd7.IDlist_indices = L2kwd7;
            end
        elseif strcmp(IDlist{cluindex,2}.clu_info.layer, 'L5/6') == 1
            if strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_16') == 1
                L5cwd0_16 = vertcat(L5cwd0_16, cluindex);
                groups.L5cwd0_16.IDlist_indices = L5cwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_21') == 1
                L5cwd0_21 = vertcat(L5cwd0_21, cluindex);
                groups.L5cwd0_21.IDlist_indices = L5cwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd2') == 1
                L5cwd2 = vertcat(L5cwd2, cluindex);
                groups.L5cwd2.IDlist_indices = L5cwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd7') == 1
                L5cwd7 = vertcat(L5cwd7, cluindex);
                groups.L5cwd7.IDlist_indices = L5cwd7;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_16') == 1
                L5kwd0_16 = vertcat(L5kwd0_16, cluindex);
                groups.L5kwd0_16.IDlist_indices = L5kwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_21') == 1
                L5kwd0_21 = vertcat(L5kwd0_21, cluindex);
                groups.L5kwd0_21.IDlist_indices = L5kwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd2') == 1
                L5kwd2 = vertcat(L5kwd2, cluindex);
                groups.L5kwd2.IDlist_indices = L5kwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd7') == 1
                L5kwd7 = vertcat(L5kwd7, cluindex);
                groups.L5kwd7.IDlist_indices = L5kwd7;
            end
        end
        %to categorize FS cells
    elseif strcmp(IDlist{cluindex,2}.clu_info.classify, 'FS') == 1
        if strcmp(IDlist{cluindex,2}.clu_info.layer, 'L4') || strcmp(IDlist{cluindex,2}.clu_info.layer, 'L2/3') == 1
            if strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_16') == 1
                L2FScwd0_16 = vertcat(L2FScwd0_16, cluindex);
                groups.L2FScwd0_16.IDlist_indices = L2FScwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_21') == 1
                L2FScwd0_21 = vertcat(L2FScwd0_21, cluindex);
                groups.L2FScwd0_21.IDlist_indices = L2FScwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd2') == 1
                L2FScwd2 = vertcat(L2FScwd2, cluindex);
                groups.L2FScwd2.IDlist_indices = L2FScwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd7') == 1
                L2FScwd7 = vertcat(L2FScwd7, cluindex);
                groups.L2FScwd7.IDlist_indices = L2FScwd7;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_16') == 1
                L2FSkwd0_16 = vertcat(L2FSkwd0_16, cluindex);
                groups.L2FSkwd0_16.IDlist_indices = L2FSkwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_21') == 1
                L2FSkwd0_21 = vertcat(L2FSkwd0_21, cluindex);
                groups.L2FSkwd0_21.IDlist_indices = L2FSkwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd2') == 1
                L2FSkwd2 = vertcat(L2FSkwd2, cluindex);
                groups.L2FSkwd2.IDlist_indices = L2FSkwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd7') == 1
                L2FSkwd7 = vertcat(L2FSkwd7, cluindex);
                groups.L2FSkwd7.IDlist_indices = L2FSkwd7;
            end
        elseif strcmp(IDlist{cluindex,2}.clu_info.layer, 'L5/6') == 1
            if strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_16') == 1
                L5FScwd0_16 = vertcat(L5FScwd0_16, cluindex);
                groups.L5FScwd0_16.IDlist_indices = L5FScwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd0_21') == 1
                L5FScwd0_21 = vertcat(L5FScwd0_21, cluindex);
                groups.L5FScwd0_21.IDlist_indices = L5FScwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd2') == 1
                L5FScwd2 = vertcat(L5FScwd2, cluindex);
                groups.L5FScwd2.IDlist_indices = L5FScwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'cwd7') == 1
                L5FScwd7 = vertcat(L5FScwd7, cluindex);
                groups.L5FScwd7.IDlist_indices = L5FScwd7;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_16') == 1
                L5FSkwd0_16 = vertcat(L5FSkwd0_16, cluindex);
                groups.L5FSkwd0_16.IDlist_indices = L5FSkwd0_16;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd0_21') == 1
                L5FSkwd0_21 = vertcat(L5FSkwd0_21, cluindex);
                groups.L5FSkwd0_21.IDlist_indices = L5FSkwd0_21;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd2') == 1
                L5FSkwd2 = vertcat(L5FSkwd2, cluindex);
                groups.L5FSkwd2.IDlist_indices = L5FSkwd2;
            elseif strcmp(IDlist{cluindex,2}.clu_info.exp_cond, 'kwd7') == 1
                L5FSkwd7 = vertcat(L5FSkwd7, cluindex);
                groups.L5FSkwd7.IDlist_indices = L5FSkwd7;
            end
        end
    end
end

end


