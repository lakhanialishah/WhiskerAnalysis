function IDlist = assignFSvRS(keep4, IDlist)
% just using trough to peak 

all_tr2pks = [];
for clu = 1: length(keep4)
   where = strfind(IDlist(:,1), keep4{clu,1});
   findindex = find(~cell2mat(cellfun( @(x) isempty(x) , where , 'UniformOutput' , false ))); 
   all_tr2pks = vertcat(all_tr2pks, IDlist{findindex,2}.clu_info.tr2pk);
end
% figure
% histogram(all_tr2pks, 'BinWidth', .04)

gm = fitgmdist(all_tr2pks, 2);
% Extract means and variances
means = gm.mu;
variances = squeeze(gm.Sigma);
% Plot the histogram
figure %('visible', 'off')
histogram(all_tr2pks, 30, 'Normalization', 'pdf');
hold on;
% Plot the Gaussian distributions
x = linspace(min(all_tr2pks), max(all_tr2pks), 1000);
pdf1 = @(x) pdf('Normal', x, means(1), sqrt(variances(1)));
pdf2 = @(x) pdf('Normal', x, means(2), sqrt(variances(2)));
y1 = pdf1(x);
y2 = pdf2(x);
min_pdf = min(y1, y2);
overlap = trapz(x, min_pdf);
% Plot the Gaussian distributions
plot(x, y1, 'r-', 'LineWidth', 2);
plot(x, y2, 'b-', 'LineWidth', 2);
% Shade the overlapping area
fill([x, fliplr(x)], [min_pdf, zeros(size(min_pdf))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Length of unit');
ylabel('Probability Density');
title(['Histogram with Gaussian Mixture Model (Overlap Area: ', num2str(overlap), ')']);
hold off;

% to assign FS vs RS 
for clu = 1:size(keep4,1)
    %code to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep4{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    %based on the intersection of curves in overlap histogram is .72,
    %eliminating middle units is FS is all units below .63 ms and RS is all units above .802 ms
    if IDlist{cluindex,2}.clu_info.tr2pk < .72
        IDlist{cluindex,2}.clu_info.classify = 'FS';
    elseif IDlist{cluindex,2}.clu_info.tr2pk > .72
        IDlist{cluindex,2}.clu_info.classify = 'RS';
    else
        IDlist{cluindex,2}.clu_info.classify = 'notincluded';
    end
end

% to put all the FS waveforms above the RS waveforms 
FSindices = [];
RSindices = [];
for clu = 1:size(keep4,1)
    %code to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep4{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    if strcmp(IDlist{cluindex,2}.clu_info.classify, 'FS') == 1
        FSindices = vertcat(FSindices, cluindex);
    elseif strcmp(IDlist{cluindex,2}.clu_info.classify, 'RS') == 1 
        RSindices = vertcat(RSindices, cluindex);
    end
end
figure
hold on
clusters_to_delete = [];
for cluindex = 1:size(RSindices,1)
    troughindex = find(IDlist{RSindices(cluindex),2}.clu_info.mwf == min(IDlist{RSindices(cluindex),2}.clu_info.mwf));
    crosses0 = find(IDlist{RSindices(cluindex),2}.clu_info.mwf(troughindex:72) > 0, 1);
    normalizedmwf = IDlist{RSindices(cluindex),2}.clu_info.mwf/-IDlist{RSindices(cluindex),2}.clu_info.mwf(troughindex);
    indicestoplot = troughindex-12: troughindex+30;
    % to delete weird ones by eye manually, find cluster no then delete from IDlist
%     if normalizedmwf(indicestoplot(31)) > 0.8
%         clusters_to_delete = vertcat(clusters_to_delete, cluindex);
%     end
    patchline(1:43, normalizedmwf(indicestoplot), 'edgecolor', 'r', 'edgealpha',0.2)
end
for cluindex = 1:size(FSindices,1)
    troughindex = find(IDlist{FSindices(cluindex),2}.clu_info.mwf == min(IDlist{FSindices(cluindex),2}.clu_info.mwf));
    crosses0 = find(IDlist{FSindices(cluindex),2}.clu_info.mwf(troughindex:72) > 0, 1);
    normalizedmwf = IDlist{FSindices(cluindex),2}.clu_info.mwf/-IDlist{FSindices(cluindex),2}.clu_info.mwf(troughindex);
    indicestoplot = troughindex-12: troughindex+30;
    % to delete weird ones by eye manually, find cluster no then delete from IDlist
%     if normalizedmwf(indicestoplot(5)) < -0.5
%         clusters_to_delete = vertcat(clusters_to_delete, cluindex);
%     end
    patchline(1:43, normalizedmwf(indicestoplot), 'edgecolor', 'b', 'edgealpha',0.2)
end

end
