%% test bed for this function
% nb_test = findNeighbors(A_or, d1, d2, 25);
% a=zeros(d1, d2);
% for i=1:size(nb, 1)
%     nb_test{i}
%     a(:) = 0;
%     mask = find(A_or(:, i));
%     a(mask) = i;
%     figure(101), hold off
%     plot(C_or(i, :))
%     for j=1:length(nb{i})
%         mask = find(A_or(:,nb{i}(j)));
%         a(mask) = nb_test{i}(j);
%         figure(101), hold on, plot(C_or(nb_test{i}(j), :));
%     end
%     figure(1), imagesc(a), axis equal, axis tight
%     pause
% end

%%
function [neighborslist, areas] = findNeighbors(spatialCompents, width, height, thresh)
%findNeighbors For each spatial component, find a list of neighboring
% components, as defined by closest spatial distance between two components
% when it is less than "thresh"
%   

nC = size(spatialCompents, 2);
neighborslist = cell(nC, 1);
areas = zeros(1, nC);

for c1 = 1 : nC
    maskLinearIdx = find(spatialCompents(:, c1));
    areas(c1) = numel(maskLinearIdx);

    [y1, x1] = ind2sub([height width], maskLinearIdx);
    for c2 = c1+1 : nC
%     for c2 = cat(2, 1:(c1-1), c1+1:nC)
        maskLinearIdx = find(spatialCompents(:, c2));
        [y2, x2] = ind2sub([height width], maskLinearIdx);
        distanceList = distancesBetweenComponents(x1, y1, x2, y2); % implemented in C-mex
        if min(distanceList) < thresh
            neighborslist{c1} = [neighborslist{c1} c2];
            neighborslist{c2} = [neighborslist{c2} c1];
        end
        
    end
end
end