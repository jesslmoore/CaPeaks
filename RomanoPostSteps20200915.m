%% Pre-step: separate all flashing and non-flashing cells

%threshold = 0.5; 
threshold = 1; %usually 0.5 except for drug treatments/videos with higher background which I do 0.8 or 1 or 2

deltaf = bsxfun(@minus, cells_mean, min(cells_mean));
deltaf_f = bsxfun(@rdivide, deltaf, min(cells_mean));

% deltaf_ft = f_f0.'-0.95;
deltaf_ft = deltaf_f.';
f_f0t_flashing = deltaf_ft;
nonflashingcells = [];
flashingcells = [];

for i=1:cell_number
    if f_f0t_flashing(i,:) < threshold
        nonflashingcells((length(nonflashingcells)+1)) = i; %determining which cells don't flash above 1 fluorescence value
    else
        flashingcells((length(flashingcells)+1)) = i; %determining which cells flash above 1 fluorescence value
    end

end

f_f0t_flashing(nonflashingcells,:) = []; %creating matrix of fluorescence values of flashing cells
f_f0t_nonflashing = deltaf_ft(nonflashingcells,:);

%assigning cells pixel information to matrix
pixel_dimensions = size(avg);
d1 = pixel_dimensions(1,1);
d2 = pixel_dimensions(1,2);
pixel_total = d1*d2;
A_in_matrix=zeros(pixel_total,cell_number);

for i=1:cell_number
    A_in_matrix(cells{i,1},i) = 1;
end

A_in = sparse(A_in_matrix); %encoding spatial information of each ROI
A_or = A_in;
A_or(:,nonflashingcells) = []; % removing spatial information of non-flashing ROIs

%% Post-step 1

% f_f0t_flashing contains dF/F% traces of all flashing components
% # of components = cell_number
N = size(deltaf_ft, 1); % # of components
xdata = 1:size(deltaf_ft, 2); % timepoints
overallEvents = [];
fprintf(1, '\nProgress:     ')
for i=1:N
    if deltaf_ft(i,:) < threshold
        continue;
    end
    [pid, tid] = Find1TracePeaksTroughs(deltaf_ft(i,:), 0.2);
    if isempty(pid)
        continue;
    end
    [events, fitfun] = MergePeaks2FitGaussians(deltaf_ft(i,:), pid, tid);
    % Lin added following three lines 2020-08-31:
    if isempty(events)
        continue;
    end
    events(:,1) = i; % fill in the component info of these events
    overallEvents = [overallEvents; events];
       figure(100), hold off
       for j=1:size(fitfun, 1)
           plot(fitfun{j}, xdata', deltaf_ft(i,:));
           hold on
       end
       
       
       plot(pid, deltaf_ft(i, pid), 'x', 'MarkerSize', 8)
       tid_flat = unique(tid(:));
       plot(tid_flat, deltaf_ft(i, tid_flat), 'o', 'MarkerSize', 8)
       ax=gca; ax.Title.String = sprintf("Trace %d", i);
%        pause
    prog = 100 * i/N;
    fprintf(1,'\b\b\b\b%3.0f%%',prog)
end
fprintf(1, '\n')

%% Post-step 2

[nbList, areaList] = findNeighbors(A_in, d2, d1, 2);
% 2 pixels is the distance threshold to consider being neighbor; used to
% be 25 but this picked up components far away

updatedEvents = overallEvents; %overallEvents is from step1 and in the workspace; includes all events
corrThresh = 0.75;  % Threshold of correlation coefficient to be considered correlated enough
N_comp = numel(deltaf_ft(:,1)); % number of components, not events
% assign area to each event based on component/ROI it belongs to
for compI=1:N_comp
    updatedEvents(updatedEvents(:,1) == compI, 5) = areaList(compI);
end


%% Post-step 3

events_info = [];
T = numel(xdata);
for i=1:N
    events_info(i,1) = i; % Component ID
    events_info(i,2) = sum(updatedEvents(:,1) == i); % Number of peaks
    events_info(i,3) = areaList(1,i); % Component area
    events_info(i,4) = events_info(i,2)/(T*2); % Events per second
    events_info(i,5) = events_info(i,2)/(T/60*2); % Events per minute
end
N = size(deltaf_ft, 1); % # of components
for i=1:N
    if ismember(i,updatedEvents(:,1)) == 1 % need to skip cells that were thrown out for messy calcium traces
        events_info(i,6) = max(updatedEvents((updatedEvents(:,1)==i),3)); % max lingering time
    else
        continue
    end
end
avg_area = mean(areaList);
N = size(updatedEvents, 1); % update N to number of events


%% Post-step 4

% Threshold for start time difference. If event J's start time is after event I but 
% within threshT, then it's considered event I's downstream event.
threshT = 5; 
% If event J's start time is <threshTe after event I's end time, then it's
% also considered event I's downstream event.
threshTe = 5;

eventDownstream = cell(size(updatedEvents, 1), 1); 
nR = numel(deltaf_ft(:,1)); % nR is number of components
adjacencymatrix = zeros(cell_number);


for compI=1:nR
    eventsI_idx = find(updatedEvents(:,1) == compI);
    
    for compJ = nbList{compI}
       
       
        eventsJ_idx = find(updatedEvents(:,1) == compJ);
        % For every event of component I
        for evOfI = eventsI_idx'
            t0I = updatedEvents(evOfI, 2); %time to event
            t1I = updatedEvents(evOfI, 2) + updatedEvents(evOfI, 3); %time to event + lingering time of event
            for evOfJ = eventsJ_idx'
                t0Diff = updatedEvents(evOfJ, 2) - t0I; 
                t1Diff = updatedEvents(evOfJ, 2) - t1I;
                if t0Diff > -threshT && t0Diff < threshT || ...
                        t1Diff > -threshTe && t1Diff < threshTe
                    % Then consider evOfJ downstream neighbor of evOfI
                    eventDownstream{evOfI} = [eventDownstream{evOfI}, evOfJ];
                    adjacencymatrix(compI,compJ) = 1;
                    adjacencymatrix(compJ,compI) = 1;
                end
            end
        end
    end
end

%figure
%spy(adjacencymatrix)
%gives spy plot of adjacency matrix

%% plots and saves graph with and without node numbers
G = graph(adjacencymatrix);
numbers = 1:cell_number;
G.Nodes.Name = string(numbers.');
H = rmnode(G,nonflashingcells);
H.Nodes.Name = string(flashingcells.');

fig1 = figure('Position', [100 20 2540 1300]);
plot(H)
cutName=strfind(fileName,'.tif');
outputFile1=[fileName(1:cutName-1) '_graph2.tif'];
fig.PaperUnits = 'normalized';
fig.PaperPosition = [0 0 27 13];
print(fig1,'-dtiff',outputFile1,'-r300')

fig2 = figure('Position', [100 20 2540 1300]);
p = plot(H,'NodeLabel',H.Nodes.Name);
outputFile1=[fileName(1:cutName-1) '_graphnumbered2.tif'];
fig.PaperUnits = 'normalized';
fig.PaperPosition = [0 0 27 13];
print(fig2,'-dtiff',outputFile1,'-r300')

%% calculating neighborhood size based on graph
nR = numel(deltaf_ft(:,1)); % nR is number of components
bins = conncomp(G); % assigns each cell/node to a numbered bin based on its connections
bins2 = bins.'; 

updatedEvents(:,6)=zeros; 
for i=1:nR
    updatedEvents((updatedEvents(:,1)==i),6) = bins2(i,1); % assigns bin number to each event
end

nB = max(bins);
bins3 = zeros(nB,1); % going to create new variable  to quantify how many cells are in each bin/neighborhood
for var2=1:nB
    bins3(var2,1) = sum(bins(1,:)==var2); % for each bin up to nB, tells total size of bin
end

bins4 = zeros(nR,1); % going to create new variable that assigns size of neighborhood each cell is in
for var2=1:nR
    bins4(var2,1) = bins3((bins2(var2,1)),1);
    updatedEvents((updatedEvents(:,1)==var2),7) = bins4(var2,1); %counts number of nodes
end
nnfc = numel(nonflashingcells);
%assign non-flashing cells neighborhood size of 0
for i = 1:nnfc
    bins4(nonflashingcells(1,i),1) = 0;
end

updatedEvents(:,8) = zeros;
updatedEvents(:,9) = zeros;
updatedEvents(:,10) = zeros;
for i=1:nR
    updatedEvents((updatedEvents(:,1)==i),8) = max(updatedEvents((updatedEvents(:,1)==i),2))-min(updatedEvents((updatedEvents(:,1)==i),2)); % timesteps between first and last calcium events' starts
    updatedEvents((updatedEvents(:,1)==i),9) = events_info(i,2);
end
for i=1:N
    updatedEvents(i,10) = updatedEvents(i,9)/(updatedEvents(i,8)/60*2); % peaks per minute for timeperiod of activity
end

 %% calculate mean intensity of ROIs
% %properly format CC2 structure
% CC2.Connectivity = CC.Connectivity;
% CC2.ImageSize = CC.ImageSize;
% 
% %load fucci2 image
%  [filename2,pathname2] = uigetfile({'*.tif';'*.TIFF'},'Open fucci2 file', 'MultiSelect', 'off');
%  fileName2=fullfile(pathname2,filename2);
%  fucci2 = imread(fileName2);
%  
%  [rows, columns, numberOfColorChannels] = size(fileName2);
%  if numberOfColorChannels > 1
%  	promptMessage = sprintf('Your image file has %d color channels.\nThis script was designed for grayscale images.\nDo you want me to convert it to grayscale for you so you can continue?', numberOfColorChannels);
%  	button = questdlg(promptMessage, 'Continue', 'Convert and Continue', 'Cancel', 'Convert and Continue');
%  	if strcmp(button, 'Cancel')
%  		fprintf(1, 'Finished running.\n');
%  		return;
%  	end
%  	% Do the conversion using standard book formula
%  	fileName2 = rgb2gray(fileName2);
%  end
%  
%  MeanIntensities = cell2mat(transpose(struct2cell(regionprops(CC2,fucci2,'MeanIntensity'))));
%  MaxIntensities = cell2mat(transpose(struct2cell(regionprops(CC2,fucci2,'MaxIntensity'))));
%  
%  MaxIntensities = double(MaxIntensities);
%  
%  for n = 1:cell_number
%  
%      events_info(n,6) = MeanIntensities(n);
%      events_info(n,7) = MaxIntensities(n);
%      
%  end

%% code if I want to show fucci values of all cells
% 
% fucci_flashing = events_info(:,2:6);
% nonflashingcells_flip = flip(nonflashingcells);
% for n=nonflashingcells_flip
%     fucci_flashing(n,:) = [];
% end
% fucci_flashing = fucci_flashing(:,5);
% fucci_flashing_log = log(fucci_flashing);
% 
% 
% fig4 = figure('Position', [100 20 2540 1300]);
% H.Nodes.NodeColors = fucci_flashing; %can set variable here to use to colorcode nodes along a gradient
% p = plot(H,'NodeCData',H.Nodes.NodeColors, 'MarkerSize',7);
% colorbar
% cutName=strfind(fileName,'.tif');
% outputFile1=[fileName(1:cutName-1) '_graphfucci1.tif'];
% fig.PaperUnits = 'normalized';
% fig.PaperPosition = [0 0 27 13];
% print(fig4,'-dtiff',outputFile1,'-r300')


%% code if I want to highlight G1 cells versus not
% fucci_flashing_red = fucci_flashing;
% listoffucci2 = [];
% fuccithresh = quantile(events_info(:,6),0.2);
% for n=1:numel(fucci_flashing_red)
%     if fucci_flashing_red(n,:) > fuccithresh
%         listoffucci2(n,1) = 1;
%     else
%         listoffucci2(n,1) = 0;
%     end
% end
% 
% for n=1:cell_number
%     if events_info(n,6) > fuccithresh
%         events_info(n,8) = 1;
%     else
%         events_info(n,8) = 0;
%     end
% end
% 
% 
% H.Nodes.Name = string(flashingcells.');
% H.Nodes.NodeColors = listoffucci2;
% 
% fig3 = figure('Position', [100 20 2540 1300]);
% p = plot(H, 'MarkerSize',7,'NodeCData', H.Nodes.NodeColors);
% p.NodeCData = H.Nodes.NodeColors;
% colorbar
% cutName=strfind(fileName,'.tif');
% outputFile1=[fileName(1:cutName-1) '_graphfucci2.tif'];
% fig.PaperUnits = 'normalized';
% fig.PaperPosition = [0 0 27 13];
% print(fig3,'-dtiff',outputFile1,'-r300')

%% save variables

% variables to save if not running fucci analysis
outputFile2=[fileName(1:cutName-1) '_ALL_CELLS_plus.mat'];
save(outputFile2,'threshold', 'A_in', 'A_or', 'adjacencymatrix', 'areaList', 'areas', 'avg', 'avg_area', 'bins', 'bins3', 'bins4', 'bkg', 'CC', 'CC2', 'cell_number', 'cell_per', 'cells', 'cells_mean', 'corrThresh', 'cutName', 'd1', 'd2', 'distances', 'eventDownstream', 'events_info', 'deltaf_ft', 'deltaf_f', 'f_f0t_flashing', 'f_f0t_nonflashing', 'fileName', 'fitfun', 'flashingcells', 'G', 'nbList', 'nonflashingcells', 'npil_mean', 'numbers', 'overallEvents', 'p', 'pid', 'pixel_total', 'pixelLengthX', 'pixelLengthY', 'prog', 'threshT', 'threshTe', 'tid', 'updatedEvents', 'xdata', 'Z');

% outputFile2=[fileName(1:cutName-1) '_ALL_CELLS_plus.mat'];
% save(outputFile2,'filename2','fileName2','fucci_flashing','fucci_flashing_red','fuccithresh','threshold', 'A_in', 'A_or', 'adjacencymatrix', 'areaList', 'areas', 'avg', 'avg_area', 'bins', 'bins3', 'bins4', 'bkg', 'CC', 'CC2', 'cell_number', 'cell_per', 'cells', 'cells_mean', 'corrThresh', 'cutName', 'd1', 'd2', 'distances', 'eventDownstream', 'events_info', 'deltaf_ft', 'deltaf_f', 'f_f0t_flashing', 'f_f0t_nonflashing', 'fileName', 'fitfun', 'flashingcells', 'G', 'nbList', 'nonflashingcells', 'npil_mean', 'numbers', 'overallEvents', 'p', 'pid', 'pixel_total', 'pixelLengthX', 'pixelLengthY', 'prog', 'threshT', 'threshTe', 'tid', 'updatedEvents', 'xdata', 'Z');

load gong
sound(y,Fs) 