%% Post-step 1

% in Workspace, A_or contains masks of all components in the form of sparse matrices,
% C_or contains dF/F% traces of all components
C_df_full = full(C_df);

N = size(C_df_full, 1); % # of components
xdata = 1:size(C_df_full, 2); % timepoints
overallEvents = [];
fprintf(1, '\nProgress:     ')
for i=1:N
    [pid, tid] = Find1TracePeaksTroughs(C_df_full(i,:), 0.2);
    if isempty(pid)
        continue;
    end
    [events, fitfun] = MergePeaks2FitGaussians(C_df_full(i,:), pid, tid);
    events(:,1) = i; % fill in the component info of these events
    overallEvents = [overallEvents; events];
     figure(100), hold off
     for j=1:size(fitfun, 1)
         plot(fitfun{j}, xdata', C_df_full(i,:));
         hold on
     end
     
     
     plot(pid, C_df_full(i, pid), 'x', 'MarkerSize', 8)
     tid_flat = unique(tid(:));
     plot(tid_flat, C_df_full(i, tid_flat), 'o', 'MarkerSize', 8)
     ax=gca; ax.Title.String = sprintf("Trace %d", i);
%      pause
    prog = 100 * i/N;
    fprintf(1,'\b\b\b\b%3.0f%%',prog)
end
fprintf(1, '\n')

%% Post-step 2

[nbList, areaList] = findNeighbors(A_or, d2, d1, 2);
% 2 pixels is the distance threshold to consider being neighbor; used to
% be 25 but this picked up components far away

updatedEvents = overallEvents; %overallEvents is from step1 and in the workspace
corrThresh = 0.75;  % Threshold of correlation coefficient to be considered correlated enough
N_comp = numel(C_df_full(:,1));
% N_comp is the number of spatial components in the workspace
% assign area to each event based on component/ROI it belongs to
for compI=1:N_comp
    updatedEvents(updatedEvents(:,1) == compI, 5) = areaList(compI);
end


%% Post-step 3

events_info = [];
N = size(C_df_full, 1); % # of components
for i=1:N
    events_info(i,1) = i; % Component ID
    events_info(i,2) = sum(updatedEvents(:,1) == i); % Number of peaks
    events_info(i,3) = areaList(1,i); % Component area
    events_info(i,4) = events_info(i,2)/(T*2); % Events per second
    events_info(i,5) = events_info(i,2)/(T/60*2); % Events per minute
end
N = size(C_df_full, 1); % # of components
for i=1:N
    if ismember(i,updatedEvents(:,1)) == 1 % need to skip cells that were thrown out for messy calcium traces
        events_info(i,6) = max(updatedEvents((updatedEvents(:,1)==i),3)); % max lingering time
    else
        continue
    end
end
avg_area = mean(areaList);
N = size(updatedEvents, 1); % update N


%% Post-step 4

% Threshold for start time difference. If event J's start time is after event I but 
% within threshT, then it's considered event I's downstream event.
threshT = 5; 
% If event J's start time is <threshTe after event I's end time, then it's
% also considered event I's downstream event.
threshTe = 5;

eventDownstream = cell(size(updatedEvents, 1), 1); 

nR = numel(C_df_full(:,1)); % nR is number of components
adjacencymatrix = zeros(nR);


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
fig1 = figure('Position', [0 0 2560 1440]);
plot(G)
cutName=strfind(fileName,'.tif');
outputFile1=[fileName(1:cutName-1) '_graphcaiman.tif'];
fig.PaperUnits = 'normalized';
fig.PaperPosition = [0 0 27 13];
print(fig1,'-dtiff',outputFile1,'-r300')

fig1 = figure('Position', [0 0 2560 1440]);
p = plot(G);
numbers = 1:nR;
p.NodeLabel= numbers;
outputFile1=[fileName(1:cutName-1) '_graphnumberedcaiman.tif'];
fig.PaperUnits = 'normalized';
fig.PaperPosition = [0 0 27 13];
print(fig1,'-dtiff',outputFile1,'-r300')

%% calculating neighborhood size based on graph

bins = conncomp(G);
bins2 = bins.'; % bins connected nodes

updatedEvents(:,6)=zeros; % assigns bin number to each event
for var=1:nR
    updatedEvents((updatedEvents(:,1)==var),6) = bins2(var,1);
end


nB = max(bins); % total number of bins
bins3 = zeros(nB,1);
updatedEvents(:,7) = zeros;
bins3(:,2) = zeros;
updatedEvents(:,8) = ones;
updatedEvents(:,9) = zeros;
for var2=1:nB
    bins3(var2,1) = sum(updatedEvents((updatedEvents(:,6)==var2),5));
    updatedEvents((updatedEvents(:,6)==var2),7) = bins3(var2,1); %combines areas of connected nodes
    bins3(var2,2) = sum(bins(1,:)==var2);
    updatedEvents((updatedEvents(:,6)==var2),9) = bins3(var2,2); %counts number of nodes
end

updatedEvents(:,8) = zeros;
for i=1:N
    updatedEvents(i,8) = updatedEvents(i,7)/avg_area; % Neighborhood size in # of cells based on average component area
end

bins4 = zeros(nR,1); % going to create new variable that assigns size of neighborhood each cell is in
for var2=1:nR
    bins4(var2,1) = bins3((bins2(var2,1)),2);
    updatedEvents((updatedEvents(:,1)==var2),7) = bins4(var2,1); %counts number of nodes
end
events_info(:,7) = bins4(:,1);

updatedEvents(:,10) = zeros;
updatedEvents(:,11) = zeros;
for i=1:nR
    updatedEvents((updatedEvents(:,1)==i),10) = max(updatedEvents((updatedEvents(:,1)==i),2))-min(updatedEvents((updatedEvents(:,1)==i),2)); % timesteps between first and last calcium events' starts
    updatedEvents((updatedEvents(:,1)==i),11) = events_info(i,2)/((max(updatedEvents((updatedEvents(:,1)==i),10))/60*2)); % peaks per minute for timeperiod of activity
end

for i=1:nR
     if ismember(i,updatedEvents(:,1)) == 1 % need to skip cells that were thrown out for messy calcium traces
        events_info(i,7) = max(updatedEvents((updatedEvents(:,1)==i),9)); % neighborhood size
    else
        continue
    end
end

%% save variables

outputFile2=[fileName(1:cutName-1) '_caiman.mat'];
save(outputFile2,'A','A2','A_keep','A_or','adjacencymatrix','Ain','Am','ans','areaList','avg_area','b','b2','bins','bins2','bins3','cutName','C','C2','C_df','C_df_full','C_keep','C_or','center','center_keep','Cin','Cm','Cn','compI','compJ','Coor','corrThresh','d','d1','d2','display_merging','eventDownstream','events_info','f','f2','file','fileName','fin','fitness','G','i','idx','ind_exc','json_file','K','K_m','keep','merged_ROIs','N','N_comp','nB','nbList','nR','numbers','options','overallEvents','p','P','P2','P_or','path','Pm','S','S2','S_or','Sm','T','tau','updatedEvents','xdata');

load gong
sound(y,Fs)