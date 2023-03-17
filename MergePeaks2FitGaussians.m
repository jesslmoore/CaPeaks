function [events, varargout] = MergePeaks2FitGaussians(trace, peakIds, srTroughs)
%MergePeaks2FitGaussians Merge the peaks based on magnitude of the troughs between
% them; then fit one Gaussian curve for each merged peak using only data points
% marked by the leftmost and rightmost troughs
%
% 
nT = length(trace);
nPeaks = length(peakIds);
thresh = 15;
intThresh = 0.5;
evConn = false(1, nPeaks-1);
for i = 1:nPeaks-1
    trough1 = trace(srTroughs(i, 2));
    % If peak_i and peak_i+1 are close enough and trough between them high
    % enough, then this two peaks should be considered connected
    if peakIds(i+1) - peakIds(i) < thresh && ...
            (trough1 > trace(peakIds(i)) * intThresh || ...
            trough1 > trace(peakIds(i+1)) * intThresh)
        evConn(i) = true;
    end
end

indivEvId = 1;
mergedPeaks = cell(500,1); % assume no more than 50 flashes at each trace
nMergedPeaks = 1;
while indivEvId < nPeaks
    mergedPeaks{nMergedPeaks} = [mergedPeaks{nMergedPeaks} indivEvId];
    while indivEvId < nPeaks && evConn(indivEvId)
        % if this peakId is supposed to be connected to peakId+1, as
        % indicated by peakConn, then include peakId+1 to the current
        % event. And keep doing till peakConn indicates a break point.
        mergedPeaks{nMergedPeaks} = [mergedPeaks{nMergedPeaks} indivEvId+1];
        indivEvId = indivEvId + 1;
    end
    % The run of connected peaks is over; advance peakId and event Id
    % to get ready for the next run
    indivEvId = indivEvId + 1;
    nMergedPeaks = nMergedPeaks + 1;
    if indivEvId == nPeaks
        mergedPeaks{nMergedPeaks} = indivEvId;
    elseif indivEvId > nPeaks
        % in this case, number of events is over counted; decrease it
        % by 1:
        nMergedPeaks = nMergedPeaks - 1;
    end
end

% There is one special case, when nPeaks==1, that the above while loop
% cannot handle properly:
if nPeaks == 1
    mergedPeaks{1} = indivEvId;
end

%% Then fit a single Gaussian for all groups of merged peaks to help defined linger time.
indivEvents = zeros(nMergedPeaks, 3);
fitfunc = cell(nMergedPeaks,1);
to_discard = [];
for i=1:nMergedPeaks
    peaksInThisGroup = mergedPeaks{i};
    t1 = srTroughs(peaksInThisGroup(1), 1);
    t2 = srTroughs(peaksInThisGroup(end), 2);
    try
        [fitfunc{i}, gof] = FitGaussians(trace, t1, t2);
        if gof.rsquare > 0.5
            sigma = fitfunc{i}.c;
            center = fitfunc{i}.b;
            amp = fitfunc{i}.a;
        else
            % not trustworthy fit
            sigma = (t2-t1)/2;
            center = (t1+t2)/2;
            amp = fitfunc{i}.a;
        end
        indivEvents(i, :) = [center sigma amp];
    catch ME
        if strcmp(ME.identifier, 'curvefit:fit:InsufficientData')
            warning("Peak discarded for having too few data points for the fit")
            to_discard = [to_discard i];
        % Lin added next four lines 2020-08-31:
        elseif strcmp(ME.identifier, 'curvefit:checkbounds:lowerBoundExceedsUpperBound')...
            || strcmp(ME.identifier, 'curvefit:checkbounds:lowerBoundsExceedsUpperBounds')
            warning("Peak discarded for lower bounds exceeding uppper bounds")
            to_discard = [to_discard i];
        end
    end
end
indivEvents(to_discard, :) = [];
fitfunc(to_discard) = [];
nMergedPeaks = size(indivEvents, 1);

%% Finally, merge these fit Gaussians if their envelops overlap
% Using the same merging algorithm as above

% Lin added next 7 lines 2020-08-31:
if nMergedPeaks < 1
    events = [];
    if nargout > 0
        varargout{1} = {};
    end
    return
end

evConn = false(1, nMergedPeaks - 1);
for i = 1:numel(evConn)
    rightEnd = indivEvents(i, 1) + indivEvents(i, 2);
    leftEnd = indivEvents(i+1, 1) - indivEvents(i+1, 2);
    if leftEnd < rightEnd  || leftEnd - rightEnd < 5 % play with this threshold
        evConn(i) = true;
    end
end

indivEvId = 1;
eventsIDs = cell(50, 1); % assume no more than 50 final events at this ROI
nEvents = 1;  % number of events
while indivEvId < nMergedPeaks
    eventsIDs{nEvents} = [eventsIDs{nEvents} indivEvId];
    while indivEvId < nMergedPeaks && evConn(indivEvId)
        % if this peakId is supposed to be connected to peakId+1, as
        % indicated by peakConn, then include peakId+1 to the current
        % event. And keep doing to peakConn indicates a break point.
        eventsIDs{nEvents} = [eventsIDs{nEvents} indivEvId+1];
        indivEvId = indivEvId + 1;
    end
    % The run of connected peaks is over; advance peakId and event Id
    % to get ready for the next run
    indivEvId = indivEvId + 1;
    nEvents = nEvents + 1;
    if indivEvId == nMergedPeaks
        eventsIDs{nEvents} = indivEvId;
    elseif indivEvId > nMergedPeaks
        % in this case, number of events is over counted; decrease it
        % by 1:
        nEvents = nEvents - 1;
    end
end

if nMergedPeaks == 1
    eventsIDs{1} = indivEvId;
end

events = zeros(nEvents, 5);
for i=1:nEvents
    evId = eventsIDs{i}(1);
    t0 = indivEvents(evId, 1) - indivEvents(evId, 2);
    evId = eventsIDs{i}(end);
    t1 = indivEvents(evId, 1) + indivEvents(evId, 2);
    t0 = max(t0, 1); t1 = min(nT, t1);
    intensityMax = max(trace(uint16(t0):uint16(t1)));
    events(i, 2:4) = [t0 t1-t0 intensityMax]; % keep getting error here: Unable to perform assignment because the size of the left side is 1-by-3 and the size of the right side is 1-by-2.
%     [eventsIDs{i}   events(i, 2:3)]
end

if nargout > 0
    varargout{1} = fitfunc;
end

end

function [fitFunc, gof] = FitGaussians(trace, t1, t2)
%Fit a Gaussian curve using trace(t1:t2)

ft = 'a*exp(-((x-b)/c)^2)+d';

xdata = ( t1:t2 )';
startGuessA = max(trace(t1:t2));
startGuessB = (t2 + t1) * 0.5;
startGuessC = (t2-t1) * 0.5;
startGuessD = 0;
lowerBounds = [startGuessA *.8, t1, 1, 0];
upperBounds = [startGuessA * 1.2, t2, startGuessC*2, mean(trace)];
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower',lowerBounds, 'Upper', upperBounds,...
    'StartPoint',[ startGuessA startGuessB startGuessC startGuessD]);
[fitFunc,gof] = fit(xdata, trace(t1:t2)', ft, fo);

end
