function [peakIds, srTroughs, eventsOfThisTrace, varargout] = Find1TracePeaksTroughs(trace, thresh)
%Find1TracePeaksTroughs Locate peaks and troughs flanking the peaks
%   Detailed explanation goes here

traceMean = median(trace);
trace = trace - traceMean;

% Find peaks by dilation and thresholding
se = strel('disk', 2);
traceD=imdilate(trace, se);
% filter peaks by percentage of the maximum of current trace:
peakIds = find(((traceD<=trace) .* trace) > max(trace) * thresh); % thresh was default to 20% 

% Remove border peakIds
peakIds = peakIds( ~(peakIds == 1 | peakIds == length(trace)) );

% Find troughs by erosion and thresholding
se = strel('disk', 2);
traceE=imerode(trace, se);
troughIds = find(traceE>=trace);

if nargout > 0
    varargout{1} = troughIds;
end

srTroughs = findSurroundingTroughs(trace, peakIds, troughIds);

% Detect the troughs that's very high and not followed by a peak. For those
% troughs, extend it towards the next deeper trough
srTroughs = detectWrongTroughs(trace, peakIds, srTroughs);

end
