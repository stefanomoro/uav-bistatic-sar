function [radar] = getCrossTalk(radar)
%GETCROSSTALK Find Crosstalk in the RC as the maximum value for each column 
%   [radar] = getCrossTalk(radar)
 
[~,cross_talk_idxs] = max(radar.RC,[],1);
% apply median filter and then moving average
cross_talk_idxs = round(movmean(medfilt1(cross_talk_idxs,1e3) ,1e3)); 

cross_talk = zeros(size(radar.RC,2),1);
for i = 1:length(cross_talk)
    cross_talk(i) = radar.RC(cross_talk_idxs(i),i);
end
radar.cross_talk.value = cross_talk;
end
