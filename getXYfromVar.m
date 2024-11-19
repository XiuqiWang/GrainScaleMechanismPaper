function [x,y]=getXYfromVar(Var, binWidth)
binEdges = min(Var):binWidth:max(Var);
[frequencies, edges] = histcounts(Var,binEdges,'Normalization', 'probability');
binCenters = (edges(1:end-1) + edges(2:end)) / 2;  % Calculate bin centers
y = [0, frequencies, 0];  % Add 0 at the beginning and end for y-axis
x = [edges(1), binCenters, edges(end)];  % Extend the x-axis
end