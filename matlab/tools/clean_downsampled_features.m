function [ cleanedFeatures, numElements ] = clean_downsampled_features( uncleanedFeatures )
%clean_downsampled_features should be used on the output features of downsample_mesh
%function
%   Because downsample_mesh is a mexed c-function, memory allocation is not
%   dynamic. So we provide arrays as input that are too large, and use this
%   function to trim the outputs.

numCols = 6;

%# Determine finalValidIndex: the index after which we should cut off
%# elements. We do that by finding the last (highest) index that is nonzero

finalValidIndex =find(uncleanedFeatures,1, 'last');

%startIndex = length(nonzeros(uncleanedFeatures));
%finalValidIndex = numel(uncleanedFeatures);
% if (startIndex < numel(uncleanedFeatures)-2)
%     for i=startIndex:numel(uncleanedFeatures)-2
%         if (uncleanedFeatures(i) == 0.0) && (uncleanedFeatures(i+1) == 0.0) && (uncleanedFeatures(i+2) == 0.0)
% 
%             finalValidIndex = i - 1;
%             break;
%         end
%     end
% end

%# If by coincidence the final valid element was a zero, we would end up
%# cutting it off. In this case, the number of elements could not be
%# divided by six (number of columns). So let's check the modulus and fix this (rare) case.
remainder = mod(finalValidIndex,numCols);
%## The value of remainder gives the amount of times that index zero
if (remainder > 0)
    finalValidIndex = finalValidIndex + remainder;
end

%# Cut off the faces matrix
cleanedFeatures = single(zeros(finalValidIndex/numCols,numCols));
for i=0:(numCols-1)
    cleanedFeatures(i*finalValidIndex/numCols + 1:(i+1)*finalValidIndex/numCols) = ...
        uncleanedFeatures(i*finalValidIndex/numCols + 1:(i+1)*finalValidIndex/numCols);
end

numElements = finalValidIndex/numCols;

end

