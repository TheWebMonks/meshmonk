function [ cleanedFaces ] = clean_downsampled_faces( uncleanedFaces )
%clean_downsample_faces should be used on the output faces of downsample_mesh
%function
%   Because downsample_mesh is a mexed c-function, memory allocation is not
%   dynamic. So we provide arrays as input that are too large, and use this
%   function to trim the outputs.

%# Determine finalValidIndex: the index after which we should cut off
%# elements. We do that by looping over the faces and stopping when we find
%# three consecutive zeros.
startIndex = length(nonzeros(uncleanedFaces));
finalValidIndex = numel(uncleanedFaces);
if (startIndex < numel(uncleanedFaces)-2)
    for i=startIndex:numel(uncleanedFaces)-2
        if (uncleanedFaces(i) == 0) && (uncleanedFaces(i+1) == 0) && (uncleanedFaces(i+2) == 0)
            finalValidIndex = i - 1;
            break;
        end
    end
end

%# If by coincidence the final valid element was a zero, we would end up
%# cutting it off. In this case, the number of elements could not be
%# divided by three. So let's check the modulus and fix this (rare) case.
remainder = mod(finalValidIndex,3);
%## The value of remainder gives the amount of times that index zero
if (remainder > 0)
    finalValidIndex = finalValidIndex + remainder;
end

%# Cut off the faces matrix
cleanedFaces = uint32(zeros(finalValidIndex/3,3));
cleanedFaces(1:finalValidIndex/3) = uncleanedFaces(1:finalValidIndex/3);
cleanedFaces(finalValidIndex/3+1:2*finalValidIndex/3) = uncleanedFaces(finalValidIndex/3+1:2*finalValidIndex/3);
cleanedFaces(2*finalValidIndex/3+1:3*finalValidIndex/3) = uncleanedFaces(2*finalValidIndex/3+1:3*finalValidIndex/3);

end

