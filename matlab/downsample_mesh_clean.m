function [ floatingFeaturesSampled,floatingFacesSampled, floatingFlagsSampled, originalIndices ]...
    = downsample_mesh_clean( floatingFeatures, floatingFaces, floatingFlags, downsampleRatio )
%DOWNSAMPLE_MESH_CLEAN Downsampling plus cleaning of output
%   The downsampling comes from an external library (meshmonk), which
%   requires static memory allocation. So a surplus of memory should be
%   allocated, and in this function is cleaned up after the initial
%   downsampling.

% Data and parameters
numDownsampledElements = uint32(round((1.0-downsampleRatio) * numFloatingElements));
numDownsampledFaces = numFloatingFaces; %can't really predict this number so we'll just use the upper limit estimation!

floatingFeaturesSampled = single(zeros(numDownsampledElements,6));
floatingFacesSampled = uint32(zeros(numDownsampledFaces,3));
floatingFlagsSampled = single(zeros(numDownsampledElements,1));
originalIndices = uint32(zeros(numDownsampledElements,1));

%# Downsample
downsample_mesh(floatingFeatures, floatingFaces, floatingFlags,...
                floatingFeaturesSampled, floatingFacesSampled, floatingFlagsSampled,...
                originalIndices, downsampleRatio);
            
%# Clean the output
[floatingFeaturesSampled,numDownsampledElements] = clean_downsampled_features(floatingFeaturesSampled);
floatingFacesSampled = clean_downsampled_faces(floatingFacesSampled);
floatingFlagsSampled = floatingFlagsSampled(1:numDownsampledElements);
originalIndices = originalIndices(1:numDownsampledElements);

end

