function [ featuresSampled, facesSampled, flagsSampled, originalIndices ]...
    = downsample_mesh_clean( features, faces, flags, downsampleRatio )
%DOWNSAMPLE_MESH_CLEAN Downsampling plus cleaning of output
%   The downsampling comes from an external library (meshmonk), which
%   requires static memory allocation. So a surplus of memory should be
%   allocated, and in this function is cleaned up after the initial
%   downsampling.

% Data and parameters
numElements = size(features,1);
numFaces = size(faces,1);
numDownsampledElements = uint32(round((1.0-downsampleRatio) * numElements));
numDownsampledFaces = numFaces; %can't really predict this number so we'll just use the upper limit estimation!

featuresSampled = single(zeros(numDownsampledElements,6));
facesSampled = uint32(zeros(numDownsampledFaces,3));
flagsSampled = single(zeros(numDownsampledElements,1));
originalIndices = uint32(zeros(numDownsampledElements,1));

%# Downsample
downsample_mesh(features, faces, flags,...
                featuresSampled, facesSampled, flagsSampled,...
                originalIndices, downsampleRatio);
            
%# Clean the output
[featuresSampled,numDownsampledElements] = clean_downsampled_features(featuresSampled);
facesSampled = clean_downsampled_faces(facesSampled);
flagsSampled = flagsSampled(1:numDownsampledElements);
originalIndices = originalIndices(1:numDownsampledElements);

end

