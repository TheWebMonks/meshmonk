clear all

%% Load Data

%Load a mesh
floatingPath = '/home/jonatan/projects/meshmonk/examples/faceTemplate.obj';
floatingPath = '/home/jonatan/projects/meshmonk/examples/ExPeter/floating.obj';
%floatingPath = '/home/jonatan/projects/meshmonk/examples/data/bunny.obj';
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing
floatingPoints = single(floatingPoints);
floatingNormals = single(zeros(size(floatingPoints)));
compute_normals(floatingPoints, floatingFaces, floatingNormals);
floatingFeatures = single([floatingPoints, floatingNormals]);
numFloatingElements = size(floatingFeatures,1);
numFloatingFaces = size(floatingFaces,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

%% Parameters and help variables
downsampleRatio = 0.7;
numDownsampledElements = uint32(round((1.0-downsampleRatio) * numFloatingElements));
numDownsampledFaces = numFloatingFaces; %can't really predict this number so we'll just use the upper limit estimation!

floatingFeaturesSampled = single(zeros(numDownsampledElements,6));
floatingFacesSampled = uint32(zeros(numDownsampledFaces,3));
floatingFlagsSampled = single(zeros(numDownsampledElements,1));
originalIndices = uint32(zeros(numDownsampledElements,1));

%% Downsample
downsample_mesh(floatingFeatures, floatingFaces, floatingFlags,...
                floatingFeaturesSampled, floatingFacesSampled, floatingFlagsSampled,...
                originalIndices, downsampleRatio);
            
%# Clean the output
[floatingFeaturesSampled,numDownsampledElements] = clean_downsampled_features(floatingFeaturesSampled);
floatingFacesSampled = clean_downsampled_faces(floatingFacesSampled);
floatingFlagsSampled = floatingFlagsSampled(1:numDownsampledElements);
originalIndices = originalIndices(1:numDownsampledElements);