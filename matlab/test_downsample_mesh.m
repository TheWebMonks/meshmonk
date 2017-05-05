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
clear floatingPoints floatingNormals;

%% Parameters and help variables
downsampleRatio = 0.7;

%% Downsample
[ floatingFeaturesSampled,floatingFacesSampled, floatingFlagsSampled, originalIndices ]...
    = downsample_mesh_clean( floatingFeatures, floatingFaces, floatingFlags, downsampleRatio );