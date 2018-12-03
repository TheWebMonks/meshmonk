%% Demo downsample mesh

% Add MeshMonk's toolbox to the working path and setup current folder
addpath(genpath('path\to\meshmonk')) % Set to location of meshmonk

studypath = 'path\to\DemoFolder\';   % Set to location of demo material
cd(studypath);

%% Load Data

%Load a mesh
floatingPath = [studypath '/demoFace.obj'];
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
 
%% Visualize
s = shape3D;
s.Vertices = floatingFeatures(:,1:3);
st = shape3D;
st.Vertices = floatingFeaturesSampled(:,1:3);
v = viewer(s);  
viewer(st,v);                          