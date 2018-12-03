%% Demo compute correspondences

% Add MeshMonk's toolbox to the working path and setup current folder
addpath(genpath('path\to\meshmonk')) % Set to location of meshmonk

studypath = 'path\to\DemoFolder\';   % Set to location of demo material
cd(studypath);

%% Load template shape (floating) and target'

floatingPath = [studypath '/Template.obj'];
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFeatures = [floatingPoints, 1/sqrt(3.0)*ones(size(floatingPoints))];
floatingFeatures = single(floatingFeatures);
floatingFaces = uint32(floatingFaces-1);    %-1 to make it compatible with C++ indexing 
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

targetPath = [studypath '/demoFace.obj'];
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetFeatures = single([targetPoints, 1/sqrt(3.0)*ones(size(targetPoints))]);
targetFaces = uint32(targetFaces-1);        %-1 to make it compatible with C++ indexing 
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Test the computation of correspondences

%# Set Parameters
correspondencesSymmetric = true;
correspondencesNumNeighbours = 3;
correspondencesFlagThreshold = 0.9;
correspondencesEqualizePushPull = false;
correspondingFeatures = floatingFeatures;
correspondingFlags = floatingFlags;

compute_correspondences(floatingFeatures, targetFeatures,...
                        floatingFlags, targetFlags,...
                        correspondingFeatures, correspondingFlags,...
                        correspondencesSymmetric, correspondencesNumNeighbours,...
                        correspondencesFlagThreshold, correspondencesEqualizePushPull);
                            
%% Visualize
s = shape3D;
s.Vertices = correspondingFeatures(:,1:3);
st = shape3D;
st.Vertices = targetFeatures(:,1:3);
v = viewer(s);  
viewer(st,v);                          
