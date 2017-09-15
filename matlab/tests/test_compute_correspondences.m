clear all

%Load a mesh
floatingPath = '/home/jonatan/projects/MIRC/meshmonk/matlab_peter/MyTestMatlab/DATA/HANNE/mappingTemplate.obj';
%floatingPath = '/home/jonatan/projects/meshmonk/examples/ExPeter/floating.obj';
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFeatures = [floatingPoints, 1/sqrt(3.0)*ones(size(floatingPoints))];
floatingFeatures = single(floatingFeatures);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing ?
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

%Load a mesh
targetPath = '/home/jonatan/projects/MIRC/meshmonk/matlab_peter/MyTestMatlab/DATA/HANNE/F069MG04.obj';
%targetPath = '/home/jonatan/projects/meshmonk/examples/ExPeter/target.obj';
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetFeatures = single([targetPoints, 1/sqrt(3.0)*ones(size(targetPoints))]);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing ?
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Test the computation of correspondences

%# Set Parameters
correspondencesSymmetric = true;
correspondencesNumNeighbours = 5;
correspondencesFlagThreshold = 0.999;
correspondencesEqualizePushPull = false;
correspondingFeatures = floatingFeatures;
correspondingFlags = floatingFlags;

compute_correspondences(floatingFeatures, targetFeatures,...
                        floatingFlags, targetFlags,...
                        correspondingFeatures, correspondingFlags,...
                        correspondencesSymmetric, correspondencesNumNeighbours,...
                        correspondencesFlagThreshold, correspondencesEqualizePushPull);

                            
%% Write Result
vertface2obj(correspondingFeatures(:,1:3),floatingFaces,'/home/jonatan/projects/meshmonk/examples/matlabResult.obj')
                            
                          