clear all

%Load a mesh
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file('/home/jonatan/projects/meshmonk/examples/faceTemplate.obj');
floatingFeatures = [floatingPoints, 1/3.0*ones(size(floatingPoints))];
floatingFeatures = single(floatingFeatures);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing ?
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

%Load a mesh
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file('/home/jonatan/projects/meshmonk/examples/faceTarget.obj');
targetFeatures = single([targetPoints, 1/3.0*ones(size(targetPoints))]);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing ?
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Try the mexed pyramid_registration
%mex pyramid_registration.cpp -lmeshmonk

%# Set Parameters
numIterations = 20;
correspondencesSymmetric = true;
correspondencesNumNeighbours = 5;
inlierKappa = 3.0;

floatingFeatures = rigid_registration(floatingFeatures, targetFeatures,...
                                floatingFaces, targetFaces,...
                                floatingFlags, targetFlags,...
                                numIterations,...
                                correspondencesSymmetric, correspondencesNumNeighbours,...
                                inlierKappa);
                            
%% Write Result
vertface2obj(floatingFeatures(:,1:3),floatingFaces,'/home/jonatan/projects/meshmonk/examples/matlabResult.obj')
                            
                          