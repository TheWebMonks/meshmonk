%% Demo rigid registration modular

% Add path to compiled MEX files (assuming mex_all.m was run in matlab/ directory)
addpath(fullfile('..', 'matlab'));

%% Load face template (floating) and face (target)

floatingPath = 'Template.obj';
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing
floatingPoints = single(floatingPoints);
floatingNormals = single(zeros(size(floatingPoints)));
% compute_normals(floatingPoints, floatingFaces, floatingNormals);
floatingFeatures = single([floatingPoints, floatingNormals]);
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
%DEBUG
floatingFlags(1:1000) = 0.0;
%END DEBUG
clear floatingPoints;

targetPath = 'demoFace.obj';
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetPoints = single(targetPoints);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing
targetNormals = single(zeros(size(targetPoints)));
% compute_normals(targetPoints, targetFaces, targetNormals);
targetFeatures = single([targetPoints, targetNormals]);
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%Set up resulting transformation matrix
transformationMatrix = single(eye(4,4));
tempTransformationMatrix = single(eye(4,4));

%% Try modular rigid registration
%# Set Parameters
numIterations = 80;
correspondencesSymmetric = true;
correspondencesNumNeighbours = 3;
correspondencesFlagThreshold = 0.9;
correspondencesEqualizePushPull = false;
inlierKappa = 12.0;
inlierUseOrientation = true;
useScaling = false;

%# Initialize data structures
correspondingFeatures = single(zeros(numFloatingElements,6));
correspondingFlags = single(ones(numFloatingElements,1));
inlierWeights = single(ones(numFloatingElements,1));

%% Registration

%# Iterative Registration process
for i=1:numIterations
    %# Compute Correspondences
    %[correspondingFeatures, correspondingFlags] = compute_correspondences(floatingFeatures, targetFeatures,...
    compute_correspondences(floatingFeatures, targetFeatures,...
                            floatingFlags, targetFlags,...
                            correspondingFeatures, correspondingFlags,...
                            correspondencesSymmetric, correspondencesNumNeighbours,...
                            correspondencesFlagThreshold, correspondencesEqualizePushPull);
    
    %# Compute Inlier Weights
    compute_inlier_weights(floatingFeatures, correspondingFeatures,...
                           correspondingFlags, inlierWeights,...
                           inlierKappa, inlierUseOrientation);
    
    %# Compute Transformation
    compute_rigid_transformation(floatingFeatures, correspondingFeatures,...
                                    inlierWeights, tempTransformationMatrix,...
                                    useScaling);
                                
    %# Update transformation matrix
    transformationMatrix = tempTransformationMatrix * transformationMatrix;
end


%% Visualize
s = shape3D;
s.Vertices = floatingFeatures(:,1:3);
st = shape3D;
st.Vertices = targetFeatures(:,1:3);
v = viewer(s);  
viewer(st,v);               
