%% Demo nonrigid registration modular
 
% Add path to compiled MEX files (assuming mex_all.m was run in matlab/ directory)
addpath(fullfile('..', 'matlab'));

%% Load face template (floating) and face (target)

floatingPath = 'Template.obj';
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing
floatingPoints = single(floatingPoints);
floatingNormals = single(zeros(size(floatingPoints)));
floatingFeatures = single([floatingPoints, floatingNormals]);
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
%DEBUG
floatingFlags(1:1000) = 0.0;
%END DEBUG
clear floatingPoints;

%Load a mesh
targetPath = 'demoFace.obj';
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetPoints = single(targetPoints);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing
targetNormals = single(zeros(size(targetPoints)));
targetFeatures = single([targetPoints, targetNormals]); 
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Try the mexed pyramid_registration
%mex pyramid_registration.cpp -lmeshmonk

%# Set Parameters
numIterations = 500;
correspondencesSymmetric = true;
correspondencesNumNeighbours = 3;
correspondencesFlagThreshold = 0.9;
correspondencesEqualizePushPull = false;
inlierKappa = 12.0;
inlierUseOrientation = true;
transformSigma = 3.0;
transformNumViscousIterationsStart = numIterations;
transformNumViscousIterationsEnd = 1;
transformNumElasticIterationsStart = numIterations;
transformNumElasticIterationsEnd = 1;

%# Derived parameters
%## Annealing
viscousAnnealingRate = exp(log(single(transformNumViscousIterationsEnd)/single(transformNumViscousIterationsStart))/(numIterations-1));
elasticAnnealingRate = exp(log(single(transformNumElasticIterationsEnd)/single(transformNumElasticIterationsStart))/(numIterations-1));
numViscousIterations = transformNumViscousIterationsStart;
numElasticIterations = transformNumElasticIterationsStart;

%# Initialize data structures
correspondingFeatures = single(zeros(numFloatingElements,6));
correspondingFlags = single(ones(numFloatingElements,1));
inlierWeights = single(ones(numFloatingElements,1));

%%

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
    compute_nonrigid_transformation(floatingFeatures, correspondingFeatures,...
                                    floatingFaces,...
                                    floatingFlags, inlierWeights,...
                                    10, transformSigma,...
                                    numViscousIterations, numElasticIterations);
                                
    %# Annealing
    numViscousIterations = uint32(round(transformNumViscousIterationsStart * viscousAnnealingRate^(i)));
    if (numViscousIterations < transformNumViscousIterationsEnd) numViscousIterations = transformNumViscousIterationsEnd;end
    numElasticIterations = uint32(round(transformNumElasticIterationsStart * elasticAnnealingRate^(i)));
    if (numElasticIterations < transformNumElasticIterationsEnd) numElasticIterations = transformNumElasticIterationsEnd;end
end

%% Visualize
s = shape3D;
s.Vertices = correspondingFeatures(:,1:3);
st = shape3D;
st.Vertices = targetFeatures(:,1:3);
v = viewer(s);  
viewer(st,v);                          