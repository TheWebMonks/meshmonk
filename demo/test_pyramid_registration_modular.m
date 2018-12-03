%% Demo pyramid registration modular

% % Add MeshMonk's toolbox to the working path and setup current folder
% addpath(genpath('path\to\meshmonk')) % Set to location of meshmonk
% 
% studypath = 'path\to\DemoFolder\';   % Set to location of demo material
% cd(studypath);

%% Load face template (floating) and face (target)

floatingPath = [studypath '/Template.obj'];
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing
floatingPoints = single(floatingPoints);
floatingNormals = single(zeros(size(floatingPoints)));
% compute_normals(floatingPoints, floatingFaces, floatingNormals);
floatingFeatures = single([floatingPoints, floatingNormals]);
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

targetPath = [studypath '/demoFace.obj'];
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetPoints = single(targetPoints);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing
targetNormals = single(zeros(size(targetPoints)));
% compute_normals(targetPoints, targetFaces, targetNormals);
targetFeatures = single([targetPoints,targetNormals]); % WARNING: we had to flip the normals in this case!
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Prepare parameters and variables
%mex pyramid_registration.cpp -lmeshmonk

%# Set Parameters
numIterations = 90;
numPyramidLayers = 3;
downsampleFloatStart = 50;
downsampleFloatEnd = 0;
downsampleTargetStart = 70;
downsampleTargetEnd = 0;
correspondencesSymmetric = true;
correspondencesNumNeighbours = 3;
correspondencesFlagThreshold = 0.999;
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

floatingFeaturesTemp = [];
floatingFacesTemp = [];
floatingFlagsTemp = [];
floatingOriginalIndices = [];
floatingFeaturesTempOld = [];
floatingOriginalIndicesOld = [];

%% Execute Pyramid Registration via submodules

%# Iterative Registration process
for i=1:numPyramidLayers
    
    %##########################
    %##### DOWNSAMPLING #######
    %##########################
    
    %# Downsample Floating Mesh
    %## Determine downsample ratio
    downsampleRatio = downsampleFloatStart;
    if (numPyramidLayers > 1)
        downsampleRatio = single(round(downsampleFloatStart - (i-1) * round((downsampleFloatStart-downsampleFloatEnd)/(numPyramidLayers-1.0))));
    end
    downsampleRatio = downsampleRatio / 100.0;
    %## Downsample
    [ floatingFeaturesTemp, floatingFacesTemp, floatingFlagsTemp, floatingOriginalIndices ]...
        = downsample_mesh_clean( floatingFeatures, floatingFaces, floatingFlags, downsampleRatio );
    
    %# Downsample Target Mesh
    %## Determine downsample ratio
    downsampleRatio = downsampleTargetStart;
    if (numPyramidLayers > 1)
        downsampleRatio = single(round(downsampleTargetStart - (i-1) * round((downsampleTargetStart-downsampleTargetEnd)/(numPyramidLayers-1.0))));
    end
    downsampleRatio = downsampleRatio / 100.0;
    %## Downsample
    [ targetFeaturesTemp, targetFacesTemp, targetFlagsTemp ]...
        = downsample_mesh_clean( targetFeatures, targetFaces, targetFlags, downsampleRatio );
    
    %##########################
    %##### SCALESHIFTING ######
    %##########################
    if (i > 1)
        scaleshift_mesh(floatingFeaturesTempOld, floatingOriginalIndicesOld,...
                        floatingFeaturesTemp, floatingOriginalIndices);
    end
    
    %##########################
    %##### REGISTRATION #######
    %##########################
    
    for j=1:numIterations
        %# Compute Correspondences
        compute_correspondences(floatingFeaturesTemp, targetFeaturesTemp,...
                                floatingFlagsTemp, targetFlagsTemp,...
                                correspondingFeatures, correspondingFlags,...
                                correspondencesSymmetric, correspondencesNumNeighbours,...
                                correspondencesFlagThreshold, correspondencesEqualizePushPull);
        
        %# Compute Inlier Weights
        compute_inlier_weights(floatingFeaturesTemp, correspondingFeatures,...
                               correspondingFlags, inlierWeights,...
                               inlierKappa, inlierUseOrientation);
        
        %# Compute Transformation
        compute_nonrigid_transformation(floatingFeaturesTemp, correspondingFeatures,...
                                        floatingFacesTemp,...
                                        floatingFlagsTemp, inlierWeights,...
                                        10, transformSigma,...
                                        numViscousIterations, numElasticIterations);
        
        %# Annealing
        numViscousIterations = uint32(round(transformNumViscousIterationsStart * viscousAnnealingRate^(j)));
        if (numViscousIterations < transformNumViscousIterationsEnd) numViscousIterations = transformNumViscousIterationsEnd;end
        numElasticIterations = uint32(round(transformNumElasticIterationsStart * elasticAnnealingRate^(j)));
        if (numElasticIterations < transformNumElasticIterationsEnd) numElasticIterations = transformNumElasticIterationsEnd;end
    end
    
    %#########################################
    %##### COPY FOR NEXT PYRAMID LAYER #######
    %#########################################
    floatingFeaturesTempOld = floatingFeaturesTemp;
    floatingOriginalIndicesOld = floatingOriginalIndices;
    
end

%#############################
%##### FINAL SCALESHIFT ######
%#############################
originalIndices = uint32(0:numFloatingElements-1);
scaleshift_mesh(floatingFeaturesTemp, floatingOriginalIndices,...
                floatingFeatures, originalIndices);
            
clear correspondingFeatures correspondingFlags floatingFacesTemp floatingFeaturesTemp floatingFeaturesTempOld...
    floatingFlagsTemp floatingOriginalIndices floatingOriginalIndicesOld originalIndices targetFacesTemp targetFeaturesTemp...
    targetFlagsTemp 
%% Write Result
s = shape3D;
s.Vertices = correspondingFeatures(:,1:3);
st = shape3D;
st.Vertices = targetFeatures(:,1:3);
v = viewer(s);  
viewer(st,v);     