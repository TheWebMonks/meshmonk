clear all

%Load a mesh
floatingPath = '/home/jonatan/projects/meshmonk/examples/faceTemplate.obj';
floatingPath = '/home/jonatan/projects/meshmonk/examples/ExPeter/floating.obj';
%floatingPath = '/home/jonatan/projects/meshmonk/examples/data/bunny.obj';
[floatingPoints,floatingFaces] = read_vertices_and_faces_from_obj_file(floatingPath);
floatingFeatures = [floatingPoints, 1/sqrt(3.0)*ones(size(floatingPoints))];
floatingFeatures = single(floatingFeatures);
floatingFaces = uint32(floatingFaces-1); %-1 to make it compatible with C++ indexing ?
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));
clear floatingPoints;

%Load a mesh
targetPath = '/home/jonatan/projects/meshmonk/examples/faceTarget.obj';
targetPath = '/home/jonatan/projects/meshmonk/examples/ExPeter/target.obj';
%targetPath = '/home/jonatan/projects/meshmonk/examples/data/bunny2.obj';
[targetPoints,targetFaces] = read_vertices_and_faces_from_obj_file(targetPath);
targetFeatures = single([targetPoints, 1/sqrt(3.0)*ones(size(targetPoints))]);
targetFaces = uint32(targetFaces-1);%-1 to make it compatible with C++ indexing ?
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));
clear targetPoints;

%% Try the mexed pyramid_registration
%mex pyramid_registration.cpp -lmeshmonk

%# Set Parameters
numIterations = 60;
correspondencesSymmetric = true;
correspondencesNumNeighbours = 5;
inlierKappa = 4.0;
transformSigma = 3.0;
transformNumViscousIterationsStart = 50;
transformNumViscousIterationsEnd = 1;
transformNumElasticIterationsStart = 50;
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
                            correspondencesSymmetric, correspondencesNumNeighbours);
    
    %# Compute Inlier Weights
    compute_inlier_weights(floatingFeatures, correspondingFeatures,...
                           correspondingFlags, inlierWeights,...
                           inlierKappa);
    
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


%% Write Result
resultPath = '/home/jonatan/projects/meshmonk/examples/matlabResult.obj';
%resultPath = '/home/jonatan/projects/meshmonk/examples/data/bunnyResult.obj';
vertface2obj(floatingFeatures(:,1:3),floatingFaces,resultPath)
                            
                          