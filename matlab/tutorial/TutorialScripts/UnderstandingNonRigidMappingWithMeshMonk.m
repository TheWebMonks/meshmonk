%% This tutorial is designed to give an understanding of nonprigid mapping 
% works in MeshMonk. Other tutorials deal more with how to use it in
% practice


tutorialPath = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/matlab/tutorial';
chdir(tutorialPath)

% add meshmonk classes to your MATLAB path

addpath(genpath('/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/matlab'))

%%
%%%%%%%%%%%%%%%%%%%%%%%  Intro to Shape3D and Viewer class %%%%%%%%%%%%%%%%%%%%%%%%%%
% Shape3D represents a mesh object. Defining properties of a mesh are the
% 'Vertices' of which it is composed, and the 'Faces' that link the points
% together into a surface. 
% 
% The normals of each vertex are derived from these properties and can also
% used in the 'MeshMonk' mapping.
%
% The 'Viewer' class is a class we provide for visualising meshes. 
% Many commercial or open software support some or all of the same
% functionality - for example 'MeshLab', 'Blender', etc...
% also the matlab function 'trisurf' can be used to plot a mesh


%% load an .obj wavefront as a Shape3D

obj = shape3D; % create instance of the class
obj.importWavefront('Template.obj','TutorialData');

%% visualise vertices using Viewer
close all;
viewer(obj);
obj.ViewMode = 'points';


%% Visualise the triangulated surface
close all;
viewer(obj);
obj.ViewMode = 'wireframe';

%% Visualise the smooth surface
close all;
v = viewer(obj);
obj.ViewMode = 'solid';
v.SceneLightVisible = true;
v.SceneLightLinked = true;


%% Visualise the surface normals
% surface normals are the direction perpendicular the surface at that
% point they convey information about the local orientation of the surface
% and are used to guide the mapping 

close all;

v = viewer(obj);
plotVectorField(obj.Vertices, obj.VertexNormals*2,v); % times vertex normals by two to make them big enough to see
v.SceneLightVisible = true;
v.SceneLightLinked = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Meshmonk nuts and bolts %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The next few sections of the tutorial ae designed to give an intuitive
% feel for how meshmonk really works under the hood. We won't cover
% everything but will cover the most important aspects of the algorithm.
% In this section we will make use of some of the core Meshmonk 'mex'
% functions directly. In practice you do not have to do this and all
% functionality is available through the more user-friendly 'ShapeMapper'
% class.


%The goal of meshmonk to effectively warp the template (floating) image
%into the shape of the target image. This is accomplished in two phases, a
%rigid-regisitration step and a non-rigid registration step. The purpose of the
%first step is to roughly align the the template to the target. In general
%this works adequately and requires little input or adjustment. Therefore
%when discussing the nuts and bolts, we focu only on the non-rigid step.

% Each iteration of the non-rigid step:
% 1. Estimates corresponding points 
% 2. Estimates and applies the non-rigid deformation from the template to the target


% We now illustrate each of these steps, but first need to assemble the
% data correctly (this is usually handled by the ShapeMapper)




%% load Floating and Target Meshes

close all; clear all;

floating = shape3D;
floating.importWavefront('Template_Rigidly_Aligned.obj','TutorialData');

target = shape3D;
target.importWavefront('demoFace.obj','TutorialData/Targets');

obj = ShapeMapper;
obj.FloatingShape = clone(floating);
obj.TargetShape = clone(target);



% assemble arrays of features - these are necessary input to the core
% meshmonk functions. In practice however this will be handled by the
% ShapeMapper class
floatingFeatures = single([floating.Vertices, floating.VertexNormals]);
floatingFaces = uint32(floating.Faces-1);
numFloatingElements = size(floatingFeatures,1);
floatingFlags = single(ones(numFloatingElements,1));

targetFeatures = single([target.Vertices, target.VertexNormals]);
targetFaces = uint32(floating.Faces-1);
numTargetElements = size(targetFeatures,1);
targetFlags = single(ones(numTargetElements,1));


%% Step 1 - compute correspondences
close all

% some settings - feel free to play with them and see what happens
correspondencesSymmetric = false; 
correspondencesNumNeighbours = 3;
correspondencesFlagThreshold = 0.9;
correspondencesEqualizePushPull = false;

% initialise array to contain corresponding features - this will be
% modified in situ when 'compute_correspondences' is run
correspondingFeatures = single(zeros(numFloatingElements,6));
correspondingFlags = single(ones(numFloatingElements,1));
inlierWeights = single(ones(numFloatingElements,1));

compute_correspondences(floatingFeatures, targetFeatures,...
                            floatingFlags, targetFlags,...
                            correspondingFeatures, correspondingFlags,...
                            correspondencesSymmetric, correspondencesNumNeighbours,...
                            correspondencesFlagThreshold, correspondencesEqualizePushPull);

                        
                        
% Visualise correspondences from template to target

v = viewer(floating);
v.SceneLightVisible = true;
v.SceneLightLinked;
v.BackgroundColor = [1.,1.,1.];
floating.SingleColor = [0.,0.,1];
floating.ViewMode = 'wireframe';
floating.Alpha = 0.5;
viewer(target,v);
target.Alpha = 1.;
target.SingleColor = [0.,0.9,0.9];

% plot vectors from temlate to corresponding points -  they provide an initial
% (but not very good) estimate of the deformation field;

% plot the deformation vectors in red
deformationVectors = double(correspondingFeatures(:,1:3))-floating.Vertices;
obj = plotVectorField(floating.Vertices,deformationVectors, v);
obj.SingleColor = [1.,.2,.2];

%% See what happens if you apply the initial estimate directly

% the correspondences give an idea of where the mesh needs to go, but we
% cannot deform straight to the corresponding points - look at what would
% happen to the floating face by executing this cell.
% You can see visually the results are not
% very smooth and the mesh is compressed and expanded haphazardly
% indicating generally poor correspondences
close all
deformedFloating = shape3D;
deformedFloating.Vertices = floating.Vertices + deformationVectors; % deform floating mesh along deformation vectors
deformedFloating.Faces = floating.Faces;

% plot solid mesh 
v = viewer(deformedFloating);
deformedFloating.ViewMode = 'Solid';
deformedFloating.SingleColor = [1.,1.,1.];
deformedFloating.Alpha = 1.

% overlay with wireframe
wireDeformedFloating = clone(deformedFloating);
wireDeformedFloating.SingleColor = [0,0,0];
wireDeformedFloating.ViewMode = 'wireframe';
viewer(wireDeformedFloating,v);



%% Step 2 - Regularise deformation Field
% So in this next step the deformation field is smoothed and then applied to the
% floating features


%%%%%Exercise: Play with the regularisation parameters and see what happens

% to reduce regularisation reduce the values of 'numViscousIterations' and
% 'numElasticIterations'. Try setting both to 200, 25, and 0 and see what
% happens
close all

% settings controlling the amount of regularisation
numViscousIterations = 0;
numElasticIterations = 0;
transformSigma = 3.;
transformNumNeighbours = 80;

% for the purposes of this demo we create a copy of floatingFeatures from the original data within the cell so
% that the cell can be run multiple times without successive deformations  being added
% to each other

floatingFeatures = single([floating.Vertices, floating.VertexNormals]);

compute_nonrigid_transformation(floatingFeatures, correspondingFeatures,...
                                    floatingFaces,...
                                    floatingFlags, inlierWeights,...
                                    transformNumNeighbours, transformSigma,...
                                    numViscousIterations, numElasticIterations);

% create shape3D of the deformed floating mesh
deformedFloating = shape3D;
deformedFloating.Vertices = floatingFeatures(:,1:3);
deformedFloating.Faces = floating.Faces;
deformedFloating.SingleColor=[1.,1.,1,];

v = viewer(deformedFloating);

% overlay with wireframe
wireDeformedFloating = clone(deformedFloating);
wireDeformedFloating.SingleColor = [0,0,0];
wireDeformedFloating.ViewMode = 'wireframe';
viewer(wireDeformedFloating,v);

% plot target in same window

viewer(target,v)
target.Alpha = 0.5;
v.SceneLightVisible = true;
V.SceneLightLinked;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Gradual Deterministic Annealing  %%%%%%%%


% You will have noticed in the last section that with an appropriate amount
% of regularisation the flaoting face can be deformed towards the target,
% but that with too little the results go 'pear-shaped'. They go
% pear-shaped because the correspondences are not very accurate. 
% the key to the successful application of meshmonk is to let the mesh
% deform gradually over many iterations. After each iteration, the
% correspondences will be incrementally better, meaning less smoothing can
% be applied on the next iteration. 


%% plot deterministic annealing
% how the regularisation (number of Viscous or Elastic Smoothing Iterations) is decreased o
% over multiple iterations of the algorithm is 'deterministic'. 
% this means it is directly a function of the numer of iterations of the algorithm iteslf and user-specified
% values of how much regularisation to start and end with. 


% In practice to use meshmonk you will use the 'ShapeMapper class' in this
% class, this process is controlled by four parameters:
close all;
NumIterations = 80; % number of iterations of the algorithm to execute

%Amount of regularisation at iteration one of the algorithm
TransformNumViscousIterationsStart = 200; 
TransformNumElasticIterationsStart = 200;
% Amount of regularisation at the end of the algorithm
TransformNumViscousIterationsEnd = 1;
TransformNumElasticIterationsEnd = 1;


% below plots the amount of regularisation (number of Viscous smoothing iterations) at each iteration of the
% algorithm,
% the function would be the same for the elastic smoothing

annealing_rate = exp(log(TransformNumViscousIterationsEnd/TransformNumViscousIterationsStart)/(NumIterations-1));
ViscousIterations = round(TransformNumViscousIterationsStart*annealing_rate.^(0:NumIterations-1));

plot(1:NumIterations,ViscousIterations);
xlabel('Algorithm Iteration')
ylabel('Amount of Smoothing Iterations')
title('Deterministic Annealing')
%%
%%%%%%%%%%%%%%%%%% Introduction to the Shape Mapper  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the whhole process is wrapped by a class called 'ShapeMapper'. This
% interacts with the core meshmonk functions and gives the user access to
% the important settings

% below we create an isntance of the ShapeMapper for nonrigid ICP and give
% the default settings for faces



SM = ShapeMapper;
SM.TransformationType = 'nonrigid';


% Settings for computing correspondences on each iteration
SM.CorrespondencesSymmetric = true;% if true correspondences are estimated from the template to target and target to template and combined - 
                                    % this can help with mapping structures such as long bones and allows the target to 'pull' the template
SM.CorrespondencesNumNeighbours = 3;% number of neighbours to use to estimate correspondences
SM.CorrespondencesUseOrientation = true; % whether to use the surface normal direction to estimate correspondences
SM.CorrespondencesEqualizePushPull =false; %






% Meshmonk will also ignore, or downweight the influence of certain points
% on the computed deformation field at each iteration. These settings
% control how these points are determined 

SM.InlierKappa = 12; % if the point is extreme within the distribution of differences between points of the floating and target face it will be ignored. This setting is z-score determining how extreme this should be  
SM.InlierUseOrientation = true; % whether to use similarity of surface normal direction is determining the 'difference' between points on the floting and target, if false this will just be the sdistnce between corresponding points. 
SM.FlagFloatingBoundary = true; % whether to ignore points that correspond to a point on the edge of the floating face
SM.FlagTargetBoundary = true; % whether to ignore points that correspond to a point on the edge of the floating face
SM.FlagTargetBadlySizedTriangles = true; % whether to ignore points that correspond to a point belonging to a traingle that is unusual in the context of the given mesh 
SM.TriangleSizeZscore = 6; % how unusual does the triangle (see above) have to be?
SM.CorrespondencesFlagThreshold = 0.9;





%%%%%%%Settings controlling deterministic annealing

% In general, if computing time is not a factor
% set the NumIterations and the TransformNumViscousIterationsStart and
% TransformNumElasticIterationsStart as high as possible. This will ensure
% a slow, gradual deformation of the template to target.
% please see the tutorial 'UnderstandingNonRigidMappingWithMeshmonk' for
% some more details on how this all works. 
% TransformNumViscousIterationsEnd, TransformNumElasticIterationsEnd should
% always be around 1. If they are higher than the floating shape may not
% actually match the shape of the target at the end of the algorithm. They
% can be 0 but this is only advisable with very high quality meshes

SM.NumIterations = 200; % number of iterations of the algorithm
SM.TransformSigma = 3; % the width of the Gaussian with which to convolve the deformation field
SM.TransformNumNeighbors = 80; % for efficiency, when convolving with the Gaussian not all deformation vectors are incorporated into the wieghted average
                                % it will only average the n closest neighbours. This setting sets the value of n
SM.TransformNumViscousIterationsStart = 200; % number of iterations of viscous smoothing applied during iteration 1
SM.TransformNumViscousIterationsEnd = 1;% number of iterations of viscous smoothing applied during iteration = numIterations
SM.TransformNumElasticIterationsStart = 200;% number of iterations of elastic smoothing applied during iteration 1
SM.TransformNumElasticIterationsEnd = 1;% number of iterations of elastic smoothing applied during iteration = numIterations

%%% settings controlling upsampling

%when meshes are irregular it can be useful to upsample the target image
%during the mapping. 

SM.UpSampleTarget = false; % whether to upsample the target or not
UpSampleMode = 'size'; %
UpSampleVal = 1.5; %





%%%%%Run and display mapping

%set floating and target face
SM.FloatingShape = clone(floating);
SM.TargetShape = clone(target);
SM.TargetShape.ViewMode = 'Solid';


% execute  and visualise mapping
SM.Display = true;
SM.map();




















% these basic meshmonk functions are a bit temperamental about the input
% they take so there needs to be some jiggery pokery with data types 
% (in parctice using the ShapeMapper, this will be handled for you. )

FloatingFeatures = single([floating.Vertices, floating.VertexNormals]);
FloatingFaces = uint32(floating.Faces-1);

TargetFeatures = single([target.Vertices, target.VertexNormals]);
TargetFaces = uint32(floating.Faces-1);






























%%