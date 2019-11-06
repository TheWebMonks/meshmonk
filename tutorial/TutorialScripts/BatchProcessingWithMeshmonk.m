%% This tutorial focuses on the practicalities of processing many files with Meshmonk
% Although we illustrate this with only one case this code will upscale to multiple cases

% before you begin you want each image in .obj format in a single folder. I
% assume this is called 'Targets'. You also want to have a few landmarks
% roughly indicated onto each target. 
% Each file in 'Targets' must have a Corresponding '.csv' file containing
% the co-ordinates of the landmarks. These should be in a separate folder
% called 'Landmarks'. You can indicate landmarks and export
% them using many software including the GeoMorph R package or MeshLab.

% The landmarks in all files must be IN
% THE SAME ORDER.
% If any mapping goes awry the landmarks are ALWAYS the first thing you
% should check.

% at present MeshMonk incorporates three steps:
% 1. initialisation based on landmarks
% 2. Rigid ICP
% 3. NonRigid ICP

% The purpose of the first two steps are to provide a good starting point
% for the thrid step. For more details on how the NonRigid ICP works see
% the tutorial 'UnderstandingNonRigidMappingWithMeshmonk'
%  In future an additional step will be implemented between 2 and 3 to
%  allow nonrigid registration constrained to by a statistical modeol of variation.

%% Set up paths

% ensure meshmonk is on your MATLAB path
clear all; close all;
addpath(genpath('/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/matlab')) % change this to the correct path on your computer

tutorialPath = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/tutorial/';
path2objs = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/tutorial/TutorialData/Targets/'; % check always that the path ends with a separator chara
path2landmarks = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/tutorial/TutorialData/Landmarks/';
path2results = strcat(tutorialPath,filesep,'MappedResults');
% get structs containing information for each file

%check if path2results exists and create it if not
if ~exist(path2results,'dir')
    mkdir(path2results)
end

objs = dir(strcat(path2objs,'*.obj'));
landmarks = dir(strcat(path2landmarks,'*.csv'));



if numel(objs)~=numel(landmarks)
    error('Check each landmark file has an equivalent obj file and vice versa')
end

%% Set up shape mappers
% The ShapeMappers handle the mapping according to the given transformation
% type. We define one for the rigid step and one for the non-rigid step
% for this demo we use the default settings we have found optimal for faces
% where the meshes are regular
% see supplementary material of White et al (2019): https://doi.org/10.1038%2Fs41598-019-42533-y

% Set up Rigid ICP step using ShapeMapper with 'rigid' transformation type
% in general this works fine and the settings don't need that much thinking
% about
RM= ShapeMapper;
RM.NumIterations = 30;

RM.TransformationType = 'rigid';
RM.UseScaling = true;

% settings for determining the correspondences
RM.CorrespondencesNumNeighbours = 3; % number of neighbours to use to estimate correspondences
RM.CorrespondencesFlagThreshold = 0.9; 
RM.CorrespondencesSymmetric = true; % if true correspondences are estimated from the template to target and target to template and combined - 
                                    % this can help with mapping structures
                                    % such as long bones and allows the
                                    % target to 'pull' the template
RM.CorrespondencesEqualizePushPull = false;

% settings that determine which points are 'outliers' not used to estimate the
% transformation. This is based on 1. whether the point is an outlier in the
% distribution of differences between floating and target and 2. Whether
% the points corresponds to a point that has been 'flagged'.
RM.InlierKappa = 3;
RM.InlierUseOrientation = true; % use surface normal direction to determine of point is inlier/outlier

% ignore points that correspond to the edges of the mesh - only applicable
% to 'open' surfaces like the face
RM.FlagFloatingBoundary = true; % ignore points that correspondences to the edge of Floating surface - only applicable if using SymmetricCorrespndences
RM.FlagTargetBoundary = true;% ignore points that correspondences to the edge of Target surface

RM.FlagTargetBadlySizedTriangles = true; % ignore points that match to regions of abnormally sized triangles
RM.TriangleSizeZscore = 6; % threshold to determine which triangles are abnormally sized

RM.UpSampleTarget = false; % will upsample the target mesh. if meshes are irregular it can help to set this to true
  

% Set up non rigid ICP with ShapeMapper with 'nonrigid' transformation type
NRM = ShapeMapper;
NRM.TransformationType = 'nonrigid';
NRM.NumIterations = 200; 
NRM.CorrespondencesSymmetric = true;
NRM.CorrespondencesNumNeighbours = 3;
NRM.CorrespondencesFlagThreshold = 0.9;
NRM.CorrespondencesUseOrientation = true;
NRM.CorrespondencesEqualizePushPull =false;
NRM.InlierKappa = 12; % basically you want this really high, so that the outlier detection is effectively turned off. in future you will be able to actually turn this off. 
NRM.InlierUseOrientation = true;
NRM.FlagFloatingBoundary = true;
NRM.FlagTargetBoundary = true;
NRM.FlagTargetBadlySizedTriangles = true;
NRM.TriangleSizeZscore = 6;
NRM.UpSampleTarget = false;




% parameters specific to non-rigid step. In general, if computing time is not a factor
% set the NumIterations and the TransformNumViscousIterationsStart and
% TransformNumElasticIterationsStart as high as possible. This will ensure
% a slow, gradual deformation of the template to target.
% please see the tutorial 'UnderstandingNonRigidMappingWithMeshmonk' for
% some more details on how this all works. 
% TransformNumViscousIterationsEnd, TransformNumElasticIterationsEnd should
% always be around 1. If they are higher than the floating shape may not
% actually match the shape of the target at the end of the algorithm. They
% can be 0 but this is only advisable with very high quality meshes


NRM.TransformSigma = 3;
NRM.TransformNumViscousIterationsStart = 200;
NRM.TransformNumViscousIterationsEnd = 1;
NRM.TransformNumElasticIterationsStart = 200;
NRM.TransformNumElasticIterationsEnd = 1;
NRM.TransformNumNeighbors = 80;


%% Batch process files

% load Template floating face
Floating = shape3D;
Floating.importWavefront('Template.obj',strcat(tutorialPath, filesep,'TutorialData/'));

% load landmarks on template
FloatingLandmarks = readTextLandmarkFile(strcat(tutorialPath,filesep,'TutorialData/Template.csv'),',');



for f = 1:numel(objs)
   % load Target obj and landmarks
   targetName = objs(f).name;
   Target = shape3D;
   Target.importWavefront(targetName, path2objs);
   
   TargetLandmarks = readTextLandmarkFile(strcat(path2landmarks,targetName(1:end-4),'.csv'),',');
   
   
   %%%% Initialisation from landmarks
   
   % calculate rigid transformation from FloatingLandmarks to
   % TargetLandmarks
   T = computeTransform(FloatingLandmarks,TargetLandmarks,true);
   
   % apply rigid transform to Floating
   forFloating = applyTransform(Floating,T);
   
   
   %check normals of floating and target consistently point inward or
   %outward. If they are inconsistent the normals will be flipped. 
   %THIS CHECK WORKS IN MOST CASES BUT IS NOT INFALLIBLE. PARTICULARLY
   % IT MAY FAIL IF THE INITIALISATION IS BAD OR FLOATING AND TARGET ARE
   % VERY DIFFERENT SHAPES
   % The only way to exactly check this is to plot
   % each shape with their normals
   % e.g v = viewer(Shape)
   %    plotVectorField(Shape.Vertices,Shape.VertexNormals,v)
   
   if ~normalsConsistent(Target, forFloating)
       Target.FlipNormals = true;
   end
   
   
   % Execute Rigid Mapping
   forRM = clone(RM);
   forRM.FloatingShape = clone(forFloating);
   forRM.TargetShape = clone(Target);
   
   forRM.map();
   
   % Execute Non-Rigid Mapping
   forNRM = clone(NRM);
   forNRM.FloatingShape = forRM.FloatingShape; % floating shape is now the floating shaoe after Rigid Mapping
   forNRM.TargetShape = Target;
   forNRM.map()
   
   
   %save result
   forNRM.FloatingShape.exportWavefront(targetName,path2results);
   
   
  
   
   
   
   

   
   
   
   
    
    
    
    
end















