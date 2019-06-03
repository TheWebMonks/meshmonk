%%%%% In this tutorial we deal with two techniques for making sense of
% high dimensional shape information. 


clear all; close all;
addpath(genpath('/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/matlab')) % change this to the correct path on your computer



tutorialPath = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/tutorial';
DataPath = '/uz/data/avalok/mic/tmp/hmatth5/Projects/meshmonk/tutorial/TutorialData/SimulatedMappedFaces';
chdir(tutorialPath);

%% Load all faces
% here we load simulated mapped data. After mapping all faces are composed
% of the same number of points and the way the points are interconnected is
% the same, so they can all be loaded into a single 3D array


% load all .objs in datapath into an n vertices
objs = dir(strcat(DataPath,filesep,'*.obj'));

% remove any files prefaced with a '._' they are not real
isreal = ~contains({objs.name},'._');
objs = objs(isreal);



for i = 1:numel(objs)
   % import as a shape 3D
   obj = shape3D;
   obj.importWavefront(objs(i).name,objs(i).folder );
   if i==1
        % initialise 3D matrix to contain vertices of all faces 
        DataMatrix3D = zeros(obj.nVertices,3,numel(objs));
        

        
   end
   
   %put vertices into matrix
   DataMatrix3D(:,:,i) = obj.Vertices;
    
end



%% Generalized Procrustes Analysis
% to use either PCA or PLSR  all shapes must be centered on their
% mean, GPA is a technique for iteratively estimating the mean shape and
% centering all shapes on it

% initialise the mean shape with the Floating face used for registration
Template = shape3D;
Template.importWavefront('Template.obj',strcat(tutorialPath,filesep,'TutorialData'));
meanPoints = Template.Vertices;

iter = 3; % number of iterations to perform
scale = true; % scale all faces to a common size? the size to which they are scaled will be the size of the template (e.g. to scale to unit centrod size first scale the mean Points to this size)



for i = 1:iter
    for f = 1:numel(objs)
       faceVerts = DataMatrix3D(:,:,f);
       
       % align to to meanPoints
       T = computeTransform(faceVerts, meanPoints,scale);
       faceVerts = applyTransform(faceVerts,T);
       
       % put back in matris
       DataMatrix3D(:,:,f) = faceVerts;
    end
    
    % update mean
    meanPoints = mean(DataMatrix3D,3);
end

%% Visualise some faces (now aligned to to the mean)

nFaces = 5; % number of faces to visualise

% plot mean face
mFace = shape3D;
mFace.Vertices = meanPoints;
mFace.Faces = Template.Faces;
mFace.SingleColor = [1.,1.,1.];
v=viewer(mFace);
v.SceneLightVisible;
v.SceneLightLinked;
for f = 1:nFaces
   
    % create shape 3D
    shp = shape3D;
    shp.Vertices = DataMatrix3D(:,:,f);
    shp.Faces = Template.Faces;
   
    % add to existing viewer
    viewer(shp,v);
    shp.ViewMode = 'wireframe';
    % make semi transparent
    shp.Alpha = 0.3;    
end


%% 
%%%%%%%%%%%%%%%%%%%Principal Components Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Components Analysis is a technique for exploring the variation
% that is present in a sample of shapes.


%% Run PCA on shapes
close all;

% convert each case to vector representation so that each face is
% represented as arow of a 2D matrix

Y = arrayStructure2Vector(DataMatrix3D);

% mean center
Y0 = Y-mean(Y,1);

% Do PCA
[PCs, Scores,Variances,~,pctExplained] = pca(Y0);


%% Plot principal component as a vector field
% a principal component is essentially a weight on each landmark
% co-ordinate, one way to interpret this is as a linear transformation of
% facial shape. For example you can think of it as a vector field
close all;
% Plot a PC as a vector field
PCnum = 1; %which PC to plot

% convert weights from a single row vector to an N(vertices) x 3 array of
% vectors at each point

Transform = arrayVector2Structure(PCs(:,PCnum)');

% scale according to 3 standard deviations
Transform = Transform*sqrt(Variances(PCnum))*3;

% plot as vector field on average face

AverageShape = shape3D;
AverageShape.Vertices = meanPoints;
AverageShape.Faces = Template.Faces;

v = viewer(AverageShape);
vf = plotVectorField(meanPoints,Transform, v);
vf.SingleColor = [1.,1.,1.];


%% Plot Principal Component as a color-map
% an alternative libraray for creating color-maps in 'R' is:
% https://github.com/Arslan-Zaidi/Facial_masculinity_MHC/tree/master/Tutorial_visualizing_facial_heatmaps
% Also in Python using the Mayavi library https://docs.enthought.com/mayavi/mayavi/auto/mlab_helper_functions.html#mayavi.mlab.triangular_mesh
close all;
% Plot a PC as a vector field
PCnum = 1; %which PC to plot

% convert weights from a single row vector to an N(vertices) x 3 array of
% vectors at each point

Transform = arrayVector2Structure(PCs(:,PCnum)');

% scale according to 3 standard deviations
Transform = Transform*sqrt(Variances(PCnum))*3;

% create mesh of average face
AverageShape = shape3D;
AverageShape.Vertices = meanPoints;
AverageShape.Faces = Template.Faces;

% get component of transformation in direction of surface normals of
% average faces

normals = AverageShape.VertexNormals;
compNormDirection = sum(Transform.*normals,2);

% set average shape to be plotted with points color-indexed
AverageShape.VertexValue = compNormDirection;
AverageShape.ColorMode = 'indexed';
v = viewer(AverageShape);
v.BackgroundColor = [1,1,1];
v.SceneLightVisible = true;
v.SceneLightLinked = true;
colorbar();



%% Visualise principal component by transforming average face


PCnum = 1; %which PC to plot

% convert weights from a single row vector to an N(vertices) x 3 array of
% vectors at each point

Transform = arrayVector2Structure(PCs(:,PCnum)');
stdDev = sqrt(Variances(PCnum));

AverageShape = shape3D;
AverageShape.Vertices = meanPoints;
AverageShape.Faces = Template.Faces;

% add transform scaled according to  z-score of +/2 onto average vertices to
% visualise PC

PCminusMorph = clone(AverageShape);
PCminusMorph.Vertices = PCminusMorph.Vertices+Transform*stdDev*-2;

PCplusMorph = clone(AverageShape);
PCplusMorph.Vertices = PCplusMorph.Vertices+Transform*stdDev*2;


% plot each in separate viewer

v1 = viewer(PCminusMorph);
v1.SceneLightLinked = true;
v1.SceneLightVisible = true;

v2 = viewer(PCplusMorph);
v2.SceneLightLinked = true;
v2.SceneLightVisible = true;


%% Relationship between PC scores and the original shapes

% principal component scores re-code the shape onto a smaller number of
% variables. One way to think about this is that each PC is a particular
% linear transformation of facial shape and each shape is essentially a sum
% of these transformations. The amount of transformation along each PC,
% that will best reconstruct their shape is coded in their scores.

close all;
reconstructFaceNum = 1; % which face to reconstruct

FaceScores = Scores(reconstructFaceNum,:);

% plot the original face in one viewer
originalFace = shape3D;
originalFace.Vertices = DataMatrix3D(:,:,reconstructFaceNum);
originalFace.Faces = Template.Faces;

originalFace.SingleColor = [0.8,0.8,0.8];
v1 = viewer(originalFace);
v1.SceneLightVisible = true;
v1.SceneLightLinked = true;
% plot reconstruction - which begins at the sample average

reconstruction = shape3D;
reconstruction.Vertices = meanPoints;
reconstruction.Faces = Template.Faces;

reconstruction.SingleColor = [0.8,0.8,0.8];
v2 = viewer(reconstruction);
v2.SceneLightVisible = true;
v2.SceneLightLinked = true;
% iteratively add transformation along each PC to the reconstruction

nPCs = 10; % to save time we will only ad up to 10 PCs
faceScores = Scores(reconstructFaceNum,:);
for c = 1:nPCs
    Transform = arrayVector2Structure(PCs(:,c)');
    % scale by individual's scores
    Transform = Transform*faceScores(c);
    
    disp(strcat('To add transform corresponding of PC ',num2str(c),' press any key'));
    pause;
    reconstruction.Vertices = reconstruction.Vertices+Transform; % obj in viewer should update automatically

end



%% Create PC bi-plots coloured according to metadata

% load simulated metadata
close all;
metadata = readtable(strcat(tutorialPath,filesep,'TutorialData',filesep,'SimulatedMetadata.xlsx'));

PCxaxis = 1; % scores on which PC to plot on the x axis.
PCyaxis = 2; % scores on which PC to plot on the y axis

xaxis_scores = Scores(:,PCxaxis);
yaxis_scores = Scores(:,PCyaxis);


% plot bi-plot scatter with markers coloured according to age
f= figure;
scatter(xaxis_scores,yaxis_scores,[],metadata.Age,'filled')

% set equal aspect ratio
daspect([1,1,1])


xlabel(strcat('PC ',num2str(PCxaxis)));
ylabel(strcat('PC ',num2str(PCyaxis)));
title(gca,'Age Effect')
cb = colorbar();
title(cb, 'Age');


% plot colored by bmi
f = figure;
scatter(xaxis_scores,yaxis_scores,[],metadata.BMI,'filled')
title(gca,'BMI Effect')
% set equal aspect ratio
daspect([1,1,1])


xlabel(strcat('PC ',num2str(PCxaxis)));
ylabel(strcat('PC ',num2str(PCyaxis)));

cb = colorbar();
title(cb, 'BMI');




% plot by sex
f = figure;
hold on

m = scatter(xaxis_scores(metadata.Sex==1),yaxis_scores(metadata.Sex==1),[],'r','filled');
f = scatter(xaxis_scores(metadata.Sex==2),yaxis_scores(metadata.Sex==2),[],'b','filled');
daspect([1,1,1])
xlabel(strcat('PC ',num2str(PCxaxis)));
ylabel(strcat('PC ',num2str(PCyaxis)));
title(gca,'Sex Effect')
legend('Males','Females');
hold off


%% Fit PLS model

close all;

% use variables in metdata as predictors of shape
X = horzcat([metadata.Age,metadata.BMI,metadata.Sex]);
% convert each case to vector representation so that each face is
% represented as arow of a 2D matrix

Y = arrayStructure2Vector(DataMatrix3D);

%
[~,~,~,Yscores,Beta,pctvar] = plsregress(X,Y);

%% plot PLS Yscores 
% In the first instance, PLS is a dimension reduction  like PCA in this case it learns
% the subspace in shape that has the most covariance with the predictors

close all;

PLSxaxis = 1; % scores on which PLS component to plot on the x axis.
PLSyaxis = 2; % scores on which P to plot on the y axis

xaxis_scores = Yscores(:,PCxaxis);
yaxis_scores = Yscores(:,PCyaxis);
% plot bi-plot scatter with markers coloured according to age
f= figure;
scatter(xaxis_scores,yaxis_scores,[],metadata.Age, 'filled')
title(gca, 'Age Effect')
% set equal aspect ratio
daspect([1,1,1])


xlabel(strcat('PLS ',num2str(PLSxaxis)));
ylabel(strcat('PLS ',num2str(PLSyaxis)));

cb = colorbar();
title(cb, 'Age');


% plot colored by bmi
f = figure;
scatter(xaxis_scores,yaxis_scores,[],metadata.BMI,'filled')
title(gca,'BMI Effect')
% set equal aspect ratio
daspect([1,1,1])


xlabel(strcat('PLS ',num2str(PCxaxis)));
ylabel(strcat('PLS ',num2str(PCyaxis)));

cb = colorbar();
title(cb, 'BMI');




% plot by sex
f = figure;
hold on

m = scatter(xaxis_scores(metadata.Sex==1),yaxis_scores(metadata.Sex==1),[],'r','filled');
f = scatter(xaxis_scores(metadata.Sex==2),yaxis_scores(metadata.Sex==2),[],'b','filled');
daspect([1,1,1])
xlabel(strcat('PLS ',num2str(PCxaxis)));
ylabel(strcat('PLS ',num2str(PCyaxis)));

title(gca, 'Sex Effect')
legend('Males','Females');
hold off

%% Visualise independent effects of age, BMI and Sex on the face



ageEffect = arrayVector2Structure(Beta(2,:)); % the first row of Beta is the model intercept, the effect of the first variable is therefore in the second row of Beta
BMIEffect = arrayVector2Structure(Beta(3,:));
SexEffect = arrayVector2Structure(Beta(4,:));

%create shape3D of average shape
AverageShape = shape3D;
AverageShape.Vertices = meanPoints;
AverageShape.Faces = Template.Faces;




titles = {'Age','BMI','Sex'};
i = 0; % count from zero
for effect = {ageEffect,BMIEffect,SexEffect}
   i = i+1;
   obj = clone(AverageShape);
   v = viewer(obj);
   v.BackgroundColor = [1,1,1];
   v.SceneLightVisible = true;
   v.SceneLightLinked = true;
   v.Tag = titles{i};
   
   % calculate effect along he surface normals of the average
   normals = obj.VertexNormals;
   normComp = sum(effect{1}.*normals,2);
   obj.VertexValue = normComp;
   
   % view obj color-indexed to the effect along the surface normals
   obj.ColorMode = "Indexed";
   colorbar()
end

%% Permutation test on partial R-squared

%REFERENCE: Shrimpton et al (2014) https://doi.org/10.1016/j.forsciint.2013.10.021
close all;
varnum = 3; % the column of the Xblock, containing the predictor for which to test the independent effect
nperms = 1000; % the number of permutations with which to derive the emprirical sampling distribution

%get variable from Xblock
Var = X(:,varnum);

% get other variables
mask = 1:size(X,2)~=varnum;
CoVars = X(:,mask);

%%%%%%%% calculate partial Rsquared

% adjust for effect of CoVars on Y
[~,~,~,~,~,~,~,STATS]=plsregress(CoVars,Y);
reducedY=STATS.Yresiduals;

% adust for effect of CoVars on Var
[~,~,~,~,~,~,~,STATS]=plsregress(CoVars,Var);
reducedVar = STATS.Yresiduals;

% fit model regressing reducedY onto reducedVar
[~,~,~,~,~,PCTVAR]=plsregress(reducedVar,reducedY);
partRsqu = sum(PCTVAR(2,:));




%%%%%% Run permutation test
empiricalNullSamplingDistribution = zeros(1,nperms);
n = size(Var,1);

% seed random number generator
rng(1);

for i =1:nperms
    
    %permute rows of reducedVar
    permInds = randperm(n);
    permVar = reducedVar(permInds);
    
    
    % fit model regressing reducedY onto permuted data
    [~,~,~,~,~,PCTVAR]=plsregress(permVar,reducedY);
    empiricalNullSamplingDistribution(i) = sum(PCTVAR(2,:));
    
end

% compute p-value

p = sum(empiricalNullSamplingDistribution>partRsqu)/nperms;


% plot null sampling distribution
hold on
hist(empiricalNullSamplingDistribution)
title('Empirical Sampling Distribution of Partial R^2')
xlabel('Partial R^2')
ylabel('Frequency')

% plot observed value of partial R-squared as red line
h = plot([partRsqu,partRsqu],[0,nperms/2],'r');
disp(strcat('p=',num2str(p)));

legend(h,'Observed partial R^2'); 
