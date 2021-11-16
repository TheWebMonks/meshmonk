classdef ShapeMapper < superHandleClass
   properties
       FloatingShape;% floating shape (template)
       TargetShape;% target shape
       FloatingFlags;% Deterministic outlier flagging on floating shape
       TargetFlags;% Deterministic outlier flagging on target shape
   end
   properties (Dependent = true)
       nFloatingFeatures;
       nTargetFeatures;
   end
   properties (Hidden = true)
       CorrespondingShape;% Shape object that holds the correspondences during mapping
       CorrespondingFeatures;% Correspondences, location + normals;
       CorrespondingFlags;% Determinstic outlier flagging on correspondences
       InlierWeights;% A weight per vertex, 0 = outlier; 1 = inlier
   end
   properties (Hidden = true, Dependent = true)
       FloatingPoints;
       FloatingFaces;
       FloatingNormals;
       FloatingFeatures;
       TargetPoints;
       TargetFaces;
       TargetNormals;
       TargetFeatures;
       CorrespondingPoints;
       CorrespondingNormals;
   end 
   properties % parameter settings
       NumIterations = 60;
       CorrespondencesSymmetric = true;
       CorrespondencesNumNeighbours = 3;
       CorrespondencesFlagThreshold = 0.9;
       CorrespondencesUseOrientation = true;
       CorrespondencesEqualizePushPull = false;
       InlierKappa = 12;
       InlierUseOrientation = true;
       FlagFloatingBoundary = true;
       FlagTargetBoundary = true;
       UseInlierWeights = true;
       FlagTargetBadlySizedTriangles = true;
       TriangleSizeZscore = 6;% maximum zscore allowed before triangles are flagged;
       UpSampleTarget = false;
       UpSampleMode = 'size';
       UpSampleVal = 1.5;
       UseScaling = true;
       TransformationType = 'rigid';% rigid (ICP) or non rigid (visco-elastic) or PCA (morphable (statistical) model) 
       TransformSigma = 3.0;
       TransformNumViscousIterationsStart = 50;
       TransformNumViscousIterationsEnd = 1;
       TransformNumElasticIterationsStart = 50;
       TransformNumElasticIterationsEnd = 1;
       TransformNumNeighbors = 32;
       TransformMorphableModel = [];
       TransformMorphableModelNetaFraction = 100;
       Display = false;
   end
   properties (Dependent = true)
   end
   properties (Hidden = true)
       Iteration = 0;
       ElapsedTime = 0;
       TransformationMatrix = single(eye(4,4));
   end
   properties (Hidden = true, Dependent = true)
       TransformMorphableModelNetaStart;
       TransformMorphableModelNetaEnd;
       ViscousAnnealingRate;
       ElasticAnnealingRate;
       NetaAnnealingRate;
       NumViscousIterations;
       NumElasticIterations;
       NetaValue;
   end
   methods %CONSTRUCTOR
       function obj = ShapeMapper(varargin)
          obj = obj@superHandleClass(varargin{:});
       end % 3D viewer Constructor       
   end
   methods % GETTING
       function out = get.FloatingShape(obj)
           out = obj.FloatingShape;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.FloatingPoints(obj)
           if isempty(obj.FloatingShape), out = []; return; end
           out = single(obj.FloatingShape.Vertices);
       end
       function out = get.FloatingFaces(obj)
           if isempty(obj.FloatingShape), out = []; return; end
           out = uint32(obj.FloatingShape.Faces-1);
       end
       function out = get.FloatingNormals(obj)
           if isempty(obj.FloatingShape), out = []; return; end
           out = single(obj.FloatingShape.VertexNormals);
       end
       function out = get.FloatingFeatures(obj)
           if isempty(obj.FloatingShape), out = []; return; end
           if obj.CorrespondencesUseOrientation, out = [obj.FloatingPoints, obj.FloatingNormals]; return; end
           out = [obj.FloatingPoints, zeros(obj.FloatingShape.nVertices,3)];
       end
       function out = get.TargetShape(obj)
           out = obj.TargetShape;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.TargetPoints(obj)
           if isempty(obj.TargetShape), out = []; return; end
           out = single(obj.TargetShape.Vertices);
       end
       function out = get.TargetFaces(obj)
           if isempty(obj.TargetShape), out = []; return; end
           out = uint32(obj.TargetShape.Faces-1);
       end
       function out = get.TargetNormals(obj)
           if isempty(obj.TargetShape), out = []; return; end
           out = single(obj.TargetShape.VertexNormals);
       end
       function out = get.TargetFeatures(obj)
           if isempty(obj.TargetShape), out = []; return; end
           if obj.CorrespondencesUseOrientation, out = [obj.TargetPoints, obj.TargetNormals]; return; end
           out = [obj.TargetPoints, zeros(obj.TargetShape.nVertices,3)];
       end
       function out = get.FloatingFlags(obj)
                if isempty(obj.FloatingShape), out = []; return; end
                out = single(obj.FloatingFlags);
                if isempty(out), out = single(ones(obj.FloatingShape.nVertices,1)); end
                if obj.FlagFloatingBoundary, out(getBoundary(obj.FloatingShape)) = single(0); end
       end
       function out = get.TargetFlags(obj)
                if isempty(obj.TargetShape), out = []; return; end
                out = single(obj.TargetFlags);
                if isempty(out), out = single(ones(obj.TargetShape.nVertices,1)); end
                if obj.FlagTargetBoundary, out(getBoundary(obj.TargetShape)) = single(0); end
                if obj.FlagTargetBadlySizedTriangles
                   [~,index] = getBadlySizedTriangles(obj.TargetShape,obj.TriangleSizeZscore);
                   out(index) = single(0);
                end
       end
       function out = get.ViscousAnnealingRate(obj)
           out = exp(log(single(obj.TransformNumViscousIterationsEnd)/single(obj.TransformNumViscousIterationsStart))/(obj.NumIterations-1));
       end
       function out = get.ElasticAnnealingRate(obj)
           out = exp(log(single(obj.TransformNumElasticIterationsEnd)/single(obj.TransformNumElasticIterationsStart))/(obj.NumIterations-1));
       end
       function out = get.NetaAnnealingRate(obj)
           out = exp(log(single(obj.TransformMorphableModelNetaEnd)/single(obj.TransformMorphableModelNetaStart))/(obj.NumIterations-1));
       end
       function out = get.NumViscousIterations(obj)
           out = uint32(round(obj.TransformNumViscousIterationsStart * obj.ViscousAnnealingRate^(obj.Iteration-1)));
           if (out < obj.TransformNumViscousIterationsEnd), out = obj.TransformNumViscousIterationsEnd;end
       end
       function out = get.NumElasticIterations(obj)
           out = uint32(round(obj.TransformNumElasticIterationsStart * obj.ElasticAnnealingRate^(obj.Iteration-1)));
           if (out < obj.TransformNumElasticIterationsEnd), out = transformNumElasticIterationsEnd;end
       end
       function out = get.NetaValue(obj)
           out = round(obj.TransformMorphableModelNetaStart * obj.NetaAnnealingRate^(obj.Iteration));
           if (out < obj.TransformMorphableModelNetaEnd), out = obj.TransformMorphableModelNetaEnd;end
           out = out./obj.TransformMorphableModelNetaStart;
       end
       function out = get.TransformMorphableModelNetaEnd(obj) %#ok<MANU>
          out = 1; 
       end
       function out = get.TransformMorphableModelNetaStart(obj)
          out = obj.TransformMorphableModelNetaFraction; 
       end
       function out = get.CorrespondingShape(obj)
           out = obj.CorrespondingShape;
           if ~superHandleClass.isH(out), out = []; end
       end
       function out = get.nFloatingFeatures(obj)
           if isempty(obj.FloatingShape), out = 0; return; end
           out = obj.FloatingShape.nVertices;
       end
       function out = get.nTargetFeatures(obj)
           if isempty(obj.TargetShape), out = 0; return; end
           out = obj.TargetShape.nVertices;
       end
       function out = get.CorrespondingPoints(obj)
          if isempty(obj.CorrespondingFeatures), out = []; return; end
          out = obj.CorrespondingFeatures(:,1:3);
       end
       function out = get.CorrespondingNormals(obj)
           if isempty(obj.CorrespondingFeatures), out = []; return; end
           out = obj.CorrespondingFeatures(:,4:6);
       end
   end
   methods % SETTING
       function obj = set.FloatingPoints(obj,in)
          if isempty(obj.FloatingShape),return;end 
          obj.FloatingShape.Vertices = double(in);
       end
       function obj = set.FloatingFaces(obj,in)
          if isempty(obj.FloatingShape),return;end 
          obj.FloatingShape.Faces = double(in); 
       end
       function obj = set.FloatingNormals(obj,in) %#ok<*INUSD>
          return;% normals are never set 
       end
       function obj = set.FloatingFeatures(obj,in)
          if isempty(obj.FloatingShape),return;end 
          obj.floatingPoints = in(:,1:3);
       end
       function obj = set.TargetPoints(obj,in)
          if isempty(obj.TargetShape),return;end 
          obj.TargetShape.Vertices = double(in);
       end
       function obj = set.TargetFaces(obj,in)
          if isempty(obj.TargetShape),return;end 
          obj.TargetShape.Faces = double(in); 
       end
       function obj = set.TargetNormals(obj,in) %#ok<*INUSD>
          return;% normals are never set 
       end
       function obj = set.TargetFeatures(obj,in)
          if isempty(obj.TargetShape),return;end 
          obj.targetPoints = in(:,1:3);
       end
       function obj = set.CorrespondingFeatures(obj,in)
          obj.CorrespondingFeatures = in;
          if~isempty(obj.CorrespondingShape), obj.CorrespondingShape.Vertices = in(:,1:3);end %#ok<*MCSUP>
       end
       function obj = set.TransformationType(obj,in)
            if ismember(in,{'nonrigid','rigid'})
                obj.TransformationType = in;
            elseif ismember(in,{'pca'})
                error('TransformationType not yet implemented')
            else
                error('Invalid TransformationType')
            end
       end 
   end
   methods % INTERFACING
       function map(obj)
           % INITIALIZATION
           floatingFeatures = obj.FloatingFeatures;
           floatingFlags = obj.FloatingFlags;
           
           targetFeatures = obj.TargetFeatures;
           targetFlags = obj.TargetFlags;
           
           if obj.UpSampleTarget
              TargetBck = clone(obj.TargetShape);
              obj.TargetShape.VertexValue = targetFlags;
              upSampleMesh(obj.TargetShape,obj.UpSampleMode,obj.UpSampleVal);
              targetFeatures = obj.TargetFeatures;
              targetFlags = single(obj.TargetShape.VertexValue);
              targetFlags(targetFlags<obj.CorrespondencesFlagThreshold) = single(0);
           end
              
           obj.CorrespondingFeatures = single(zeros(obj.nFloatingFeatures,6));
           obj.CorrespondingFlags = single(ones(obj.nFloatingFeatures,1));
           obj.InlierWeights = single(ones(obj.nFloatingFeatures,1));
           
           obj.CorrespondingShape = clone(obj.FloatingShape);
           obj.CorrespondingShape.VertexValue = floatingFlags;
           
           if obj.Display % Setting up viewers and then press enter to continue;
              v = viewer(obj.FloatingShape);
              v.Figure.Position = [150 24 100 30];
              obj.TargetShape.RenderAxes = v.RenderAxes;
              obj.TargetShape.Visible = true;
              obj.TargetShape.ViewMode = 'Wireframe';
              obj.TargetShape.SingleColor = [0.8 0.8 0.8];
              obj.FloatingShape.ViewMode = 'Wireframe';
              obj.FloatingShape.SingleColor = [0.1 0.5 0.9];
              v2 = viewer(obj.CorrespondingShape);
              v2.Figure.Position = [40 24 100 30];
              obj.CorrespondingShape.ColorMode = 'Indexed';
              obj.CorrespondingShape.ViewMode = 'Wireframe';
              str = input('Are you ready? Y/N [Y]: ','s');
              if isempty(str), str = 'y'; end
              switch lower(str)
                  case 'y'
                      disp('Started...');
                  case 'n'
                      return;
              end
               h=waitbar(0,'Shape Mapping...');
           end
           
           obj.Iteration = 0;
           while obj.Iteration <= obj.NumIterations
               % ESTABLISH CORRESPONDENCES
               compute_correspondences(floatingFeatures, targetFeatures,...
                        floatingFlags, targetFlags,...
                        obj.CorrespondingFeatures, obj.CorrespondingFlags,...
                        obj.CorrespondencesSymmetric, obj.CorrespondencesNumNeighbours,...
                        obj.CorrespondencesFlagThreshold,obj.CorrespondencesEqualizePushPull);
               % UPDATE INLIER WEIGHTS
               compute_inlier_weights(floatingFeatures, obj.CorrespondingFeatures,...
                           obj.CorrespondingFlags, obj.InlierWeights,...
                           obj.InlierKappa, obj.InlierUseOrientation);
                           
               if ~obj.UseInlierWeights % if not using inlier weights - still make use of flags
                   %set any points that are not flagged to have
                   %inlierWeights of 1
                   obj.InlierWeights(obj.CorrespondingFlags==1) = 1;
                   
               end               
               if obj.Display||strcmpi(obj.TransformationType,'pca')
                  obj.CorrespondingShape.Vertices = obj.CorrespondingFeatures(:,1:3);
                  obj.CorrespondingShape.VertexValue = obj.InlierWeights;
                  if obj.Display, drawnow expose;end
               end
               % UPDATE TRANSFORMATION
               switch lower(obj.TransformationType)
                   case 'rigid'
                       compute_rigid_transformation(floatingFeatures, obj.CorrespondingFeatures,...
                                   obj.InlierWeights, obj.TransformationMatrix, obj.UseScaling);
                   case 'nonrigid'
                       compute_nonrigid_transformation(floatingFeatures, obj.CorrespondingFeatures,...
                                    obj.FloatingFaces,...
                                    floatingFlags, obj.InlierWeights,...
                                    obj.TransformNumNeighbors, obj.TransformSigma,...
                                    obj.NumViscousIterations, obj.NumElasticIterations);
                   case 'pca'
                       obj.TransformMorphableModel.Neta = double(obj.NetaValue);
                       update = compute_morphable_transformation(obj.TransformMorphableModel,obj.CorrespondingShape,obj.InlierWeights');
                       if obj.CorrespondencesUseOrientation
                          floatingFeatures = single([update.Vertices, update.VertexNormals]);
                       else
                          floatingFeatures = single([update.Vertices, zeros(obj.FloatingShape.nVertices,3)]);
                       end
               end
               obj.Iteration = obj.Iteration+1;
               obj.FloatingShape.Vertices = double(floatingFeatures(:,1:3));
               if obj.Display
                  waitbar(obj.Iteration/obj.NumIterations,h);
                  %pause;
               else
                  if obj.FloatingShape.Visible, drawnow expose;end
               end
           end        
           if obj.Display
              delete(h); 
           end
           if obj.UpSampleTarget, obj.TargetShape = TargetBck; end
           obj.FloatingShape.UserData.InlierWeights = double(obj.InlierWeights);
       end
   end
end
