classdef shape3D < superHandleClass
    % Object-oriented Handle class to manipulate and visualize 3D shapes,
    % given as point-lists or meshes
    % Key feature(s), include automatic interaction with matlab patch object
    % for real time visualisation; import and export with .obj wavefront
    % files, link with matlab triangulation object for mesh queries, 
    % interaction with MeshMonk 3D shape mapping software
    % and a full support into the viewer3D class which is a matlab
    % figure class tuned for 3D visualizations.
    %
    % Written by Peter Claes, peter.claes@kuleuven.be, copyright 2017
    properties % Main properties
        Vertices = [];% List of Nx3 Vertices
        Faces = [];% List of Nx3 Faces, accepting triangles only
    end
    properties (Dependent = true) % Main dependent properties
        VertexNormals;% Normal in each Vertex
        FaceNormals;% Normal in each Face
        Adjacency;% Adjacenct matrix coding the edge connectivities between vertices
        CentroidSize;% Size of the form
        nVertices;% Number of Vertices
        nFaces;% Number of Faces
    end
    properties
        FlipNormals = false;% Some meshes have a different rotation defined in the triangles, use this to correct normal orientation
        VertexValue = [];% A value per vertex, can be displayed in ColorMode = 'Indexed';
        VertexRGB = [];% A Red Green Blue color per vertex, can be given or extracted from the texturemap if present.
    end
    properties % COLOR options
        SingleColor; % Static color for vertices and edges
        TextureMap; % 2D image
        UV;% 3D -> 2D texture mapping, Nx2 coordinates 
    end
    properties % DisplayModes
        ViewMode = 'solid'; % solid, wireframe, points, solidwireframe, dottedwireframe
        ColorMode = 'single';% single, indexed, texture
        ReflectionMode = 'gouraud';% flat, phong, gouraud from low to high reflection quality
        Material = 'default';% default, facial, metal, shiny, image3D, dull
        Alpha = 1;% transparancy factor
        VertexSize = 5;% size of vertices when rendering
    end
    properties (Hidden = true, Transient = true)
        PatchHandle;% Matlab based Rendering Object;
    end
    properties (Dependent = true, Hidden = true)
        Centroid;% Point of gravity, average of all vertices
        TriHandle;% Matlab based Triangulation Object, used to obtain particular interactions;
        Level;% 0 = empty vertices and Faces, 1 empty faces (point list), 2 mesh
        TextureMapChannels;% number of channels in texture map, typically, 3 (RGB) or 1 (Gray Image)
        TextureMapSize;% first two dimensions of the texture map
        MeshMonkVertices;% A data convertor towards MeshMonk mapping routines
        MeshMonkFaces;% A data convertor towards MeshMonk mapping routines
        MeshMonkFeatures;% A data convertor towards MeshMonk mapping routines
    end
    properties (Transient = true)
        Selected = false; % When selected, viewers can change viewing properties, otherwise not
        Visible = false;% When visible, a patch handle is created, and visualizes any changes
        RenderAxes = [];% A patch handle requires an axes object for creation
    end
    properties % closing properties
        UserData;% Property to assign data by user
        Tag = 'Name Tag';% Tag is much like a name
    end
    properties (Dependent = true)
        Content;% Empty object, Point Cloud, Mesh 
    end
    methods %Constructor
        function obj = shape3D(varargin)
            obj = obj@superHandleClass(varargin{:});
            obj.SingleColor = 0.2 + (1-0.2).*rand(1,3);% Initialize a random color
        end
    end
    methods % GETTING
        function out = get.TriHandle(obj)
            if obj.Level<2, out = []; return; end
            out = triangulation(obj.Faces,obj.Vertices);% triangulation is a matlab construct that helps in a lot of computations
        end
        function out = get.Level(obj)
           if isempty(obj.Vertices), out = 0; return; end % Detecting that fields are empty
           if isempty(obj.Faces), out = 1; return; end% Detecting that we have a point list only
           out = 2;% Full mesh information, vertices and faces, is available
        end
        function out = get.Content(obj)
            switch obj.Level
                case 0
                    out = 'Empty';
                case 1
                    out = 'Point List';
                case 2
                    out = 'Mesh';
            end
        end
        function out = get.nVertices(obj)
            out = size(obj.Vertices,1);
        end
        function out = get.nFaces(obj)
            out = size(obj.Faces,1);
        end
        function out = get.FaceNormals(obj)
            if obj.Level<2, out = []; return; end
            out = faceNormal(obj.TriHandle);
            if obj.FlipNormals, out = -1*out;end
        end
        function out = get.VertexNormals(obj)
            if obj.Level<2, out = []; return; end
            out = vertexNormal(obj.TriHandle);
            if obj.FlipNormals, out = -1*out;end
        end
        function out = get.VertexRGB(obj)
            if isempty(obj.TextureMap)||isempty(obj.UV), out = double(obj.VertexRGB)/255;return;end
            out = textureMapValues(obj);
        end
        function out = get.Adjacency(obj)
            f = double(obj.Faces);
            out = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
                       [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
                                                                   1.0);
            % avoid double links
            out = double(out>0);
        end
        function out = get.Faces(obj)
            out = double(obj.Faces);
        end
        function out = get.RenderAxes(obj)
            if isempty(obj.RenderAxes)||~ishandle(obj.RenderAxes), out = []; return; end 
            out = obj.RenderAxes;
        end
        function out = get.PatchHandle(obj)
            if isempty(obj.PatchHandle)||~ishandle(obj.PatchHandle), out = []; return; end
            out = obj.PatchHandle;
        end
        function out = get.ViewMode(obj)
           if obj.Level<2, out = 'points'; return; end 
           out = obj.ViewMode;
        end
        function out = get.Centroid(obj)
            out = mean(obj.Vertices,1); 
        end
        function out = get.CentroidSize(obj)
            if obj.nVertices==0, out = 0; return; end
            Differences = repmat(obj.Centroid,obj.nVertices,1)-obj.Vertices;
            Distances = sqrt(sum(Differences.^2,2));
            out = sqrt(sum(Distances.^2)/obj.nVertices);
        end
        function out = get.Visible(obj)
            if isempty(obj.PatchHandle), out = false; return; end
            out = obj.Visible;
        end
        function out = get.TextureMapSize(obj)
            if isempty(obj.TextureMap), out = [0 0]; return; end
            out = [size(obj.TextureMap,2) size(obj.TextureMap,1)];
        end
        function out = get.TextureMapChannels(obj)
            if isempty(obj.TextureMap), out = [0 0]; return; end
            out = size(obj.TextureMap,3);
        end
        function out = get.MeshMonkVertices(obj)
                 out = single(obj.Vertices);% MeshMonk expects single precision data
        end
        function out = get.MeshMonkFaces(obj)
                 out = obj.Faces-1;% MeshMonk is C++ code, iterating from 0 instead of 1
                 out = uint32(out);% MeshMonk required unsigned integer data;
        end
        function out = get.MeshMonkFeatures(obj)
                 out = [obj.Vertices, obj.Gradient];% MeshMonk features are Vertices + Gradients, 6 dimensional
        end
    end
    methods % SETTING
        function obj = set.Vertices(obj,in)
            [a,b] = size(in);
            if a==3&&b~=a, in = in';end
            obj.Vertices = in;
            checkVerticesColorUpdate(obj,'update vertices');
        end
        function obj = set.Faces(obj,in)
            [a,b] = size(in);
            if a==3&&b~=a, in = in';end
            obj.Faces = superHandleClass.convertUInt(in);
            setPatch(obj,'Faces');setPatch(obj,'ViewMode');
        end
        function obj = set.VertexRGB(obj,in)
            [a,b] = size(in);
            if a==3&&b~=a, in = in';end
            %warning('off','All');
            if max(in(:))<=3
               in = uint8(in*255);
            else
               in = uint8(in);
            end
            obj.VertexRGB = in;
            %warning('on','All');
            checkVerticesColorUpdate(obj,'update rgbs');
        end
        function obj = set.VertexValue(obj,in)
            if size(in,1)==1, in = in'; end
            obj.VertexValue = double(in);
            checkVerticesColorUpdate(obj,'update values');
        end
        function obj = set.Visible(obj,in)
            if ~islogical(in), error('input must be logical true or false'); end
            switch in
                case true
                    obj.Visible = true;
                    createPatch(obj);setPatch(obj,'All');
                case false
                    obj.Visible = false;
                    deletePatch(obj);
                otherwise
                    error('Visible must be logical true or false');
            end
        end
        function obj = set.ViewMode(obj,in)
              if ~(strcmpi(in,'Solid') ||... 
                      strcmpi(in,'Wireframe') ||... 
                      strcmpi(in,'Points') ||...
                      strcmpi(in,'DottedWireframe') ||...
                      strcmpi(in,'SolidWireframe'))
                      error('ViewMode must be Solid, Wireframe, SolidWireframe, DottedWireframe  or Points');
              end
              obj.ViewMode = in;
              setPatch(obj,'ViewMode');
        end
        function obj = set.VertexSize(obj,in)
            obj.VertexSize = in;
            setPatch(obj,'markersize');
        end
        function obj = set.Alpha(obj,in)
            if (in<0) || (in>1)
               error('Alpha Value must be in the range of [0 1]');
            end
            obj.Alpha = in;
            setPatch(obj,'Alpha');
        end
        function obj = set.ColorMode(obj,in)
            switch lower(in)
                case 'texture'
                   if ~(size(obj.VertexRGB,1)==obj.nVertices), return; end %#ok<*MCSUP>
                case 'indexed'
                   if ~(length(obj.VertexValue)==obj.nVertices), return; end
                case 'single'
                otherwise
                   error('ColorMode must be Texture, Indexed (Values) or Single');   
            end
            obj.ColorMode = in;
            setPatch(obj,'ColorMode');
        end
        function obj = set.ReflectionMode(obj,in)
                 if ~(strcmpi(in,'none') ||... 
                      strcmpi(in,'flat') ||... 
                      strcmpi(in,'gouraud') ||... 
                      strcmpi(in,'phong'))
                      error('Lighting must be none, flat, gouraud or phong');
                 end
                 obj.ReflectionMode = in;
                 setPatch(obj,'ReflectionMode');
        end
        function obj = set.Material(obj,in)
                 if ~(strcmpi(in,'Facial') ||... 
                      strcmpi(in,'Dull') ||... 
                      strcmpi(in,'Shiny') ||... 
                      strcmpi(in,'Metal') ||...
                      strcmpi(in,'Image3D') ||...
                      strcmpi(in,'Default'))
                      error('Material must be Facial, Dull, Shiny, Metal or Default');
                 end
                 obj.Material = in;
                 setPatch(obj,'Material');
        end
        function obj = set.Selected(obj,in)
            if ~islogical(in), error('input must be logical true or false'); end
            obj.Selected = in;
        end
        function obj = set.RenderAxes(obj,in)
            if isempty(in), obj.RenderAxes = []; return; end
            if ~ishandle(in), obj.RenderAxes = []; return; end
            if ~strcmp(get(in,'Type'),'axes'),obj.RenderAxes = []; return; end
            obj.RenderAxes = in;
            setPatch(obj,'Axes');
        end
        function obj = set.UV(obj,in)
            [a,b] = size(in);
            if a==2&&b~=a, in = in';end
            obj.UV = in;
            checkVerticesColorUpdate(obj,'update rgbs');
        end
        function obj = set.SingleColor(obj,Color)
                 if ~(numel(Color)==3), error('Invalid single Color argument'); end
                 if ~(size(Color,2)==3), Color = Color'; end
                 obj.SingleColor = Color;
                 setPatch(obj,'ColorMode');
        end
        function obj = set.TextureMap(obj,in)
            if max(in(:))<=1, in= in*255;end% convert [0 1] to [0 255];
            in = uint8(in);% convert to unsigned integer class;
            obj.TextureMap = in;
            checkVerticesColorUpdate(obj,'update rgbs');
        end
    end
    methods % SIMPLE INTERFACING
        function out = centerVertices(obj)
            if nargout==1, out = clone(obj); else out = obj; end %#ok<*SEPEX>
            center = out.Centroid;
            out.Vertices = out.Vertices - repmat(center,out.nVertices,1);
        end
        function out = textureMapValues(obj,uv)
            if nargin<2, uv = obj.UV; end
            [X,Y] = textureMapGrid(obj);
            out = zeros(size(uv,1),obj.TextureMapChannels);
            for i=1:1:obj.TextureMapChannels
                out(:,i) = interp2(X,Y,double(obj.TextureMap(:,:,i)),uv(:,1),uv(:,2),'linear');
            end
            if obj.TextureMapChannels==1, out = repmat(out,1,3);end% Gray values
            out = double(out)/255;% convert to double precision, range [0 1], needed for matlab rendering 
        end
        function [X,Y,I,J] = textureMapGrid(obj)
            % function to retrieve the grid underlying the image of the
            % texture map, X and Y are expressed in UV dimensions, I and J
            % are indices, for easy reference into 3D matrices like the
            % textureMap
            n = obj.TextureMapSize(1);
            m = obj.TextureMapSize(2);
            if max(obj.UV(:))<=1 && min(obj.UV(:))>=0
               x = [0 1];
               y = [0 1];
            else
               x = [1 n];
               y = [1 m];
            end
            xmin = min(x(:)); ymin = min(y(:));
            xmax = max(x(:)); ymax = max(y(:));
            dx = (xmax-xmin)/max(n-1,1);
            dy = (ymax-ymin)/max(m-1,1);
            xx = xmin:dx:xmax;yy = ymin:dy:ymax;
            [X,Y] = meshgrid(xx,yy);
            ii=1:length(xx);jj=1:length(yy);
            [I,J] = meshgrid(ii,jj);
        end
        function importObject(obj,filename,varargin)
            % changing directory
                 Input = find(strcmp(varargin,'Path'));
                 if ~isempty(Input), cd(varargin{Input+1});end
             % loading    
               in = load(filename);
               loaded = fields(in);
               copy(in.(loaded{1}),obj);
               delete(in.(loaded{1}));
        end
        function out = Findex2Vindex(obj,in)
           if obj.Level<2, out = []; return; end
           % converts a Face Index (selection of mesh faces) into Vertex index
           % (selection of mesh vertices)
           faces = obj.Faces(in,:);
           out = unique(faces(:));
        end
        function out = Vindex2Findex(obj,in)
           if obj.Level<2, out = []; return; end
           % converts a Vertex Index (selection of mesh vertices) into Face index
           % (selection of mesh faces)
           tmp = ismember(obj.Faces,in);
           [j,~] = find(tmp==0);
           out = setdiff(1:obj.nFaces,j);
        end
        function [Vindex,Findex] = getVindexFindex(obj,varargin)
           % tests the input arguments of a function in varargin whether Face index or
           % Vertex index are given
           % using a given Vertex index
           Input = find(strcmp(varargin, 'VertexIndex'));
           if isempty(Input)
               Vindex = 1:obj.nVertices;
           else
               Vindex = varargin{Input+1};
           end
           if obj.Level<2, Findex = []; return; end
           Findex = Vindex2Findex(obj,Vindex);
           % using a given face index (OVERIDES A GIVEN VERTEX INDEX!!!!!)
           Input = find(strcmp(varargin, 'FaceIndex'));
           if ~isempty(Input)
               Findex = varargin{Input+1};
               Vindex = Findex2Vindex(obj,Findex);
           end
        end
        function out = intraDistances(obj,varargin)
            if obj.Level<2, out = []; return; end% A Mesh is required
            [Vindex,~] = getVindexFindex(obj,varargin{:});
            Input = find(strcmp(varargin, 'Type'));
            if ~isempty(Input)
               distance_type = varargin{Input+1};
            else
               distance_type = 'triangle';
            end
            Input = find(strcmp(varargin, 'Max'));
            if ~isempty(Input)
               max_dist = varargin{Input+1};
            else
               max_dist = inf;
            end
            switch distance_type
                case 'triangle'
                    A = obj.Adjacency;
                case 'edge'
                    A = vertexAdjacency(obj);
                otherwise
                    return;
            end
            out = nan*zeros(1,obj.nVertices);
            To_Do = 1:obj.nVertices;
            To_Do = setdiff(To_Do,Vindex);
            current = Vindex;
            out(current) = 0;  
            Neighbors = 1;
            counter = 0; 
            while ~isempty(Neighbors)&&(counter<obj.nVertices)
                 counter = counter + 1;
                 [current_index,Neighbors,dist] = find(A(current,:));
                 dist = dist(:)' + out(current(current_index));
                 [~,I,~] = unique(Neighbors,'first');
                 dist_first = dist(I);
                 [Neighbors,I,~] = unique(Neighbors,'last');
                 dist_last = dist(I);
                 dist = min([dist_first;dist_last]);
                 [Neighbors,~,J] = intersect(To_Do,Neighbors);
                 dist = dist(J);
                 out(Neighbors) = dist;
                 if isempty(find(out(Neighbors)<=max_dist,1)), break; end
                 current = Neighbors;
                 To_Do = setdiff(To_Do,Neighbors);
            end
            out(find(out>max_dist,1)) = nan;
        end
        function out = vertexAdjacency(obj)
           % generates Adjancency matrix with distances between connected points
            out = obj.Adjacency;
            [i,j] = find(out);
            ind = sub2ind(size(out),i,j);
            distances = sqrt(sum((obj.Vertices(i,:)-obj.Vertices(j,:)).^2,2));
            out(ind) = distances;
            % make sure that all edges are symmetric
            out = (out+out')/2;
        end
        function delete(obj)
            if ~isempty(obj.PatchHandle), delete(obj.PatchHandle); end
            delete@superHandleClass(obj);
        end
        function out = crop(obj,varargin)
           %reduces 'shape' so that it only contains the vertices indicated in 'indices'
           %	: indices : vector containing the indices of the vertices you want to
           %	keep or delete dependent on the action in varargin
           if obj.Level<1, out = []; return; end% minimum level is landmarks
           % example usage: newshape = crop(shape,'VertexIndex',index);
           if nargout == 1, obj = clone(obj);out = obj;end% create a copy if output is requested
           Input = find(strcmp(varargin, 'Action'));
           if isempty(Input), action = 'crop'; else, action = varargin{Input+1};end
           [Vindex,Findex] = getVindexFindex(obj,varargin{:});
           if strcmpi(action,'delete')% turn into crop data
              if isempty(Vindex), return; end
              fullindex = 1:obj.nVertices;
              Vindex = setdiff(fullindex,Vindex);
              Findex = Vindex2Findex(obj,Vindex);
           end
           %throw away all quads containing an index not in the indices we want to keep
           if ~isempty(obj.Faces),[~,LOC] = ismember(obj.Faces,Vindex); obj.Faces = LOC(Findex,:);end
           obj.Vertices = obj.Vertices(Vindex,:);
           if~isempty(obj.VertexValue), obj.VertexValue = obj.VertexValue(Vindex);end
           if~isempty(obj.VertexRGB), obj.VertexRGB = obj.VertexRGB(Vindex,:);end
           if~isempty(obj.UV), obj.UV = obj.UV(Vindex,:);end
        end
        function out = getBoundary(obj)
           if obj.Level<2, out = []; return; end% requires a mesh
           out = unique(freeBoundary(obj.TriHandle));
        end
        function out = cleanManifold(obj)
           if nargout == 1, obj = clone(obj);out = obj; end 
           % remove nan vertices
           [i,~] = find(isnan(obj.Vertices));
           if ~isempty(i), crop(obj,'VertexIndex',unique(),'Action','delete'); end
           % remove vertices not part in triangles
           crop(obj,'VertexIndex',unique(obj.Faces(:)));
        end
        function out = getTriangleMetrics(obj)
          % Generates a series of quality metrics for the triangles in a mesh
          % On can subsequently opt ot remove or downscale their importance
          if obj.Level<2, out = []; return; end% requires a mesh
          faces = obj.Faces';
          vertices = obj.Vertices';
          nrF = size(faces,2);
          LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
          for i=1:1:3
              LOC(:,:,i) = reshape(vertices(i,faces(:)),3,size(faces,2));
              AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
              AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
          end
          areas = 0.5*sqrt(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2);
          Z_areas = (areas-mean(areas))/std(areas);
          out.TriangleAreas = areas;
          out.AreaZscores = Z_areas;
          out.TriangleQuality = pdetriq(vertices,faces);
        end
        function out = numberTrianglesPerVertex(obj)
          out = hist(obj.Faces(:),1:obj.nVertices);% nr of triangles a vertex participates in    
        end
        function out = numberEdgesPerVertex(obj)
          out = full(sum(obj.Adjacency)); % nr of adjecant elements per vertex  
        end
        function [outTriangles,outVertices] = getBadlySizedTriangles(obj,zscore)
            % extracts bad triangles, bad is defined in the context of
            % relative triangle size, given all triangles in the mesh
            if obj.Level<2, outTriangles = [];outVertices = []; return; end% requires a mesh
            out = getTriangleMetrics(obj);
            outTriangles = find(abs(out.AreaZscores)>zscore);
            if isempty(outTriangles), outVertices = []; return; end
            outVertices = Findex2Vindex(obj,outTriangles);
        end
        function out = edgeLengths(obj)
               DA = vertexAdjacency(obj);
               DA = DA(:);
               out = full(DA(DA>0));
        end
        function out = flipFaces(obj)
            % function to rearrange vertex order in faces, which influences
            % the orientation of normal computations
            if nargout==1, out = clone(obj); else out = obj; end
            tmp = obj.Faces;
            tmp(:,2) = obj.Faces(:,3);
            tmp(:,3) = obj.Faces(:,2);
            obj.Faces = tmp;
        end
        function [out,degrees] = getNormalOrientation(obj)
            warning off;
            obj = cleanManifold(obj);
            Vec2Centroid = repmat(obj.Centroid,obj.nVertices,1)-obj.Vertices;
            Normals = obj.VertexNormals; 
            degrees = zeros(1,obj.nVertices);               
            parfor i=1:1:obj.nVertices
                [~,degrees(i)] = shape3D.vectorangle(Vec2Centroid(i,:)',Normals(i,:)');  
            end
            if nanmean(degrees)>90
               out = 1;% outward pointing normals
            else
               out = 0;% inward pointing normals
            end
            warning on;
        end
        function out = getResolution(obj)
            out = nanmedian(edgeLengths(obj));
        end
    end
    methods % RENDERING & VISUALIZATION
        function out = viewer(obj,in)
            if nargin<2, in = viewer3D;end
            out = in;
            obj.RenderAxes = out.RenderAxes;
            obj.Visible = true;
            obj.Selected = true;
        end
        function out = figure(obj,varargin)
            if isempty(varargin)
               out = figure;axis equal;
            elseif ~strcmp(get(varargin{1},'Type'),'figure')
               out = figure;axis equal; 
            else
               out = varargin{1};
               figure(out);axis equal;
            end
            obj.RenderAxes = gca;
            obj.Visible = true;
        end
        function createPatch(obj)
            if isempty(obj.RenderAxes), obj.RenderAxes = gca; end % a patch object needs an axis system
            deletePatch(obj); % clear any previous patch object
            obj.PatchHandle = patch('Parent',obj.RenderAxes,'Visible','off','UserData',obj);
        end
        function deletePatch(obj)
            if ~isempty(obj.PatchHandle), delete(obj.PatchHandle); end
        end
        function setPatch(obj,action)
            if isempty(obj.PatchHandle), return; end% There is no active patch handle
            if obj.Level<1, return; end% There is nothing to display
            switch lower(action)
                case 'all'
                    setPatch(obj,'Vertices');setPatch(obj,'ColorMode');setPatch(obj,'ViewMode');
                    setPatch(obj,'Faces');setPatch(obj,'Viewer');setPatch(obj,'Alpha');
                    setPatch(obj,'Material');setPatch(obj,'ReflectionMode');setPatch(obj,'MarkerSize');
                    setPatch(obj,'Visible');
                    return;
                case 'shape'
                    setPatch(obj,'Vertices');setPatch(obj,'Faces');
                    return;
                case 'verticescolor'
                    setPatch(obj,'Vertices');setPatch(obj,'Faces');setPatch(obj,'ColorMode');
                    return;
                case 'axes'
                    prop = {'Parent'};val{1} = obj.RenderAxes;
                case 'vertices'
                    prop = {'Vertices'}; val{1} = obj.Vertices;
                case 'faces'
                    if obj.Level==1 % point list
                       val{1} = [(1:1:size(obj.Vertices,1));nan*ones(1,size(obj.Vertices,1))];
                    else % mesh
                       val{1} = obj.Faces;
                    end
                    prop = {'Faces'};
                case 'alpha'
                    alpha(obj.PatchHandle,obj.Alpha);
                    return;
                case 'material'
                    prop={'AmbientStrength';'DiffuseStrength';'SpecularStrength';...
                          'SpecularExponent';'SpecularColorReflectance'};
                      switch lower(obj.Material)
                          case 'facial'
                              val = {0.5; 0.4; 0.1; 'default'; 'default'};
                          case 'dull'
                              val = {0.3; 0.8; 0.0; 10; 1.0};
                          case 'metal'
                              val = {0.3; 0.3; 1.0; 25; 0.5};
                          case 'shiny'
                              val = {0.3; 0.6; 0.9; 20; 1.0};
                          case 'image3d'
                              val = {0.5; 0; 0; 10000000000; 0};
                          case 'default'
                              val = {'default'; 'default'; 'default'; 'default'; 'default'};
                          otherwise
                              return
                      end
                case 'colormode'
                    prop = {'FaceVertexCData'};
                    switch lower(obj.ColorMode)
                        case 'texture'
                           val = obj.VertexRGB;
                        case 'indexed'
                           val = obj.VertexValue;
                        case 'single'
                             val = repmat(obj.SingleColor,size(get(obj.PatchHandle,'Vertices'),1),1);
                        otherwise
                             return;
                    end
                case 'viewmode'
                    prop = {'EdgeColor'; 'FaceColor'; 'LineStyle';'Marker';'MarkerEdgeColor';'MarkerFaceColor'};
                    switch lower(obj.ViewMode)
                         case 'solid'
                             val{1,1} = 'none'; val{2,1} = 'interp'; val{3,1} = 'none';
                             val{4,1} = 'none'; val{5,1} = 'auto';   val{6,1} = 'none';
                         case 'wireframe'
                             val{1,1} = 'interp'; val{2,1} = 'none'; val{3,1} = '-';
                             val{4,1} = 'none';   val{5,1} = 'auto'; val{6,1} = 'none';
                         case 'dottedwireframe'
                             val{1,1} = 'interp'; val{2,1} = 'none'; val{3,1} = ':';
                             val{4,1} = 'none';   val{5,1} = 'auto'; val{6,1} = 'none';
                         case 'solidwireframe'
                             val{1,1} = [0 0 0.35]; val{2,1} = 'interp'; val{3,1} = '-';
                             val{4,1} = 'none';     val{5,1} = 'auto';   val{6,1} = 'none';
                         case 'points'
                             val{1,1} = 'none'; val{2,1} = 'none'; val{3,1} = '-';
                             val{4,1} = '.';    val{5,1} = 'flat'; val{6,1} = 'flat';
                         otherwise
                             return;
                    end        
                case 'reflectionmode'
                    prop = {'FaceLighting'};
                    val{1} = obj.ReflectionMode;
                case 'visible'
                    prop = {'Visible'};
                    if obj.Visible, val = {'on'}; else val = {'off'}; end
                case 'markersize'
                    prop = 'MarkerSize';
                    val{1} = obj.VertexSize;
                otherwise
                    return;
            end
            pvPairs = [prop,val]';
            set(obj.PatchHandle,pvPairs{:});
        end
        function out = checkVerticesColorUpdate(obj,action)
            out = false;
            switch lower(action)
                case 'update rgbs'
                    if strcmpi(obj.ColorMode,'texture')&& (size(obj.VertexRGB,1)==obj.nVertices), setPatch(obj,'VerticesColor'); out = true; end
                case 'update values'
                    if strcmpi(obj.ColorMode,'indexed')&& (length(obj.VertexValue)==obj.nVertices), setPatch(obj,'VerticesColor'); out = true;end
                case 'update vertices'
                    switch lower(obj.ColorMode)
                           case 'single'
                                setPatch(obj,'VerticesColor'); out = true;          
                           case 'indexed'
                                if (obj.nVertices==length(obj.VertexValue)), setPatch(obj,'VerticesColor'); out = true; end
                           case 'texture'
                                if (obj.nVertices==size(obj.VertexRGB,1)), setPatch(obj,'VerticesColor'); out = true; end
                    end
                otherwise
            end
        end
    end
    methods (Static = true) % STATIC functions
        function varargout = importGetFile
             [filename, pathname,filterindex] = uigetfile({'*.mat','Matlab Object';...
                                                           '*.obj','Obj Wavefront';...
                                                           '*.obj','Obj MevisLab';...
                                                            },'MultiSelect','on');
             if isequal([filename,pathname],[0,0]),varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4} = [];return; end
             if ~iscell(filename)
                varargout{1} = 1;
                varargout{2} = {filename};
             else
                varargout{1} = size(filename,2);
                varargout{2} = filename;
             end
             varargout{3} = pathname;
             switch filterindex
                 case 1
                     varargout{4} = 'Matlab Object';
                 case 2
                     varargout{4} = 'Wavefront';
                 case 3
                     varargout{4} = 'WavefrontMevisLab';
                 otherwise
             end
        end
        function obj = import(filename,path,mevis)
                 if nargin == 0
                    [nrfiles,filename,path,filterindex] = shape3D.importGetFile;
                    if isempty(filename), obj = shape3D; return; end% create an empty object
                    obj = cell(1,nrfiles);
                    for i=1:1:nrfiles
                        obj{i} = shape3D;
                        obj{i}.Tag = filename{i}(1:end-4);
                        switch filterindex
                            case 'Matlab Object'
                                importObject(obj{i},filename{i},path);
                            case 'Wavefront'
                                importWavefront(obj{i},filename{i},path,false);
                            case 'WavefrontMevisLab'
                                importWavefront(obj{i},filename{i},path,true);    
                            otherwise
                        end
                    end
                    if nrfiles == 1, obj = obj{1}; end% no need to return a cell array
                 else
                    if nargin<3, mevis = false; end
                    obj = shape3D;
                    obj.Tag = filename(1:end-4);
                    switch filename(end-3:end)
                        case '.obj'                        
                            importWavefront(obj,filename,path,mevis);
                        case '.mat'
                            importObject(obj,filename,path);
                        otherwise
                            error('Unknown file format');
                    end
                 end
        end
        function [out,outd] = vectorangle(v1,v2)
              T = v1'*v2;
              N = sqrt((v1'*v1)*(v2'*v2));
              out = T/N;outd = acosd(out);
        end
    end
end