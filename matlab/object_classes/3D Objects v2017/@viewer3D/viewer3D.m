classdef viewer3D < superHandleClass
    % this class is a class of 3D viewers
    % It builds on a matlab figure, axes and light handle, and interacts
    % with handles from the shape3D class. 
    properties
        % figure object
        Figure = []; % Matlab figure handle
        % Mode status
        Mode = 'camera'% the current mode of the viewer (none, camera, light)
        SelectionMode = 'none'% the selection mode (none, landmarks, brush, fill, full)
        LandmarkSelection = [];
        AreaSelection = [];
        % Light objects
        SceneLight = [];% Matlab light handle
        SceneLightRotMode = 'vertical';
        SceneLightLinked = false;
        % Axes object
        RenderAxes = [];% Matlab axes handle
        % Toolbar objects
        Toolbar = [];
        % Colorbar objects
        Colorbar = [];
        % Context menu objects
        ContextMenu = [];
        % General
        MotionData = [];
        Action = 'none';
        SelectionActive = false;
        Selected = true;
        Status = 'Ready';
        Tag = '3D Viewer'; 
        UserData= [];% User specific data;
    end
    properties (Dependent = true)% Dependent on the Figure, Axes or Light object(s)
        BackgroundColor;
        Renderer;
        AxesVisible;
        AxesGrid;
        AxesBox;
        AxesWallColor;
        AxesXColor;
        AxesYColor;
        AxesZColor;
        CameraPosition;
        CameraTarget;
        CameraUpVector;
        CameraViewAngle;
        DataAspectRatio;
        SceneLightMode;
        SceneLightPosition;
        SceneLightColor;
        SceneLightVisible;
        Visible;
        CamProjection;
        ColorbarVisible;
    end
    properties (Hidden = true)
        OldBackgroundColor = [];
        SelectionSphereRadius = 10;
        SelectionSphereCenter = [0 0 0];
        SelectionTH = 200;
        SelectionDistances = [];
        AreaSelectionDummy = [];
        AreaSelectionIndex = [];
        AreaSelectionDummyIndex = [];
        MuteCallbacks = false;
    end
    properties (Hidden = true, Dependent = true)
        ShapeChildren;
        CurrentShape;
        NoCurrentShapeMsg;
    end
    methods %Constructor
        function obj = viewer3D(varargin)
          obj = obj@superHandleClass(varargin{:});
          constructViewer(obj);
          updateToolbar(obj);
          updateContextMenu(obj);
          setMousePointer(obj);
        end % 3D viewer Constructor
    end 
    methods % SETTING & GETTING FIGURE HANDLE
        function out = get.Figure(obj)
            out = obj.Figure;
            if ~ishandle(out), out = []; end
        end
        function out = get.BackgroundColor(obj)
            if isempty(obj.Figure), out = []; return; end
            out = get(obj.Figure,'Color');
        end
        function obj = set.BackgroundColor(obj,in)
            if isempty(obj.Figure), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.Figure,'Color',in);
        end
        function out = get.Renderer(obj)
            if isempty(obj.Figure), out = []; return; end
            out = get(obj.Figure,'renderer');
        end
        function obj = set.Renderer(obj,in)
            if isempty(obj.Figure), return; end
            if ~(strcmpi(in,'opengl') ||strcmpi(in,'zbuffer'))
                error('renderer must be opengl or zbuffer');
            end
            set(obj.Figure,'renderer',in);
            updateContextMenu(obj);
        end        
        function out = get.Visible(obj)
            if isempty(obj.Figure), out = []; return; end
            switch get(obj.Figure,'Visible')
                case 'on'
                    out = true; 
                case 'off'
                    out = false;
            end
        end
        function obj = set.Visible(obj,in)
            if isempty(obj.Figure), return; end
            switch in
                case true
                    set(obj.Figure,'Visible','on');
                case false
                    set(obj.Figure,'Visible','off');
                otherwise
            end
        end
    end
    methods % SETTING & GETTING AXES HANDLE
        function out = get.RenderAxes(obj)
            out = obj.RenderAxes;
            if ~ishandle(out), out = []; end
        end
        function out = get.AxesVisible(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            switch get(obj.RenderAxes,'Visible')
                case 'on'
                    out = true;
                case 'off'
                    out = false;
            end
        end
        function obj = set.AxesVisible(obj,in)
            if isempty(obj.RenderAxes), return; end
            switch in
                case true
                    set(obj.RenderAxes,'Visible','on');
                case false
                    set(obj.RenderAxes,'Visible','off');
                otherwise
                    return;
            end
            updateContextMenu(obj)
        end
        function out = get.AxesGrid(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            switch get(obj.RenderAxes,'XGrid')
                case 'on'
                    out = true;
                case 'off'
                    out = false;
            end
        end
        function obj = set.AxesGrid(obj,in)
            if isempty(obj.RenderAxes), return; end
            switch in
                case true
                    grid(obj.RenderAxes,'on');
                case false
                    grid(obj.RenderAxes,'off');
                otherwise
                    return;
            end
            updateContextMenu(obj);
        end
        function out = get.AxesBox(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            switch get(obj.RenderAxes,'Box')
                case 'on'
                    out = true;
                case 'off'
                    out = false;
            end
        end
        function obj = set.AxesBox(obj,in)
            if isempty(obj.RenderAxes), return; end
            switch in
                case true
                    set(obj.RenderAxes,'Box','on');
                case false
                    set(obj.RenderAxes,'Box','off');
                otherwise
                    return;
            end
            updateContextMenu(obj);
        end
        function out = get.AxesWallColor(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'Color');
        end
        function obj = set.AxesWallColor(obj,in)
            if isempty(obj.RenderAxes), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.RenderAxes,'Color',in);
        end
        function out = get.AxesXColor(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'XColor');
        end
        function obj = set.AxesXColor(obj,in)
            if isempty(obj.RenderAxes), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.RenderAxes,'XColor',in);
        end
        function out = get.AxesYColor(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'YColor');
        end
        function obj = set.AxesYColor(obj,in)
            if isempty(obj.RenderAxes), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.RenderAxes,'YColor',in);
        end
        function out = get.AxesZColor(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'ZColor');
        end
        function obj = set.AxesZColor(obj,in)
            if isempty(obj.RenderAxes), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.RenderAxes,'ZColor',in);
        end
        function obj = set.CamProjection(obj,in)
            if isempty(obj.RenderAxes), return; end
            if ~(strcmp(in,'orthographic')||strcmp(in,'perspective'))
               error('Projection be must either orthographic or perspective');
            end
            set(obj.RenderAxes,'Projection',in);
            updateToolbar(obj);
            updateContextMenu(obj);
        end
        function out = get.CamProjection(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'Projection');
        end
        function out = get.CameraPosition(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'cameraposition');
        end
        function obj = set.CameraPosition(obj,in)
            if isempty(obj.RenderAxes), return; end
            set(obj.RenderAxes,'cameraposition',in);
        end
        function out = get.CameraTarget(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'cameratarget');
        end
        function obj = set.CameraTarget(obj,in)
            if isempty(obj.RenderAxes), return; end
            set(obj.RenderAxes,'cameratarget',in);
        end
        function out = get.CameraViewAngle(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'CameraViewAngle');
        end
        function obj = set.CameraViewAngle(obj,in)
            if isempty(obj.RenderAxes), return; end
            set(obj.RenderAxes,'CameraViewAngle',in);
        end
        function out = get.DataAspectRatio(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'dataaspectratio');
        end
        function obj = set.DataAspectRatio(obj,in)
            if isempty(obj.RenderAxes), return; end
            set(obj.RenderAxes,'DataAspectRatio',in);
        end
        function out = get.CameraUpVector(obj)
            if isempty(obj.RenderAxes), out = []; return; end
            out = get(obj.RenderAxes,'cameraupvector');
        end
        function obj = set.CameraUpVector(obj,in)
            if isempty(obj.RenderAxes), return; end
            set(obj.RenderAxes,'cameraupvector',in);
        end
    end
    methods % SETTING & GETTING LIGHT HANDLE
        function out = get.SceneLight(obj)
            out = obj.SceneLight;
            if ~ishandle(out), out = []; end
        end
        function out = get.SceneLightMode(obj)
            if isempty(obj.SceneLight), out = []; return; end
            out = get(obj.SceneLight,'Style');
        end
        function obj = set.SceneLightMode(obj,in)
            if isempty(obj.SceneLight), return; end
            if ~(strcmp(in,'infinite')||strcmp(in,'local'))
               error('Viewer scene light mode must be infinite or local');
            end
            set(obj.SceneLight,'Style',in);
        end
        function obj = set.SceneLightRotMode(obj,in)
            if isempty(obj.SceneLight), return; end
            if ~(strcmp(in,'vertical')||strcmp(in,'horizontal'))
               error('Viewer scene light rotation mode must be vertical or horizontal');
            end
            obj.SceneLightRotMode = in;
        end
        function out = get.SceneLightPosition(obj)
            if isempty(obj.SceneLight), out = []; return; end
            out = get(obj.SceneLight,'Position');
        end
        function obj = set.SceneLightPosition(obj,in)
            if isempty(obj.SceneLight), return; end
            set(obj.SceneLight,'Position',in);
        end
        function out = get.SceneLightColor(obj)
            if isempty(obj.SceneLight), out = []; return; end
            out = get(obj.SceneLight,'Color');
        end
        function obj = set.SceneLightColor(obj,in)
            if isempty(obj.SceneLight), return; end
            if ~numel(in)==3, return; end
            if size(in,1)==3, in = in'; end
            set(obj.SceneLight,'Color',in);
        end
        function out = get.SceneLightVisible(obj)
            if isempty(obj.SceneLight), out = []; return; end
            switch get(obj.SceneLight,'Visible')
                case 'on'
                    out = true;
                case 'off'
                    out = false;
            end
        end
        function obj = set.SceneLightVisible(obj,in)
            if isempty(obj.SceneLight), return; end
            switch in
                case true
                    set(obj.SceneLight,'Visible','on');
                case false
                    set(obj.SceneLight,'Visible','off');
                otherwise
                    return;
            end
            updateContextMenu(obj);
            updateToolbar(obj);
        end
        function obj = set.SceneLightLinked(obj,in)
            %if isempty(obj.SceneLight), return; end
            if ~islogical(in), return; end
            obj.SceneLightLinked = in;
            updateSceneLightPosition(obj);
            updateContextMenu(obj);
            updateToolbar(obj);
        end  
    end
    methods % SETTING & GETTING COLORBAR HANDLE
        function out = get.ColorbarVisible(obj)
            out = obj.Colorbar.Visible;
        end
        function obj = set.ColorbarVisible(obj,in)
            obj.Colorbar.Visible = in;
            updateToolbar(obj);
        end   
    end
    methods % GENERAL SETTING & GETTING
        function obj = set.Mode(obj,in)
            if ~(strcmp(in,'camera')||strcmp(in,'light'))
               error('Viewer mode must be, camera or light');
            end
            obj.Mode = in;
            updateContextMenu(obj);
            updateToolbar(obj);        
            setMousePointer(obj);
        end
        function obj = set.SelectionMode(obj,in)
            if ~(strcmp(in,'none')||strcmp(in,'landmark')||strcmp(in,'brush')||strcmp(in,'fill')||strcmp(in,'full'))
               error('Viewer selection mode must be none, landmark, brush, fill or full');
            end
            obj.SelectionActive = false;
            closeSelectionMode(obj);
            obj.SelectionMode = in;
            openSelectionMode(obj);
            updateToolbar(obj);
            updateContextMenu(obj);           
            setMousePointer(obj);
        end
        function obj = set.SelectionActive(obj,in)
                 obj.SelectionActive = in;
                 switch in
                     case true
                         obj.OldBackgroundColor = obj.BackgroundColor; %#ok<*MCSUP>
                         obj.BackgroundColor = [0.5 0.5 0.5];
                     case false
                         if ~isempty(obj.OldBackgroundColor)
                            obj.BackgroundColor = obj.OldBackgroundColor;
                            obj.OldBackgroundColor = [];
                         end
                 end
                 setMousePointer(obj);
        end          
        function obj = set.Tag(obj,in)
            if ~ischar(in), return; end
            obj.Tag = in;
            if ishandle(obj.Figure), set(obj.Figure,'Name',in); end
        end
        function out = get.ShapeChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            out = {};
            if isempty(patches), return; end
            counter = 0;
            nPatches = length(patches);
            out = cell(1,nPatches);
            good = zeros(1,nPatches);
            for p=1:1:nPatches
                UD = get(patches(p),'UserData');
                if ~isa(UD,'shape3D'),continue; end
                if strcmp(UD.Tag,'private'), continue; end
                if UD.Level<2, continue;end% you can only select on meshes
                counter = counter + 1;
                out{counter} = UD;
                good(counter) = 1;
            end
            if sum(good)==0, return;end
            out = out(find(good)); %#ok<*FNDSB>
        end
        function obj = set.ShapeChildren(obj,in)
            if isempty(in), return; end
            for c=1:1:length(in)
                if isa(in{c},'shape3D')
                   in{c}.RenderAxes = obj.RenderAxes;
                end
            end
        end
        function out = get.CurrentShape(obj)
                 out = [];
                 meshes = obj.ShapeChildren;
                 if isempty(meshes), out = []; return; end
                 for m=1:1:length(meshes)
                     if (meshes{m}.Selected) && (meshes{m}.Visible)
                         out = meshes{m};
                         break;% take the first in the list that is selected and visible
                     end
                 end            
        end
        function out = get.NoCurrentShapeMsg(obj) %#ok<*MANU>
                 out = 'No active shape! To activate a shape Select AND Display it. When Multiple scans are selected and displayed the first in the list of obj.ShapeChildren is made active, it is adviced to select and display only one scan whitin this viewer';
        end
        function obj = set.Status(obj,in)
                 obj.Status = in;
                 setMousePointer(obj);
        end
        function obj = set.Action(obj,in)
                 obj.Action = in;
                 setMousePointer(obj);
        end
        function out = get.AreaSelection(obj)
                 out = obj.AreaSelection;
                 if ~superHandleClass.isH(out), out = []; end
        end
        function out = get.AreaSelectionDummy(obj)
                 out = obj.AreaSelectionDummy;
                 if ~superHandleClass.isH(out), out = []; end
        end
        function out = get.LandmarkSelection(obj)
                 out = obj.LandmarkSelection;
                 if ~superHandleClass.isH(out), out = []; end
        end
        function obj = set.AreaSelectionIndex(obj,in)
            obj.AreaSelectionIndex = in;
            updateArea(obj,obj.AreaSelection,obj.AreaSelectionIndex);
            
        end
        function obj = set.AreaSelectionDummyIndex(obj,in)
            obj.AreaSelectionDummyIndex = in;
            updateArea(obj,obj.AreaSelectionDummy,obj.AreaSelectionDummyIndex);
        end
    end
    methods % CAMERA MANIPULATIONS
       function resetCamera(obj)
          camlookat(obj.RenderAxes);
          updateSceneLightPosition(obj);
       end
       function rotateCamera(obj,xy)
          xy = -xy;
          camorbit(obj.RenderAxes,xy(1),xy(2),'none')%unconstrained rotation
          updateSceneLightPosition(obj);
          drawnow nocallbacks;
       end
       function panCamera(obj,xy)
          xy = -xy;
          xy = xy*camva(obj.RenderAxes)/500;
          campan(obj.RenderAxes,xy(1),xy(2),'none')%unconstrained panning
          updateSceneLightPosition(obj);
          drawnow nocallbacks;
       end
       function zoomCamera(obj,xy,q)
          if isempty(q), q = max(-.9, min(.9, sum(xy)/70)); q = 1+q; end
          % hueristic avoids small view angles which will crash on solaris
          MIN_VIEW_ANGLE = .001;
          MAX_VIEW_ANGLE = 75;
          vaOld = camva(obj.RenderAxes);
          camzoom(obj.RenderAxes,q);
	      va = camva(obj.RenderAxes);
          %If the act of zooming puts us at an extreme, back the zoom out
          if ~((q>1 || va<MAX_VIEW_ANGLE) && (va>MIN_VIEW_ANGLE))
             set(obj.RenderAxes,'CameraViewAngle',vaOld);
          end
          drawnow nocallbacks;
       end
       function resetLight(obj)
          obj.SceneLightPosition = get(obj.RenderAxes,'CameraPosition');
       end
       function rotateLight(obj,xy)
          [az, el] = lightangle(obj.SceneLight);
          az = mod(abs(az),360);% Check if the light is on the other side of the object
          if az > 90 && az < 270
             xy(2) = -xy(2);
          end
          az = mod(az + xy(1), 360);
          el = mod(el + xy(2), 360);
          if abs(el) > 90
             el = 180 - el;
             az = 180 + az;
          end
          lightangle(obj.SceneLight, az, el);
          drawnow nocallbacks;
       end
       function screenPoint2currentPoint(obj)
          scrn_pt = get(0, 'PointerLocation');
          set(obj.Figure,'units','pixels')
          loc = get(obj.Figure, 'Position');
          % We need to compensate for an off-by-one error:
          pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
          set(obj.Figure,'CurrentPoint',pt);
       end
       function updateSceneLightPosition(obj)
           if obj.SceneLightLinked, set(obj.SceneLight,'Position',get(obj.RenderAxes,'CameraPosition'));end
       end
    end
    methods % INTERFACING
        function im = captureImage(obj)
           f = getframe(obj.Figure);
           [im,map] = frame2im(f);
           if isempty(map), return; end
           im = ind2rgb(im,map);
        end
        function openSelectionMode(obj)
            switch obj.SelectionMode
                case 'none'
                    return;
                case 'landmark'
                    obj.LandmarkSelection = shape3D;
                    obj.LandmarkSelection.Visible = false;
                    obj.LandmarkSelection.RenderAxes = obj.RenderAxes;
                    obj.LandmarkSelection.SingleColor = [1 0 0];
                    obj.LandmarkSelection.ViewMode = 'Points';
                    obj.LandmarkSelection.VertexSize = 20;
                    obj.LandmarkSelection.Tag = 'private';
                    
                    obj.CamProjection = 'orthographic';
                case 'brush'
                    obj.AreaSelection = shape3D;
                    obj.AreaSelection.Visible = false;
                    obj.AreaSelection.SingleColor = [1 0 0];
                    obj.AreaSelection.RenderAxes = obj.RenderAxes;
                    obj.AreaSelection.ViewMode = 'wireframe';
                    obj.AreaSelection.Tag = 'private';
                    obj.AreaSelectionIndex = [];
                    
                    obj.AreaSelectionDummy = shape3D;
                    obj.AreaSelectionDummy.Visible = false;
                    obj.AreaSelectionDummy.SingleColor = [0 1 0];
                    obj.AreaSelectionDummy.RenderAxes = obj.RenderAxes;
                    obj.AreaSelectionDummy.ViewMode = 'wireframe';
                    obj.AreaSelectionDummy.Tag = 'private';
                    obj.AreaSelectionDummyIndex = [];
                    
                    obj.CamProjection = 'orthographic';
                case 'fill'
                    obj.AreaSelection = shape3D;
                    obj.AreaSelection.Visible = false;
                    obj.AreaSelection.SingleColor = [1 0 0];
                    obj.AreaSelection.RenderAxes = obj.RenderAxes;
                    obj.AreaSelection.ViewMode = 'wireframe';
                    
                    obj.CamProjection = 'orthographic';
                case 'full'
                    return;
                otherwise
                    return;
            end        
        end
        function closeSelectionMode(obj)
            switch obj.SelectionMode
                case 'none'
                    return;
                case 'landmark'
                    if ~isempty(obj.LandmarkSelection), delete(obj.LandmarkSelection);end
                case 'brush'
                    if ~isempty(obj.AreaSelection), delete(obj.AreaSelection);end
                    if ~isempty(obj.AreaSelectionDummy), delete(obj.AreaSelectionDummy);end
                    obj.AreaSelectionIndex = [];
                    obj.AreaSelectionDummyIndex = [];
                case 'fill'
                    if ~isempty(obj.AreaSelection), delete(obj.AreaSelection);end
                    obj.AreaSelectionIndex = [];
                case 'full'
                    return;
                otherwise
                    return;
            end
        end
        function updateArea(obj,area,Vindex)
            if isempty(area), return; end
            if length(Vindex)<2, area.Visible = false; return; end
            shape = obj.CurrentShape;
            Findex = Vindex2Findex(shape,Vindex);
            [~,LOC] = ismember(shape.Faces,Vindex); 
            area.Faces = LOC(Findex,:);
            area.Vertices = shape.Vertices(Vindex,:);
            if~area.Visible, area.Visible = true;end
        end
    end
    methods % DELETE, HIDE and SHOW
        function delete(obj)
           if ishandle(obj.Figure), delete(obj.Figure);end
        end
        function deletePatchChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            if isempty(patches), return; end
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                if superHandleClass.isH(UD), delete(UD); end
            end 
        end
        function hidePatchChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            if isempty(patches), return; end
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                if isa(UD,'shape3D'), UD.Visible = false;end
            end
        end
        function showPatchChildren(obj)
            patches = findobj(get(obj.RenderAxes,'Children'),'flat','Type','patch');
            if isempty(patches), return; end
            for p=1:1:length(patches)
                UD = get(patches(p),'UserData');
                if isa(UD,'shape3D'), UD.Visible = true;end
            end
        end
    end
end % classdef
