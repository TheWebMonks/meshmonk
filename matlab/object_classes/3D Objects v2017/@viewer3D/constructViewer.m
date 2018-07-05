function constructViewer(obj)
%% Creation
% --- FIGURE CREATION -------------------------------------
obj.Figure = figure(	'Tag', '3D Viewer', ...
	'Units', 'characters', ...
	'Position', [103 24 100 30], ...
	'Name', obj.Tag, ...
	'MenuBar', 'none', ...
    'renderer','opengl',...
	'NumberTitle', 'off', ...
	'Color', [0.15    0.15    0.15], ...
	'Resize', 'on',...
    'CloseRequestFcn',@figure_close_Callback,...
    'WindowButtonDownFcn',@BtnDown,...
    'WindowButtonMotionFcn', @BtnMotion,...
    'WindowButtonUpFcn',@BtnUp,...
    'WindowKeyPressFcn',@KeyPress,...
    'WindowScrollWheelFcn',@ScrollWeel,...
    'UserData',obj);% important is that the viewer3D object is stored in the Figure Userdata, to reconnect;
% --- AXES CREATION -------------------------------------
obj.RenderAxes = axes('Parent', obj.Figure, ...
	'Tag', 'render_axes', ...
	'UserData', obj, ...% important is that the viewer3D object is stored in the Axes Userdata, to reconnect;
	'Units', 'normalized', ...
    'Position',[0.0946428571428572 0.0761904761904762 0.805357142857143 0.835714285714286],...
    'Color',[0.15 0.15 0.15],...
    'XColor',[1 1 1],...
    'YColor',[1 1 1],...
    'ZColor',[1 1 1]);
xlabel(obj.RenderAxes,'X');
ylabel(obj.RenderAxes,'Y');
zlabel(obj.RenderAxes,'Z');
axis(obj.RenderAxes,'off','equal');
hold(obj.RenderAxes,'on');
% --- COLORBAR CREATION ------------------------------------
obj.Colorbar = colorbar('peer',obj.RenderAxes,'Visible','off','Color',[1 1 1],'Location','east');
% --- SCENE LIGHT CONTRUCTION ------------------------------
obj.SceneLight = camlight('headlight','infinite');
set(obj.SceneLight,'Visible','off','Color',[1 1 1]);
% --- TOOLBAR CREATION -------------------------------------
constructToolbar(obj);
% --- CONTEXT MENU CREATION --------------------------------
constructContextMenu(obj);
set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
end

%% CALLBACKS

function figure_close_Callback(hObject,eventdata) %#ok<*INUSD>
         obj = get(hObject,'UserData');
         delete(obj);
end

function KeyPress(varargin)
 obj = get(varargin{1},'UserData');
 if obj.MuteCallbacks, return;end
 if strcmp(obj.Status,'Busy'), return; end
 switch obj.SelectionMode
     case 'landmark'
         switch varargin{2}.Key
             case 'space'
                 obj.SelectionActive = ~obj.SelectionActive;
                 switch obj.SelectionActive
                     case true
                          if isempty(obj.CurrentShape) 
                             b = msgbox(obj.NoCurrentShapeMsg,'Warning');
                             waitfor(b);
                             obj.SelectionActive = ~obj.SelectionActive;
                             return; 
                          end
                          screenPoint2currentPoint(obj)
                          set(obj.Figure,'UIContextMenu',[]);
                     case false
                          set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
                 end
             case 'backspace'
                 updateLMSelection(obj,'Clear'); 
             case {'u' '1'}
                 updateLMSelection(obj,'Store2UserData');
             case 'delete'
                 updateLMSelection(obj,'delete last lm');
             otherwise
                 otherKeys(varargin{2},obj);
         end
     case {'brush' 'fill'}
         switch varargin{2}.Key
             case 'space'
                 obj.SelectionActive = ~obj.SelectionActive;
                 switch obj.SelectionActive
                     case true
                         if isempty(obj.CurrentShape) 
                             b = msgbox(obj.NoCurrentShapeMsg,'Warning');
                             waitfor(b);
                             obj.SelectionActive = ~obj.SelectionActive;
                             return; 
                         end
                         set(obj.Figure,'UIContextMenu',[]);
                         % special for brush mode
                         if strcmp(obj.SelectionMode,'brush') 
                            obj.Action = 'Brush Dummy';
                            updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
                            set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
                         end
                     case false
                         set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
                         % special for brush mode
                         if strcmp(obj.SelectionMode,'brush') 
                             updateAreaSelectionRadius(obj,'Clear');
                             obj.Action = 'none';
                             set(obj.Figure,'WindowButtonMotionFcn',[]);
                         end
                 end
             case 'delete'
                 tmp = varargin{2}.Modifier;
                 if isempty(tmp)
                    mode = 'normal';
                 else
                    mode = tmp{1};
                 end
                 switch mode
                     case 'shift'
                         updateAreaSelection(obj,'Crop');
                     otherwise
                         updateAreaSelection(obj,'Delete');
                 end
             case {'return' '1'}
                 updateAreaSelection(obj,'Crop');
             case 'backspace'
                 updateAreaSelection(obj,'Clear');
             case 'i'
                 updateAreaSelection(obj,'Invert');
             case 'u'
                 updateAreaSelection(obj,'Store2UserData');
             otherwise
                 otherKeys(varargin{2},obj);
         end
     case 'full'
         switch varargin{2}.Key
             case 'space'
                 obj.SelectionActive = ~obj.SelectionActive;
                 switch obj.SelectionActive
                     case true
                         set(obj.Figure,'UIContextMenu',[]);
                     case false
                         set(obj.Figure,'UIContextMenu',obj.ContextMenu.h);
                 end
             case 'delete'
                 tmp = varargin{2}.Modifier;
                 if isempty(tmp)
                    mode = 'normal';
                 else
                    mode = tmp{1};
                 end
                 switch mode
                     case 'shift'
                         updateFullSelection(obj,'Crop');
                     otherwise
                         updateFullSelection(obj,'Delete');
                 end
             case 'backspace'
                 updateFullSelection(obj,'Clear');
             case 'i'
                 updateFullSelection(obj,'Invert');
             otherwise
                 otherKeys(varargin{2},obj);
         end
     otherwise
         otherKeys(varargin{2},obj);
 end
end

function otherKeys(input,obj)
       switch input.Key
           case 'l'
               obj.SceneLightVisible = ~obj.SceneLightVisible;
           case 'k'
               obj.SceneLightLinked = ~obj.SceneLightLinked;
           case 'f1'
               updateFullSelection(obj,'Rendering Surface','Solid');
           case 'f2'
               updateFullSelection(obj,'Rendering Surface','Wireframe');
           case 'f3'
               updateFullSelection(obj,'Rendering Surface','Points');
           case 'f4'
               updateFullSelection(obj,'Rendering Surface','SolidWireframe');
           case 'f5'
               updateFullSelection(obj,'Rendering Color','Texture');
           case 'f6'
               updateFullSelection(obj,'Rendering Color','Indexed');
           case 'f7'
               updateFullSelection(obj,'Rendering Color','Single');
           case 'f8'
               updateFullSelection(obj,'Rendering Material','Facial');
           case 'f9'
               updateFullSelection(obj,'Rendering Material','Dull');
           case 'f10'
               updateFullSelection(obj,'Rendering Material','Shiny');
           case 'f11'
               updateFullSelection(obj,'Rendering Material','Metal');
           case 'f12'
               updateFullSelection(obj,'Rendering Material','Default');
           case {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'numpad0' 'numpad1' 'numpad2' 'numpad3' 'numpad4' 'numpad5' 'numpad6' 'numpad7' 'numpad8' 'numpad9'}
                updateFullSelection(obj,'Rendering Transparancy',(10-str2double(input.Key(end)))/10);
           case 'escape'
               delete(obj);% close the viewer
           case 'v' % print index of closest vertex
               [p,v,vi]  = select3DPoint(obj);
               disp(num2str(vi));
           otherwise
               return
       end
end

function BtnDown(varargin)
     obj = get(varargin{1},'UserData');
     pt = hgconvertunits(obj.Figure,[0 0 get(obj.Figure,'CurrentPoint')],get(obj.Figure,'Units'),'pixels',0);
     obj.MotionData.CurrentPoint = pt(3:4);
     switch obj.SelectionActive
         case true
             BtnDownSelectionMode(obj);
         case false
             if strcmp(obj.Mode,'none'), return; end
             BtnDownMode(obj);
     end
end

function BtnDownSelectionMode(obj)
     switch obj.SelectionMode
        case 'landmark'
                 p  = select3DPoint(obj);
                 if isempty(p), return; end
                 switch get(obj.Figure, 'selectiontype')
                     case 'normal'
                         updateLMSelection(obj,'Add LM',p);
                     case 'alt'
                         updateLMSelection(obj,'Delete LM',p);
                     case 'extend'
                         updateLMSelection(obj,'Show ID',p);
                 end
        case 'fill'
                switch get(obj.Figure, 'selectiontype')
                     case 'normal'% Connected add selection
                          obj.Status = 'Busy';drawnow;
                          areaindex = select3DConnectedArea(obj);
                          if ~isempty(areaindex), updateAreaSelection(obj,'Add Area',areaindex); end
                          obj.Status = 'Ready';
                     case 'alt'% Connected remove selection
                          obj.Status = 'Busy';drawnow;
                          areaindex = select3DConnectedArea(obj);
                          if ~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
                          obj.Status = 'Ready';
                     case 'extend'
                         updateAreaSelection(obj,'Invert'); 
                     case 'open'
                         updateAreaSelection(obj,'Clear');
                     otherwise
                         return;
                end  
        case 'brush'
                 switch get(obj.Figure, 'selectiontype')
                     case 'normal'% Radius add selection
                          updateAreaSelectionRadius(obj,'Clear');
                          areaindex = select3DRadiusArea(obj);
                          if ~isempty(areaindex), updateAreaSelection(obj,'Add Area',areaindex); end
                          obj.Action = 'Brush Selection';
                          set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
                     case 'alt'% Radius remove selection
                          areaindex = select3DRadiusArea(obj);
                          if ~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
                          obj.Action = 'Brush Deselection';
                          obj.AreaSelectionDummy.SingleColor = [0 0 0.502];
                          set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);                      
                     case 'extend'
                          updateAreaSelection(obj,'Invert');
                     case 'open'
                          updateAreaSelection(obj,'Clear');
                     otherwise
                         return;
                 end
         case 'full'
                 selection = get(obj.Figure, 'CurrentObject');
                 if ~strcmp(get(selection,'Type'),'patch'), return; end
                 shapeobject = get(selection,'UserData');
                 if ~isa(shapeobject,'shape3D'), return; end
                 switch get(obj.Figure, 'Selectiontype')
                     case 'normal'% full scan Selection
                         shapeobject.Selected = true;
                         obj.Action = 'Full Scan Selected';                                 
                     case 'alt'% full scan deselection
                         shapeobject.Selected = false;          
                     case 'extend'% ?
                     case 'open'% ? 
                     otherwise
                 end
        otherwise
            return
    end
end

function BtnDownMode(obj)
      switch obj.Mode
         case 'camera'
             axis(obj.RenderAxes,'vis3d');
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'
                     obj.Action = 'rotate camera';
                 case 'alt'
                     obj.Action = 'pan camera';
                 case 'extend'
                     obj.Action = 'zoom camera';
                 case 'open'
                     resetCamera(obj);
                     return;
                 otherwise
                     return;
             end
         case 'light'
             axis(obj.RenderAxes,'vis3d');
             if ~obj.SceneLightVisible, return; end
             switch get(obj.Figure, 'selectiontype')
                 case 'normal'% free light rotation
                     obj.Action = 'rotate light';
                 case 'alt'% change light mode
                     switch obj.SceneLightMode
                         case 'infinite'
                             obj.SceneLightMode = 'local';
                         case 'local'
                             obj.SceneLightMode = 'infinite';                         
                         otherwise
                             return;
                     end
                     return;
                 case 'extend'
                     switch obj.SceneLightRotMode
                         case 'vertical'
                             obj.SceneLightRotMode = 'horizontal';
                         case 'horizontal'
                             obj.SceneLightRotMode = 'vertical';
                         otherwise
                             return;
                     end
                     return
                 case 'open'
                     resetLight(obj);
                     return;
                 otherwise
                     return;
             end
         otherwise
             return;
     end
     set(obj.Figure,'WindowButtonMotionFcn',@BtnMotion);
end

function BtnUp(varargin)
         obj = get(varargin{1},'UserData');
         switch obj.Action
             case 'Full Scan Selected'
                 if obj.SelectionActive
                    scrn_pt = get(0, 'PointerLocation');
                    if ~PointerOnFigure(obj.Figure,scrn_pt)
                       figure_list = findobj('Type','figure');
                       for f=1:1:length(figure_list)
                           if ~PointerOnFigure(figure_list(f),scrn_pt), continue; end
                           UD = get(figure_list(f),'UserData');
                           if ~isa(UD,'viewer3D'), continue; end
                           if (UD==obj), break; end
                           selection = get(obj.Figure, 'CurrentObject');
                           shapeobject = get(selection,'UserData');
                           shapeobject.RenderAxes = UD.RenderAxes;
                           UD.SelectionActive = false;
                           obj.SelectionActive = false;
                           break;
                       end
                    end                        
                 end
             otherwise
         end
         % special for active brush mode
         if strcmp(obj.SelectionMode,'brush')&&obj.SelectionActive
            obj.Action = 'Brush Dummy';
            obj.AreaSelectionDummy.SingleColor = [0 0.7 0];
            updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
         else
             obj.Action = 'none';
             set(obj.Figure,'WindowButtonMotionFcn','');
         end   
         drawnow;
end

function BtnMotion(varargin)
         obj = get(varargin{1},'UserData');
         if strcmp(obj.Action,'none'), return; end
         pt = hgconvertunits(obj.Figure,[0 0 get(obj.Figure,'CurrentPoint')],get(obj.Figure,'Units'),'pixels',0);
         CurrentPoint = pt(3:4);
         try
           deltaPix  = CurrentPoint-obj.MotionData.CurrentPoint;
         catch
           deltaPix = 0;
         end
         obj.MotionData.CurrentPoint = CurrentPoint;
         switch obj.Action
             case 'rotate camera'
                 rotateCamera(obj,deltaPix);
             case 'pan camera'
                 panCamera(obj,deltaPix);
             case 'zoom camera'
                 zoomCamera(obj,deltaPix,[]);
             case 'rotate light'
                 rotateLight(obj,deltaPix);
             case 'Brush Selection'
                 if ~obj.SelectionActive, return; end
                 updateAreaSelection(obj,'Add Area',select3DRadiusArea(obj));
             case 'Brush Deselection'
                 if ~obj.SelectionActive, return; end
                 areaindex = select3DRadiusArea(obj);
                 updateAreaSelectionRadius(obj,'Add Area', areaindex); 
                 if~isempty(areaindex), updateAreaSelection(obj,'Remove Area',areaindex); end
             case 'Brush Dummy'
                 if ~obj.SelectionActive, return; end
                 updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
             otherwise
                 return;
         end
end

function ScrollWeel(varargin)
         obj = get(varargin{1},'UserData');
         if strcmp(obj.SelectionMode,'brush')&& obj.SelectionActive
            if varargin{2}.VerticalScrollCount < 0
                 obj.SelectionSphereRadius  = 1.3*obj.SelectionSphereRadius;
            else
                 obj.SelectionSphereRadius  = 0.7*obj.SelectionSphereRadius;
            end
            updateAreaSelectionRadius(obj,'Add Area', select3DRadiusArea(obj));
            return; 
         end
         if strcmp(obj.SelectionMode,'fill')&& obj.SelectionActive
            if isempty(obj.SelectionDistances), return;end 
            if varargin{2}.VerticalScrollCount < 0
                 obj.SelectionTH  = 1.1*obj.SelectionTH;
            else
                 obj.SelectionTH   = 0.9*obj.SelectionTH ;
            end
            updateAreaSelection(obj,'Set Area', find(obj.SelectionDistances<=obj.SelectionTH));
            return; 
         end
         if obj.SelectionActive, return; end
         switch obj.Mode
             case 'camera'
                 if varargin{2}.VerticalScrollCount < 0
                    q = 1.2;
                    obj.Action = 'zoom in';
                 else
                    q = 0.8;
                    obj.Action = 'zoom out';
                 end
                 zoomCamera(obj,[],q);
             case 'light'
                 if ~obj.SceneLightVisible, return; end
                 xy = [0 0];
                 switch obj.SceneLightRotMode
                     case 'vertical'
                         xy(2) = varargin{2}.VerticalScrollCount*10;
                     case 'horizontal'
                         xy(1) = varargin{2}.VerticalScrollCount*10;
                     otherwise
                         return;
                 end
                 rotateLight(obj,xy);
             otherwise
                 return
         end
         obj.Action = 'none';
end

function val = PointerOnFigure(fig,scrn_pt)
         set(fig,'Units','pixels')
         pos = get(fig,'Position');
         if (scrn_pt(1)>=pos(1))&&(scrn_pt(2)>=pos(2))&&(scrn_pt(1)-pos(1)<=pos(3))&&(scrn_pt(2)-pos(2)<=pos(4))
             val = true;
         else           
             val = false;
         end
end
