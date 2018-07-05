function constructContextMenu(obj)
    obj.ContextMenu.h = uicontextmenu;
    item1 = uimenu(obj.ContextMenu.h, 'Label', 'Movement', 'Callback', []);
        obj.ContextMenu.CameraMode = uimenu(item1,'Label','Camera', 'Callback', {@ModeCallback,'camera'},'UserData',obj);
        obj.ContextMenu.LightMode = uimenu(item1,'Label','Light', 'Callback', {@ModeCallback,'light'},'UserData',obj);
    item2 = uimenu(obj.ContextMenu.h, 'Label', 'Camera', 'Callback', []);
        obj.ContextMenu.Ortho = uimenu(item2, 'Label', 'Orthographic', 'Callback', {@ProjectionCallback, 'orthographic'},'UserData',obj);
        obj.ContextMenu.Persp = uimenu(item2, 'Label', 'Perspective', 'Callback', {@ProjectionCallback 'perspective'},'UserData',obj);
        obj.ContextMenu.ResetCamera = uimenu(item2, 'Label', 'Reset', 'Callback', @ResetCameraCallback,'UserData',obj);
    item3 = uimenu(obj.ContextMenu.h, 'Label', 'Light', 'Callback', []);
        obj.ContextMenu.LightVisible = uimenu(item3, 'Label', 'Visible (l)', 'Callback', @LightVisibleCallback,'UserData',obj);
        obj.ContextMenu.LinkLight = uimenu(item3, 'Label', 'Link (k)', 'Callback', @LinkLightCallback,'UserData',obj);
        obj.ContextMenu.LightColor = uimenu(item3, 'Label', 'Color', 'Callback', @LightColorCallback,'UserData',obj);
    item4 = uimenu(obj.ContextMenu.h, 'Label', 'Axes', 'Callback', []);
        obj.ContextMenu.AxesVisible = uimenu(item4, 'Label', 'Visible', 'Callback', @AxesVisibleCallback,'UserData',obj);
        obj.ContextMenu.AxesGrid = uimenu(item4, 'Label', 'Grid', 'Callback', @AxesGridCallback,'UserData',obj);
        obj.ContextMenu.AxesBox = uimenu(item4, 'Label', 'Box', 'Callback', @AxesBoxCallback,'UserData',obj);
        item41 = uimenu(item4, 'Label', 'Color', 'Callback', []);
            obj.ContextMenu.Background = uimenu(item41, 'Label', 'Background', 'Callback', @BackgroundCallback,'UserData',obj);
            obj.ContextMenu.Wall = uimenu(item41, 'Label', 'Wall', 'Callback', @WallCallback,'UserData',obj);
            obj.ContextMenu.XColor = uimenu(item41, 'Label', 'XColor', 'Callback', @XColorCallback,'UserData',obj);
            obj.ContextMenu.YColor = uimenu(item41, 'Label', 'YColor', 'Callback', @YColorCallback,'UserData',obj);
            obj.ContextMenu.ZColor = uimenu(item41, 'Label', 'ZColor', 'Callback', @ZColorCallback,'UserData',obj);
        item42 = uimenu(item4, 'Label', 'Renderer', 'Callback', []);
            obj.ContextMenu.Zbuffer = uimenu(item42, 'Label', 'Zbuffer', 'Callback', {@RendererCallback,'zbuffer'}, 'checked', 'on','UserData',obj);
            obj.ContextMenu.OpenGl = uimenu(item42, 'Label', 'OpenGl', 'Callback', {@RendererCallback,'opengl'},'UserData',obj);
    item5 = uimenu(obj.ContextMenu.h, 'Label', 'Selection', 'Callback', []);
        item51 = uimenu(item5, 'Label', 'Mode', 'Callback', []);
            obj.ContextMenu.FullMode = uimenu(item51, 'Label', 'Full Scan', 'Callback', {@SelectionModeCallback,'full'},'UserData',obj,'Checked','off');
            obj.ContextMenu.LandmarkMode = uimenu(item51, 'Label', 'Landmark', 'Callback', {@SelectionModeCallback,'landmark'},'UserData',obj);
            obj.ContextMenu.BrushMode = uimenu(item51, 'Label', 'Brush', 'Callback', {@SelectionModeCallback,'brush'},'UserData',obj);
            obj.ContextMenu.FillMode = uimenu(item51, 'Label', 'Fill', 'Callback', {@SelectionModeCallback,'fill'},'UserData',obj);
        item52 = uimenu(item5, 'Label', 'Edit', 'Callback', []);
            obj.ContextMenu.ClearSelection = uimenu(item52, 'Label', 'Clear (Backspace)', 'Callback', @ClearSelectionCallback,'UserData',obj);
            obj.ContextMenu.InvertSelection = uimenu(item52, 'Label', 'Invert (i)', 'Callback', @InvertSelectionCallback,'UserData',obj);
            obj.ContextMenu.DeleteSelection = uimenu(item52, 'Label', 'Delete (Del)', 'Callback', @DeleteSelectionCallback,'UserData',obj);
            obj.ContextMenu.CropSelection = uimenu(item52, 'Label', 'Crop (alt+Del)', 'Callback', @CropSelectionCallback,'UserData',obj);
        obj.ContextMenu.StoreIntoUserData = uimenu(item5, 'Label', 'Store2UserData (u)', 'Callback', @StoreIntoUserDataCallback,'UserData',obj);

    item6 = uimenu(obj.ContextMenu.h, 'Label', 'Rendering', 'Callback', []);
        item61 = uimenu(item6, 'Label', 'Surface', 'Callback', []);
            obj.ContextMenu.SurfaceSolid = uimenu(item61,'Label','Solid (F1)','Callback',{@SurfaceCallback, 'Solid'},'UserData',obj);
            obj.ContextMenu.SurfaceWireframe = uimenu(item61,'Label','Wireframe (F2)','Callback',{@SurfaceCallback, 'Wireframe'},'UserData',obj);
            obj.ContextMenu.SurfacePoints = uimenu(item61,'Label','Points (F3)','Callback',{@SurfaceCallback, 'Points'},'UserData',obj);
            obj.ContextMenu.SurfaceSolidWireframe = uimenu(item61,'Label','Solid/Wireframe (F4)','Callback',{@SurfaceCallback, 'Solid/Wireframe'},'UserData',obj);
        item62 = uimenu(item6, 'Label', 'Color', 'Callback', []);
            obj.ContextMenu.ColorTexture = uimenu(item62,'Label','Texture (F5)','Callback',{@ColorCallback, 'Texture'},'UserData',obj);
            obj.ContextMenu.ColorIndexed = uimenu(item62,'Label','Indexed (F6)','Callback',{@ColorCallback, 'Indexed'},'UserData',obj);
            obj.ContextMenu.ColorSingle = uimenu(item62,'Label','Single (F7)','Callback',{@ColorCallback, 'Single'},'UserData',obj);
            obj.ContextMenu.SetColorSingle = uimenu(item62,'Label','Single Color','Callback',{@SetColorCallback},'UserData',obj);
        item63 = uimenu(item6, 'Label', 'Material', 'Callback', []);
            obj.ContextMenu.MaterialFacial = uimenu(item63,'Label','Facial (F8)','Callback',{@MaterialCallback,'Facial'},'UserData',obj);
            obj.ContextMenu.MaterialDull = uimenu(item63,'Label','Dull (F9)','Callback',{@MaterialCallback,'Dull'},'UserData',obj);
            obj.ContextMenu.MaterialShiny = uimenu(item63,'Label','Shiny (F10)','Callback',{@MaterialCallback,'Shiny'},'UserData',obj);
            obj.ContextMenu.MaterialMetal = uimenu(item63,'Label','Metal (F11)','Callback',{@MaterialCallback,'Metal'},'UserData',obj);
            obj.ContextMenu.MaterialDefault = uimenu(item63,'Label','Default (F12)','Callback',{@MaterialCallback,'Default'},'UserData',obj);
        item64 = uimenu(item6, 'Label', 'Lighting', 'Callback', []);
            obj.ContextMenu.LightingNone = uimenu(item64,'Label','None','Callback',{@LightingCallback,'none'},'UserData',obj);
            obj.ContextMenu.LightingFlat = uimenu(item64,'Label','Flat','Callback',{@LightingCallback,'flat'},'UserData',obj);
            obj.ContextMenu.LightingGouraud = uimenu(item64,'Label','Gouraud','Callback',{@LightingCallback,'gouraud'},'UserData',obj);
            obj.ContextMenu.LightingPhong = uimenu(item64,'Label','Phong','Callback',{@LightingCallback,'phong'},'UserData',obj);
        item65 = uimenu(item6, 'Label', 'Transparancy', 'Callback', []);
            obj.ContextMenu.Transparancy0 = uimenu(item65,'Label','0 (0)','Callback',{@TransparancyCallback,1},'UserData',obj);
            obj.ContextMenu.Transparancy1 = uimenu(item65,'Label','1 (1)','Callback',{@TransparancyCallback,0.9},'UserData',obj);
            obj.ContextMenu.Transparancy2 = uimenu(item65,'Label','2 (2)','Callback',{@TransparancyCallback,0.8},'UserData',obj);
            obj.ContextMenu.Transparancy3 = uimenu(item65,'Label','3 (3)','Callback',{@TransparancyCallback,0.7},'UserData',obj);
            obj.ContextMenu.Transparancy4 = uimenu(item65,'Label','4 (4)','Callback',{@TransparancyCallback,0.6},'UserData',obj);
            obj.ContextMenu.Transparancy5 = uimenu(item65,'Label','5 (5)','Callback',{@TransparancyCallback,0.5},'UserData',obj);
            obj.ContextMenu.Transparancy6 = uimenu(item65,'Label','6 (6)','Callback',{@TransparancyCallback,0.4},'UserData',obj);
            obj.ContextMenu.Transparancy7 = uimenu(item65,'Label','7 (7)','Callback',{@TransparancyCallback,0.3},'UserData',obj);
            obj.ContextMenu.Transparancy8 = uimenu(item65,'Label','8 (8)','Callback',{@TransparancyCallback,0.2},'UserData',obj);
            obj.ContextMenu.Transparancy9 = uimenu(item65,'Label','9 (9)','Callback',{@TransparancyCallback,0.1},'UserData',obj);
end

%% CALLBACKS
function ModeCallback(hObject,eventdata,arg) %#ok<*INUSL>
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Mode = arg;
end

function ProjectionCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.CamProjection = arg;
end

function ResetCameraCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         resetCamera(obj)
end

function LightVisibleCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightVisible = ~obj.SceneLightVisible;
end

function LinkLightCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightLinked = ~obj.SceneLightLinked;
end

function LightColorCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.SceneLightColor = color; end
         
end

function AxesVisibleCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.AxesVisible = ~obj.AxesVisible;
end

function AxesGridCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.AxesGrid = ~obj.AxesGrid;  
end

function AxesBoxCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.AxesBox = ~obj.AxesBox;
end

function BackgroundCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.BackgroundColor = color; end
end

function WallCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.AxesWallColor = color; end
end

function XColorCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.AxesXColor = color; end
end

function YColorCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.AxesYColor = color; end
end

function ZColorCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), obj.AxesZColor = color; end
end

function RendererCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Renderer = arg;
end

function SelectionModeCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = lower(arg);   
end

function ClearSelectionCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         switch obj.SelectionMode
             case 'landmark'
                 updateLMSelection(obj,'Clear');
             case {'area' 'brush' 'fill'}
                 updateAreaSelection(obj,'Clear');
             case 'full'
                 updateFullScanSelection(obj,'Clear')
             otherwise
         end
end

function InvertSelectionCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         switch obj.SelectionMode
             case 'landmark'
             case {'area' 'brush' 'fill'}
                 updateAreaSelection(obj,'Invert');
             case 'full'
                 updateFullScanSelection(obj,'Invert')
             otherwise
         end
end

function DeleteSelectionCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         switch obj.SelectionMode
             case 'landmark'
                 updateLMSelection(obj,'Clear');
             case {'area' 'brush' 'fill'}
                 updateAreaSelection(obj,'Delete');
             case 'full'
                 updateFullScanSelection(obj,'Delete')
             otherwise
         end
end

function CropSelectionCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         switch obj.SelectionMode
             case 'landmark'
             case {'area' 'brush' 'fill'}
                 updateAreaSelection(obj,'Crop');
             case 'full'
                 updateFullScanSelection(obj,'Crop')
             otherwise
         end
end

function StoreIntoUserDataCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         % DO SOMETHING MORE
         %updateLMSelection(obj,'Link Pose LM');         
end

function SurfaceCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         updateFullSelection(obj,'Rendering Surface',arg);
end

function ColorCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         updateFullSelection(obj,'Rendering Color',arg);
end

function SetColorCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         color = uisetcolor;
         if ~(length(color) == 1), updateFullSelection(obj,'Rendering Single Color',color); end
end

function MaterialCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         updateFullSelection(obj,'Rendering Material',arg);
end

function LightingCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         updateFullSelection(obj,'Rendering Lighting',arg);
end

function TransparancyCallback(hObject,eventdata,arg)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         updateFullSelection(obj,'Rendering Transparancy',arg);
end
