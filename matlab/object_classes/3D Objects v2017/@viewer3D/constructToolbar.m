function constructToolbar(obj)
%% CREATION

load viewertoolbarimagesOBJ.mat;
obj.Toolbar.h = uitoolbar(obj.Figure);
obj.Toolbar.open = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.open,...
                'TooltipString','Open 3D',...
                'Separator','off',...
                'ClickedCallback',@open_Callback,...
                'UserData',obj);
obj.Toolbar.save = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.save,...
                'TooltipString','Save 3D',...
                'Separator','off',...
                'ClickedCallback',@save_Callback,...
                'UserData',obj);           
if strcmpi(obj.Mode,'camera'), state = 'on'; else, state = 'off';end            
obj.Toolbar.cam_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.orbit,...
                'TooltipString','Toggle Camera Mode',...
                'Separator','on',...
                'OnCallback',@camera_mode_on_Callback,...
                'OffCallback',@camera_mode_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.CamProjection,'orthographic'), state = 'on'; else, state = 'off';end             
obj.Toolbar.ortho = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.ortho,...
                'TooltipString','Toggle Orthographic Projection',...
                'Separator','off',...
                'OnCallback',@camera_ortho_on_Callback,...
                'OffCallback',@camera_ortho_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.CamProjection,'perspective'), state = 'on'; else, state = 'off';end
obj.Toolbar.persp = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.perspective,...
                'TooltipString','Toggle Perspective Projection',...
                'Separator','off',...
                'OnCallback',@camera_persp_on_Callback,...
                'OffCallback',@camera_persp_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.Mode,'light'), state = 'on'; else, state = 'off';end            
obj.Toolbar.light_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.orbitlight,...
                'TooltipString','Toggle Light Mode',...
                'Separator','on',...
                'OnCallback',@light_mode_on_Callback,...
                'OffCallback',@light_mode_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmp(obj.SceneLight.Visible,'on'), state = 'on'; else, state = 'off';end
obj.Toolbar.light_toggle = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.light,...
                'TooltipString','Toggle Light',...
                'Separator','off',...
                'OnCallback',@light_on_Callback,...
                'OffCallback',@light_off_Callback,...
                'State',state,...
                'UserData',obj);
if obj.SceneLightLinked, state = 'on'; else, state = 'off';end            
obj.Toolbar.link_toggle = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.link,...
                'TooltipString','Link Light',...
                'Separator','off',...
                'OnCallback',@link_on_Callback,...
                'OffCallback',@link_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.SelectionMode,'full'), state = 'on'; else, state = 'off';end            
obj.Toolbar.full_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.objectcube,...
                'TooltipString','Full Scan Mode',...
                'Separator','on',...
                'OnCallback',@full_mode_on_Callback,...
                'OffCallback',@full_mode_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.SelectionMode,'landmark'), state = 'on'; else, state = 'off';end             
obj.Toolbar.landmark_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.LMindicator,...
                'TooltipString','Landmark Mode',...
                'Separator','off',...
                'OnCallback',@landmark_mode_on_Callback,...
                'OffCallback',@landmark_mode_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.SelectionMode,'brush'), state = 'on'; else, state = 'off';end 
obj.Toolbar.brush_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.brush,...
                'TooltipString','Brush Mode',...
                'Separator','off',...
                'OnCallback',@brush_mode_on_Callback,...
                'OffCallback',@brush_mode_off_Callback,...
                'State',state,...
                'UserData',obj);
if strcmpi(obj.SelectionMode,'fill'), state = 'on'; else, state = 'off';end             
obj.Toolbar.fill_mode = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.fill,...
                'TooltipString','Fill Mode',...
                'Separator','off',...
                'OnCallback',@fill_mode_on_Callback,...
                'OffCallback',@fill_mode_off_Callback,...
                'State',state,...
                'UserData',obj);            
if strcmp(obj.Colorbar.Visible,'on'), state = 'on'; else, state = 'off';end
obj.Toolbar.colorbar = uitoggletool(obj.Toolbar.h,'Cdata',toolbarimages.colorbar,...
                'TooltipString','Toggle Colorbar',...
                'Separator','on',...
                'OnCallback',@colorbar_on_Callback,...
                'OffCallback',@colorbar_off_Callback,...
                'State',state,...
                'UserData',obj);
obj.Toolbar.colorbareditor = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.colorbareditor,...
                'TooltipString','Edit Colorbar',...
                'Separator','off',...
                'ClickedCallback',@colorbar_editor_Callback,...
                'UserData',[]);
obj.Toolbar.imagecapture = uipushtool(obj.Toolbar.h,'Cdata',toolbarimages.camera,...
                'TooltipString','Edit Colorbar',...
                'Separator','on',...
                'ClickedCallback',@imageCaptureCallback,...
                'UserData',obj);
end

%% CALLBACKS
function open_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         shapes = shape3D.import;
         if isempty(shapes), return; end
         if iscell(shapes)
             for i=1:1:length(shapes)
                shapes{i}.RenderAxes = obj.RenderAxes;
                shapes{i}.Visible = true;
                shapes{i}.Selected = true; 
             end
         else
             shapes.RenderAxes = obj.RenderAxes;
             shapes.Visible = true;
             shapes.Selected = true;
         end                        
end

function save_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if isempty(obj.ShapeChildren), return; end
         for i=1:1:length(obj.ShapeChildren)
             if obj.ShapeChildren{i}.Selected
                 export(obj.ShapeChildren{i});
             end
         end                     
end

function camera_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Mode = 'camera';                   
end

function camera_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Mode = 'light';        
end

function camera_ortho_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.CamProjection = 'orthographic';
end

function camera_ortho_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.CamProjection = 'perspective';
end

function camera_persp_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.CamProjection = 'perspective';
end

function camera_persp_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.CamProjection = 'orthographic';
end

function light_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Mode = 'light';              
end

function light_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.Mode = 'camera';
end

function light_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightVisible = true;
end

function light_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightVisible = false;
end

function link_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightLinked = true;
end

function link_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SceneLightLinked = false;
end

function full_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'full';                    
end

function full_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'none';
end

function landmark_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'landmark';         
end

function landmark_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'none';
end

function fill_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'fill';      
end

function fill_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'none';
end

function brush_mode_on_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'brush';          
         
end

function brush_mode_off_Callback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         obj.SelectionMode = 'none';
end

function colorbar_on_Callback(hObject,eventdata)
           obj = get(hObject,'UserData');
           if obj.MuteCallbacks, return;end
           obj.Colorbar.Visible = 'on';
end

function colorbar_off_Callback(hObject,eventdata)
           obj = get(hObject,'UserData');
           if obj.MuteCallbacks, return;end
           obj.Colorbar.Visible = 'off';
end

function colorbar_editor_Callback(hObject,eventdata)
         colormapeditor;
end

function imageCaptureCallback(hObject,eventdata)
         obj = get(hObject,'UserData');
         if obj.MuteCallbacks, return;end
         im = captureImage(obj);
         pause(0.1);
         [filename, pathname, filterindex] = uiputfile( ...
                                             {'*.tiff','Print Images (*.tiff)'; ...
                                              '*.tiff','Images (*.tiff)'; ...
                                              '*.png','Images (*.png)'; ...
                                              '*.tiff','Images (*.tiff)'; ...
                                              '*.bmp','Images (*.bmp)'; ...
                                              '*.jpg','Images (*.jpg)'; ...
                                              '*.gif','Images (*.gif)'},...
                                              'Save as', get(obj.Figure,'Name'));
          if filename == 0; return; end
          cd(pathname);             
          switch filterindex
              case 1
                  print(obj.Figure,'-dtiffn','-r600',filename(1:end-5));
              case 2
                  imwrite(im,filename,'tiff','Compression','none','Resolution',300);
              case 3
                  imwrite(im,filename,'png');
              case 4
                  imwrite(im,filename,'bmp');
              case 5
                  imwrite(im,filename,'jpeg');
              case 6
                  imwrite(im,filename,'gif');
          end
          disp('Image Saved');
end