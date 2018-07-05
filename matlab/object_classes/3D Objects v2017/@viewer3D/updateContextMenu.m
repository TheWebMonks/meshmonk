function updateContextMenu(obj)
    obj.MuteCallbacks = true;% disable callbacks to avoid loops
    if strcmpi(obj.Mode,'camera'), obj.ContextMenu.CameraMode.Checked = 'on'; else, obj.ContextMenu.CameraMode.Checked = 'off';end
    if strcmpi(obj.CamProjection,'orthographic'), obj.ContextMenu.Ortho.Checked = 'on'; else, obj.ContextMenu.Ortho.Checked = 'off';end
    if strcmpi(obj.CamProjection,'perspective'), obj.ContextMenu.Persp.Checked = 'on'; else, obj.ContextMenu.Persp.Checked = 'off';end
    if strcmpi(obj.Mode,'light'), obj.ContextMenu.LightMode.Checked = 'on'; else, obj.ContextMenu.LightMode.Checked = 'off';end
    if strcmp(obj.SceneLight.Visible,'on'), obj.ContextMenu.LightVisible.Checked = 'on'; else, obj.ContextMenu.LightVisible.Checked = 'off';end
    if obj.SceneLightLinked, obj.ContextMenu.LinkLight.Checked = 'on'; else, obj.ContextMenu.LinkLight.Checked = 'off';end
    if obj.AxesVisible, obj.ContextMenu.AxesVisible.Checked = 'on'; else, obj.ContextMenu.AxesVisible.Checked = 'off';end
    if obj.AxesGrid, obj.ContextMenu.AxesGrid.Checked = 'on'; else, obj.ContextMenu.AxesGrid.Checked = 'off';end
    if obj.AxesBox, obj.ContextMenu.AxesBox.Checked = 'on'; else, obj.ContextMenu.AxesBox.Checked = 'off';end
    if strcmpi(obj.Renderer,'opengl'), obj.ContextMenu.OpenGl.Checked = 'on'; else, obj.ContextMenu.OpenGl.Checked = 'off';end
    if strcmpi(obj.Renderer,'zbuffer'), obj.ContextMenu.Zbuffer.Checked = 'on'; else, obj.ContextMenu.Zbuffer.Checked = 'off';end
    if strcmpi(obj.SelectionMode,'full'), obj.ContextMenu.FullMode.Checked = 'on'; else, obj.ContextMenu.FullMode.Checked = 'off';end
    if strcmpi(obj.SelectionMode,'landmark'), obj.ContextMenu.LandmarkMode.Checked = 'on'; else, obj.ContextMenu.LandmarkMode.Checked = 'off';end
    if strcmpi(obj.SelectionMode,'brush'), obj.ContextMenu.BrushMode.Checked = 'on'; else, obj.ContextMenu.BrushMode.Checked = 'off';end
    if strcmpi(obj.SelectionMode,'fill'), obj.ContextMenu.FillMode.Checked = 'on'; else, obj.ContextMenu.FillMode.Checked = 'off';end
    if strcmpi(obj.SelectionMode,'none')
        set(obj.ContextMenu.ClearSelection,'Enable','off');
        set(obj.ContextMenu.InvertSelection,'Enable','off');
        set(obj.ContextMenu.DeleteSelection,'Enable','off');
        set(obj.ContextMenu.CropSelection,'Enable','off');
    elseif strcmpi(obj.SelectionMode,'landmark')
        set(obj.ContextMenu.ClearSelection,'Enable','on');
        set(obj.ContextMenu.InvertSelection,'Enable','off');
        set(obj.ContextMenu.DeleteSelection,'Enable','off');
        set(obj.ContextMenu.CropSelection,'Enable','off');
    else
        set(obj.ContextMenu.ClearSelection,'Enable','on');
        set(obj.ContextMenu.InvertSelection,'Enable','on');
        set(obj.ContextMenu.DeleteSelection,'Enable','on');
        set(obj.ContextMenu.CropSelection,'Enable','on');
    end   
    obj.MuteCallbacks = false;% enalble callbacks again
end