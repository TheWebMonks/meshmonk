function updateToolbar(obj)
    obj.MuteCallbacks = true;% disable callbacks to avoid loops
    if strcmpi(obj.Mode,'camera'), obj.Toolbar.cam_mode.State = 'on'; else, obj.Toolbar.cam_mode.State = 'off';end
    if strcmpi(obj.CamProjection,'orthographic'), obj.Toolbar.ortho.State = 'on'; else, obj.Toolbar.ortho.State = 'off';end
    if strcmpi(obj.CamProjection,'perspective'), obj.Toolbar.persp.State = 'on'; else, obj.Toolbar.persp.State = 'off';end
    if strcmpi(obj.Mode,'light'), obj.Toolbar.light_mode.State = 'on'; else, obj.Toolbar.light_mode.State = 'off';end
    if strcmp(obj.SceneLight.Visible,'on'), obj.Toolbar.light_toggle.State = 'on'; else, obj.Toolbar.light_toggle.State = 'off';end
    if obj.SceneLightLinked, obj.Toolbar.link_toggle.State = 'on'; else, obj.Toolbar.link_toggle.State = 'off';end
    if strcmpi(obj.SelectionMode,'full'), obj.Toolbar.full_mode.State = 'on'; else, obj.Toolbar.full_mode.State = 'off';end
    if strcmpi(obj.SelectionMode,'landmark'), obj.Toolbar.landmark_mode.State = 'on'; else, obj.Toolbar.landmark_mode.State = 'off';end
    if strcmpi(obj.SelectionMode,'brush'), obj.Toolbar.brush_mode.State = 'on'; else, obj.Toolbar.brush_mode.State = 'off';end 
    if strcmpi(obj.SelectionMode,'fill'), obj.Toolbar.fill_mode.State = 'on'; else, obj.Toolbar.fill_mode.State = 'off';end
    if strcmp(obj.Colorbar.Visible,'on'), obj.Toolbar.colorbar.State = 'on'; else, obj.Toolbar.colorbar.State = 'off';end
    obj.MuteCallbacks = false;% enalble callbacks again
end