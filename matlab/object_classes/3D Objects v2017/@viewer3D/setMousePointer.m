function setMousePointer(obj)
  switch obj.Status
      case 'Ready'
        switch obj.SelectionActive
            case true
                 switch obj.SelectionMode
                     case 'landmark'
                         set(obj.Figure,'Pointer','crosshair');
                     case {'area' 'fill' 'brush' 'full'}
                         set(obj.Figure,'Pointer','arrow');
                     otherwise
                 end
            case false
                switch obj.Mode
                    case 'camera'
                        switch obj.Action
                            case 'none'
                                set(obj.Figure,'Pointer','hand');
                            case 'rotate camera'
                                set(obj.Figure,'Pointer','fleur');
                            case 'pan camera'
                                set(obj.Figure,'Pointer','fleur');
                            case 'zoom camera'
                                set(obj.Figure,'Pointer','fleur');
                            case 'zoom in'
                                set(obj.Figure,'Pointer','fleur');
                            case 'zoom out'
                                set(obj.Figure,'Pointer','fleur');
                            otherwise
                                setptr(obj.Figure,'hand');
                        end
                    case 'light'
                        switch obj.Action
                            case 'none'
                                set(obj.Figure,'Pointer','hand');
                            case 'rotate light'
                                set(obj.Figure,'Pointer','fleur');
                            otherwise
                                set(obj.Figure,'Pointer','hand');
                        end                           
                    otherwise
                        set(obj.Figure,'Pointer','arrow');
                end
        end
      case 'Busy'
          set(obj.Figure,'Pointer','watch');
      otherwise
  end
end