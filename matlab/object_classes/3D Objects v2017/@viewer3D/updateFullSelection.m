function updateFullSelection(obj,action,varargin)
    if ~isempty(find(strcmp(varargin,'Busy'),1)),obj.Status = 'Busy';drawnow;end
    children = obj.ShapeChildren;
    for i=1:1:length(children)
        switch action
            case 'Delete'
                if ~children{i}.Selected, continue; end
                delete(children{i});
            case 'Clear'
                if ~children{i}.Selected, continue; end
                children{i}.Selected = false;
            case 'Invert'
                switch children{i}.Selected
                    case true
                       children{i}.Selected = false; 
                    case false
                        children{i}.Selected = true;
                    otherwise
                end
            case 'Crop'
                if children{i}.Selected, continue; end
                delete(children{i});
            case 'Rendering Surface'
                if ~children{i}.Selected, continue; end
                children{i}.ViewMode = varargin{1};
            case 'Rendering Color'
                if ~children{i}.Selected, continue; end
                children{i}.ColorMode = varargin{1};
            case 'Rendering Material'
                if ~children{i}.Selected, continue; end
                children{i}.Material = varargin{1};
            case 'Rendering Lighting'
                if ~children{i}.Selected, continue; end
                children{i}.LightMode = varargin{1};
            case 'Rendering Transparancy'
                if ~children{i}.Selected, continue; end
                children{i}.Alpha = varargin{1};
            case 'Rendering Single Color'
                if ~children{i}.Selected, continue; end
                children{i}.SingleColor = varargin{1};
            otherwise
        end
    end
    obj.Status = 'Ready';
end