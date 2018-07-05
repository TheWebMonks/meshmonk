function updateAreaSelectionRadius(obj, action, varargin)
        switch action
            case 'Add Area'
                if isempty(varargin{1}), updateAreaSelectionRadius(obj,'Clear'), return; end
                obj.AreaSelectionDummyIndex = varargin{1};
            case 'Clear'
                obj.AreaSelectionDummyIndex = [];
            otherwise
        end
end