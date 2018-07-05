function updateAreaSelection(obj,action,varargin)
        switch lower(action)
            case 'add area'
                if isempty(varargin{1}),return; end
                obj.AreaSelectionIndex = union(obj.AreaSelectionIndex, varargin{1});
            case 'remove area'
                if isempty(varargin{1}),return; end
                obj.AreaSelectionIndex = setdiff(obj.AreaSelectionIndex, varargin{1});
            case 'intersect area'
                if isempty(varargin{1}),return; end
                obj.AreaSelectionIndex = intersect(obj.AreaSelectionIndex, varargin{1});
            case 'set area'
                if isempty(varargin{1}),return; end
                obj.AreaSelectionIndex = varargin{1};
            case 'delete'
                if ~isempty(obj.CurrentShape)
                    crop(obj.CurrentShape,'VertexIndex',obj.AreaSelectionIndex,'Action','delete');
                end
                obj.AreaSelectionIndex = [];
            case 'crop'
                if isempty(obj.AreaSelectionIndex), return;end% nothing is selected
                if ~isempty(obj.CurrentShape)
                   crop(obj.CurrentShape,'VertexIndex',obj.AreaSelectionIndex,'Action','crop');
                end
                obj.AreaSelectionIndex = [];
            case 'clear'
                obj.AreaSelectionIndex = [];
            case 'invert'
                if ~isempty(obj.CurrentShape)
                   fullindex = (1:obj.CurrentShape.nVertices);
                   obj.AreaSelectionIndex = setdiff(fullindex,obj.AreaSelectionIndex);
                end
            case 'store2userdata'
                obj.CurrentShape.UserData.AreaSelectionIndex = obj.AreaSelectionIndex;
                disp('Area Indexed stored in UserData');
                %updateAreaSelection(obj,'clear');
                %obj.SelectionActive = false;
            otherwise
                return;
        end
end
