function updateLMSelection(obj,action,varargin)
         switch lower(action)
             case 'add lm'
                 point = varargin{1};
                 obj.LandmarkSelection.Vertices = [obj.LandmarkSelection.Vertices; point'];
             case 'delete lm'
                 if obj.LandmarkSelection.Level<1, return; end
                 point = varargin{1};
                 index = knnsearch(obj.LandmarkSelection.Vertices,point');
                 if isempty(index), return; end
                 full_index = 1:obj.LandmarkSelection.nVertices;
                 rest_index = setdiff(full_index,index);
                 obj.LandmarkSelection.Vertices = obj.LandmarkSelection.Vertices(rest_index,:);
             case 'delete last lm'
                 if obj.LandmarkSelection.Level<1, return; end
                 if obj.LandmarkSelection.nVertices == 1,updateLMSelection(obj,'clear'); return;end %patch is going to be invalid
                 obj.LandmarkSelection.Vertices = obj.LandmarkSelection.Vertices(1:end-1,:);
             case 'clear'
                 obj.LandmarkSelection.Vertices = [];
             case 'store2userdata'
                 if isfield(obj.CurrentShape.UserData,'LandmarkSelection')
                    tmp = obj.CurrentShape.UserData.LandmarkSelection;
                    if superHandleClass.isH(tmp), delete(tmp); end
                 end
                 obj.CurrentShape.UserData.LandmarkSelection = clone(obj.LandmarkSelection);
                 obj.CurrentShape.UserData.LandmarkSelection.SingleColor = [0 1 0];
                 obj.CurrentShape.UserData.LandmarkSelection.RenderAxes = obj.RenderAxes;
                 obj.CurrentShape.UserData.LandmarkSelection.Visible = true;
                 updateLMSelection(obj,'clear');
                 obj.SelectionActive = false;
             case 'show id'
                 if obj.LandmarkSelection.Level<1, return; end
                 point = varargin{1};
                 index = knnsearch(obj.LandmarkSelection.Vertices,point');
                 disp(['ID = ' num2str(index)]);
             case 'updated shape'
                 % nothing yet
             otherwise
                 return;
         end
         if ~isempty(obj.LandmarkSelection)
             if obj.LandmarkSelection.Level>0
                obj.LandmarkSelection.Visible = true;
             else
                obj.LandmarkSelection.Visible = false;
             end
         end
         drawnow;
end