function areaindex = select3DConnectedArea(obj)
    % selects a connected in terms of triangles area
    areaindex = [];
    [p,~,~,~,fi]  = select3DPoint(obj);
    if isempty(p), return; end
    distances = intraDistances(obj.CurrentShape,'VertexIndex',obj.CurrentShape.Faces(fi,:));
    areaindex = find(~isnan(distances));
    obj.SelectionDistances = distances;
    obj.SelectionTH = max(distances);
end
