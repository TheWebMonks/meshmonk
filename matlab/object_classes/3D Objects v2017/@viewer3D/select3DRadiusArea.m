function [areaindex,p] = select3DRadiusArea(obj)
    % selects a connected in terms of triangles area
    areaindex = [];
    p  = select3DPoint(obj);
    if isempty(p), return; end
    distances = sqrt(sum((repmat(p',obj.CurrentShape.nVertices,1)-obj.CurrentShape.Vertices).^2,2));
    areaindex = find(distances<=obj.SelectionSphereRadius);
end
