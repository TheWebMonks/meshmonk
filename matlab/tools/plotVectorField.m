function obj = plotVectorField(positions,vectors,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% create a shape3D with points at the start and end of the vector

obj = shape3D;
obj.Vertices = [positions;positions+vectors];

numEdges = size(positions,1);
obj.Faces = nan(numEdges,3);
obj.Faces(:,1)= 1:numEdges;
obj.Faces(:,2) = (numEdges+1):numEdges*2;
obj.Faces(:,3) = obj.Faces(:,1); % redundant extra edge

obj.ViewMode = 'wireframe';
viewer(obj,v);

end

