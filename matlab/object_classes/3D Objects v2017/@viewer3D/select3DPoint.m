function [pout, vout, viout, facevout, faceiout]  = select3DPoint(obj)
%SELECT3D(H) Determines the selected point in 3-D data space.
%  P = SELECT3D determines the point, P, in data space corresponding 
%  to the current selection position. P is a point on the first 
%  patch or surface face intersected along the selection ray. If no 
%  face is encountered along the selection ray, P returns empty.
%
%  P = SELECT3D(H) constrains selection to graphics handle H and,
%  if applicable, any of its children. H can be a figure, axes, 
%  patch, or surface object.
%
%  [P V] = SELECT3D(...), V is the closest face or line vertex 
%  selected based on the figure's current object. 
%
%  [P V VI] = SELECT3D(...), VI is the index into the object's 
%  x,y,zdata properties corresponding to V, the closest face vertex 
%  selected.
%
%  [P V VI FACEV] = SELECT3D(...), FACE is an array of vertices 
%  corresponding to the face polygon containing P and V. 
%  
%  [P V VI FACEV FACEI] = SELECT3D(...), FACEI is the row index into 
%  the object's face array corresponding to FACE. For patch 
%  objects, the face array can be obtained by doing 
%  get(mypatch,'faces'). For surface objects, the face array 
%  can be obtained from the output of SURF2PATCH (see 
%  SURF2PATCH for more information).
%
%  RESTRICTIONS:
%  SELECT3D supports surface, patch, or line object primitives. For surface 
%  and patches, the algorithm assumes non-self-intersecting planar faces. 
%  For line objects, the algorithm always returns P as empty, and V will
%  be the closest vertex relative to the selection point. 
% Output variables
pout = [];
vout = [];
viout = [];
facevout = [];
faceiout = [];

% other variables
isline = false;
isperspective = false;

screenPoint2currentPoint(obj);
fig = obj.Figure;
ax = obj.RenderAxes;
axchild = obj.CurrentShape.PatchHandle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get vertex, face, and current point data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp = get(ax,'currentpoint')';

vert = get(axchild,'vertices');
faces = get(axchild,'faces');

% Add z if empty
if size(vert,2)==2
   vert(:,3) = zeros(size(vert(:,2)));
end

% NaN and Inf check 
nan_inf_test1 = isnan(faces) | isinf(faces);
nan_inf_test2 = isnan(vert) | isinf(vert);
if any(nan_inf_test1(:)) | any(nan_inf_test2(:))
    warning(sprintf('%s does not support NaNs or Infs in face/vertex data.',mfilename));
end

%% COMMENTEDOUT
% xvert = local_Data2PixelTransform(ax,vert)';
% xcp = local_Data2PixelTransform(ax,cp')';
%% NEW CODE
persistent convertDataSpaceCoordsToViewerCoords;

if verLessThan('matlab','8.4.0')
    xvert = local_Data2PixelTransform(ax,vert)';
    xcp = local_Data2PixelTransform(ax,cp')';
else
    if isempty(convertDataSpaceCoordsToViewerCoords)
        convertDataSpaceCoordsToViewerCoords = specgraphhelper('convertDataSpaceCoordsToViewerCoords');
    end
    xvert = convertDataSpaceCoordsToViewerCoords(axchild, vert');
    xcp = convertDataSpaceCoordsToViewerCoords(ax, cp);
end

% Translate vertices so that the selection point is at the origin.
xvert(1,:) = xvert(1,:) - xcp(1,2);
xvert(2,:) = xvert(2,:) - xcp(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform 2-D crossing test (Jordan Curve Theorem) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all vertices that have y components less than zero
vert_with_negative_y = zeros(size(faces));
face_y_vert = xvert(2,faces);
ind_vert_with_negative_y = find(face_y_vert<0); 
vert_with_negative_y(ind_vert_with_negative_y) = true;

% Find all the line segments that span the x axis
is_line_segment_spanning_x = abs(diff([vert_with_negative_y, vert_with_negative_y(:,1)],1,2));

% Find all the faces that have line segments that span the x axis
ind_is_face_spanning_x = find(any(is_line_segment_spanning_x,2));

% Ignore data that doesn't span the x axis
candidate_faces = faces(ind_is_face_spanning_x,:);
vert_with_negative_y = vert_with_negative_y(ind_is_face_spanning_x,:);
is_line_segment_spanning_x = is_line_segment_spanning_x(ind_is_face_spanning_x,:);

% Create line segment arrays
pt1 = candidate_faces;
pt2 = [candidate_faces(:,2:end), candidate_faces(:,1)];

% Point 1
x1 = reshape(xvert(1,pt1),size(pt1));
y1 = reshape(xvert(2,pt1),size(pt1));

% Point 2
x2 = reshape(xvert(1,pt2),size(pt2));
y2 = reshape(xvert(2,pt2),size(pt2));

% Cross product of vector to origin with line segment
cross_product_test = -x1.*(y2-y1) > -y1.*(x2-x1);

% Find all line segments that cross the positive x axis
crossing_test = (cross_product_test==vert_with_negative_y) & is_line_segment_spanning_x;

% If the number of line segments is odd, then we intersected the polygon
s = sum(crossing_test,2);
s = mod(s,2);
ind_intersection_test = find(s~=0);

% Bail out early if no faces were hit
if isempty(ind_intersection_test)
   return;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plane/ray intersection test %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform plane/ray intersection with the faces that passed 
% the polygon intersection tests. Grab the only the first 
% three vertices since that is all we need to define a plane).
% assuming planar polygons.
candidate_faces = candidate_faces(ind_intersection_test,1:3);
candidate_faces = reshape(candidate_faces',1,numel(candidate_faces));
vert = vert';
candidate_facev = vert(:,candidate_faces);
candidate_facev = reshape(candidate_facev,3,3,length(ind_intersection_test));

% Get three contiguous vertices along polygon 
v1 = squeeze(candidate_facev(:,1,:));
v2 = squeeze(candidate_facev(:,2,:));
v3 = squeeze(candidate_facev(:,3,:));

% Get normal to face plane
vec1 = [v2-v1];
vec2 = [v3-v2];
crs = cross(vec1,vec2);
mag = sqrt(sum(crs.*crs));
nplane(1,:) = crs(1,:)./mag;
nplane(2,:) = crs(2,:)./mag;
nplane(3,:) = crs(3,:)./mag;

% Compute intersection between plane and ray
cp1 = cp(:,1);  
cp2 = cp(:,2);  
d = cp2-cp1;
dp = dot(-nplane,v1);

%A = dot(nplane,d);
A(1,:) = nplane(1,:).*d(1);
A(2,:) = nplane(2,:).*d(2);
A(3,:) = nplane(3,:).*d(3);
A = sum(A,1);

%B = dot(nplane,pt1) 
B(1,:) = nplane(1,:).*cp1(1);
B(2,:) = nplane(2,:).*cp1(2);
B(3,:) = nplane(3,:).*cp1(3);
B = sum(B,1);

% Distance to intersection point
t = (-dp-B)./A;

% Find "best" distance (smallest)
[tbest ind_best] = min(t);

% Determine intersection point
pout = cp1 + tbest .* d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign additional output variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1
    
    % Get face index and vertices
    faceiout = ind_is_face_spanning_x(ind_intersection_test(ind_best));
    facevout = vert(:,faces(faceiout,:));
    
    % Determine index of closest face vertex intersected 
    facexv = xvert(:,faces(faceiout,:));
    dist = sqrt(facexv(1,:).*facexv(1,:) +  facexv(2,:).*facexv(2,:));
    min_dist = min(dist);
    min_index = find(dist==min_dist);
    
    % Get closest vertex index and vertex
    viout = faces(faceiout,min_index);
    vout = vert(:,viout);
end

%--------------------------------------------------------%
function [p] = local_Data2PixelTransform(ax,vert)
% Transform vertices from data space to pixel space.

% Get needed transforms
xform = get(ax,'x_RenderTransform');
offset = get(ax,'x_RenderOffset');
scale = get(ax,'x_RenderScale');

% Equivalent: nvert = vert/scale - offset;
nvert(:,1) = vert(:,1)./scale(1) - offset(1);
nvert(:,2) = vert(:,2)./scale(2) - offset(2);
nvert(:,3) = vert(:,3)./scale(3) - offset(3);

% Equivalent xvert = xform*xvert;
w = xform(4,1) * nvert(:,1) + xform(4,2) * nvert(:,2) + xform(4,3) * nvert(:,3) + xform(4,4);
xvert(:,1) = xform(1,1) * nvert(:,1) + xform(1,2) * nvert(:,2) + xform(1,3) * nvert(:,3) + xform(1,4);
xvert(:,2) = xform(2,1) * nvert(:,1) + xform(2,2) * nvert(:,2) + xform(2,3) * nvert(:,3) + xform(2,4);

% w may be 0 for perspective plots 
ind = find(w==0);
w(ind) = 1; % avoid divide by zero warning
xvert(ind,:) = 0; % set pixel to 0

p(:,1) = xvert(:,1) ./ w;
p(:,2) = xvert(:,2) ./ w;
