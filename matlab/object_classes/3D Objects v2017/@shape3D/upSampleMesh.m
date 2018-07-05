function varargout = upSampleMesh(obj,mode,val,varargin)
    if nargout == 1
           obj = clone(obj);
           obj.Visible = false;
           varargout{1} = obj;
    end
    [~,Findex] = getVindexFindex(obj,varargin{:});
% switching mode    
    switch mode
        case 'runs'
            if (val<1), return; end
            for r=1:1:val
                Findex = performSubdivision(obj,Findex);
            end
        case 'size'
            if val<=0, return; end
            for r=1:1:50% safety maximum number of iterations
                face = obj.Faces(Findex,:)';
                nrF = length(Findex);
                LOC =  zeros(3,nrF,3);AB = zeros(3,nrF);AC = zeros(3,nrF);
                for i=1:1:3
                    LOC(:,:,i) = reshape(obj.Vertices(face(:),i)',3,size(face,2));
                    AB(i,:) = LOC(1,:,i)-LOC(2,:,i);
                    AC(i,:) = LOC(1,:,i)-LOC(3,:,i);
                end
                areas = 0.5*sqrt(dot(AB,AB).*dot(AC,AC)-dot(AB,AC).^2);
                Aindex = find(areas>val);
                if isempty(Aindex), break; end% no triangles to subdivide
                Findex = Findex(Aindex);
                Findex = performSubdivision(obj,Findex);              
            end
            if(r==50), disp('maximum number of iterations');end
        otherwise
            return;
    end
end

function Findex = performSubdivision(obj,Findex)
    if isempty(Findex), return; end
    % one run of subdivision
    face = obj.Faces(Findex,:)';
    keepFindex = setdiff((1:1:obj.nFaces),Findex);
    vertex = obj.Vertices';
    n = size(vertex,2);
    % TRIANGLES
    i = [face(1,:) face(2,:) face(3,:) face(2,:) face(3,:) face(1,:)];
    j = [face(2,:) face(3,:) face(1,:) face(1,:) face(2,:) face(3,:)];
    I = find(i<j);
    i = i(I); j = j(I);
    [~,I] = unique(i + 1234567*j);
    i = i(I); j = j(I);
    ne = length(i); % number of edges
    s = n+(1:ne);
    A = sparse([i;j],[j;i],[s;s],n,n);
    try
        v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
        v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
        v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );
    catch %#ok<CTCH>
        disp('Maximum number of Triangles obtained');
        Findex = [];
        return;
    end  
    face = [    cat(1,face(1,:),v12,v31),...
                cat(1,face(2,:),v23,v12),...
                cat(1,face(3,:),v31,v23),...
                cat(1,v12,v23,v31)   ];
    % VERTICES
    obj.Vertices = ([vertex, (vertex(:,i)+vertex(:,j))/2 ])';
    clear vertex;
    obj.Faces = [obj.Faces(keepFindex,:); face'];
    Findex = (1:size(face,2))+length(keepFindex);
    clear face;
    % VALUE INFORMATION
    if ~isempty(obj.VertexValue)
        val = obj.VertexValue';
        val = [val (val(i)+val(j))/2];
        obj.VertexValue = val';
    end
    % COLOR INFORMATION
    if ~isempty(obj.VertexRGB)
       rgb = obj.VertexRGB';
       rgb = [rgb, (rgb(:,i)+rgb(:,j))/2 ];
       obj.TextureColor = rgb';
       clear rgb;
    end
    if ~isempty(obj.UV)
        uv = obj.UV';
        uv = [uv (uv(:,i)+uv(:,j))/2 ];
        obj.UV = uv';
        clear uv;
    end
end