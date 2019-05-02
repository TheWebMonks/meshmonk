function out = applyTransform(obj,T)
    if isobject(obj)
        vertices = obj.Vertices;
    else
        vertices=obj;
    end
    
    if size(vertices,1)~=3
        vertices = vertices';
    end
    
    assert(size(vertices,1)==3);
    nverts = size(vertices,2);

    
    tVertices = (T.Scale*(T.Rotation*(vertices)+repmat(T.Translation,1,nverts)))';
    
    if isobject(obj)
        out = clone(obj);
        out.Vertices = tVertices;
    else
        out = tVertices;
    end
    
    
end