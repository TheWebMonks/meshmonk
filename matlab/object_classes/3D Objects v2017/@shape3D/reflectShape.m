function out = reflectShape(obj,ReflectIndex)
         out = clone(obj);
         out.Vertices(:,1) = -1*out.Vertices(:,1);
         out.Vertices = out.Vertices(ReflectIndex,:);
         out.FlipNormals = true;% normals are to be flipped to maintain same inside versus outside standard
end