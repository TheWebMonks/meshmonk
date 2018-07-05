function out = reMesh(obj,dist,res)
         if nargin<3, res = 2; end
         if nargin<2, dist = res*3;end
         
          res = 1.5;
          dist = res*3;
         
         % create voxel grid
          range = [min(obj.Vertices);max(obj.Vertices)];
          [X,Y,Z] = meshgrid(range(1,1)-res:res:range(2,1)+res,range(1,2)-res:res:range(2,2)+res,range(1,3)-res:res:range(2,3)+res);
          PL = [X(:),Y(:),Z(:)];
          n = size(PL,1);
         % fast KNN approximiation           
          [~,D] = knnsearch(obj.Vertices,PL);
         % find closeby points for an exact computation of the distance
          index = find(D<=dist);
          remindex = setdiff(1:n,index);
          subD = point2trimesh('Faces',obj.Faces,'Vertices',obj.Vertices,'QueryPoints',PL(index,:),'MaxDistance',1,'Algorithm','parallel');
         % generating the whole map into a signed distance transform 
          D(index) = subD;
          ind = knnsearch(PL(index,:),PL(remindex,:));
          signs = ones(1,n);
          signs(remindex) = sign(subD(ind));
          D = signs'.*D;
         % filling in the grid 
          C = 0*X;
          C(1:n) = D;
         % running marching cubes 
          out = shape3D;
          [out.Faces,out.Vertices] = MarchingCubes(single(X),single(Y),single(Z),single(C),0);
         % cutting of border extrapolations 
          bindex = getBoundary(obj);
          ind = knnsearch(obj.Vertices,out.Vertices);
          goodindex = find(~ismember(ind,bindex));
          crop(out,'VertexIndex',goodindex);
          %viewer(out);
          
          
          test = meshObj;
          test.Faces = out.Faces';
          test.Vertices = out.Vertices';
          
          new = subdivideTriangles(test,'runs',1);
          
          
end