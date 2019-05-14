function out = arrayStructure2Vector(in)
                [x,y,z] = size(in);
                % check second dimension is length 3
                if y~=3
                    error('Vertices must be shaped [no vertices, 3, no shapes] or a 2D matrix');
                end 
                
                %reshape array - there must be a smpler way to do this
                perm = permute(in,[2,1,3]); out = reshape(perm,[y*x,z])';
                
 end
  
