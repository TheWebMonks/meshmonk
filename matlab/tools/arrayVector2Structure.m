function out = arrayVector2Structure(in)
            [ncases, dim2] = size(in);
            npoints = dim2/3;
            out = reshape(in',[3,npoints,ncases]);
            out = permute(out, [2,1,3]);
end
 
