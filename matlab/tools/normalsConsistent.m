function [out] = normalsConsistent(Target,Floating)

    n = knnsearch(Target.Vertices,Floating.Vertices,'K',3);
    
    % get normals of three closest points
    correspondingNormals = zeros(size(n,1),3,size(n,2));
    for i = 1:size(n,2)
        correspondingNormals(:,:,i) = Target.VertexNormals(n(:,i),:);                       
    end

    %get average of them
    correspondingNormals = mean(correspondingNormals,3);

    %get cosine between normals
    cosine = sum(Floating.VertexNormals.*correspondingNormals,2);

    %check if mostly positive (consistent) or negative (inconsistent)
    if median(cosine)<0.
        out = false;
    else
        out = true;
    end
end

