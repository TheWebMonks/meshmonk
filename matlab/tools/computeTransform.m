function T = computeTransform(floating, target,scale,w)
                
             p = floating;
             q = target;

             if isobject(q)
                q = q.Vertices';
             end
             if isobject(p)
                 p = p.Vertices';
             end
             
             
            if size(p,1)~=3
                p = p';
            end
            if size(q,1)~=3
                q = q';
            end
            assert(size(q,1)==3);
            assert(size(p,1)==3);
                

             nbpts = size(p,2);
             
             if nargin<4
                 w = ones(1,nbpts);
             end
             
             if nargin<3
                 scale = true;
             end
             
             
             index = find(w);% find points not having weights == 0;
             if length(index)<nbpts
                 p = p(:,index);
                 q = q(:,index);
                 w = w(:,index);
             end
             totalw = sum(w,2);
             rw = repmat(w,3,1);
             nbpts = size(p,2);
             % align centers
             centerq = sum(rw.*q,2)./repmat(totalw,3,1);
             centerp = sum(rw.*p,2)./repmat(totalw,3,1);
             q = q-repmat(centerq,1,nbpts);
             p = p-repmat(centerp,1,nbpts);
             % Scale both by making mean length 1
             if scale
                 lengths = sqrt(sum((rw.*p).^2,1));
                 meanLength1 = sum(lengths,2)/totalw;
                 p = p./meanLength1;
                 lengths = sqrt(sum((rw.*q).^2,1));
                 meanLength2 = sum(lengths,2)/totalw;
                 q = q./meanLength2;
                 scalefactor = meanLength2/meanLength1;
             else
                 scalefactor = 1;% keep scale fixed
             end
             % Rotate
             [U,S,V] = svd(rw.*q*p');
             H = V*sign(S)*U';
             H = H';
             % putting it all together
             T.Scale=scalefactor;
             T.Rotation = H;
             transform = T.Scale*H;
             Ta=eye(4);
             Ta(1:3,4)=centerq;
             Tb=eye(4);
             Tb(1:3,4)=-centerp;
             R=eye(4);
             R(1:3,1:3)=transform;
             Tout=Ta*R*Tb;
             T.Translation = Tout(1:3,4)/T.Scale;
end

