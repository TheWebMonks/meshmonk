% function importWavefront(obj,filename,path,usePrefix)
%                if nargin<4, usePrefix= []; end
%                if nargin<3, path = pwd; end
%                if ~strcmp(path(end),'/'), path = [path '/'];end
%                if ~strcmp(filename(end-3:end),'.obj'), filename = [filename '.obj'];end
%                obj.Tag = filename(1:end-4);
%                filename = [path filename];
%              % quick inreading of vertices, normals and texturecoordinates;
%              try
%                fid = fopen(filename,'r');
%                names = textscan(fid,'%s%*[^\n]');names = names{1};
%                fclose(fid);
%                fid = fopen(filename,'r');
%                data = textscan(fid,'%[^\n\r]');data = data{1};
%                fclose(fid);
%              catch
%                  error('File does not exist, check path and filename');
%              end
%              % reading Vertices
%                list = strcmp('v',names);index = find(list==1);
%                if isempty(index), return; end %nothing in the file
%                % small test to preallocate memory
%                Lyn = data{index(1)};tmp = sscanf(Lyn(2:end),'%f');
%                Location = nan*zeros(length(index),size(tmp,1));
%                fordata = data(index);
%                parfor i=1:1:length(index)
%                    Lyn = fordata{i};Location(i,:) = sscanf(Lyn(2:end),'%f');
%                end
%                obj.Vertices = Location(:,1:3);
%                if size(Location,2)>3,  obj.VertexRGB = Location(:,4:end)./255; else, obj.VertexRGB = []; end 
%                clear list index Location;
%              % reading normals 
%                % not required because computed in shape3D class on the fly    
%              % reading texture coordinates
%                list = strcmp('vt',names);index = find(list==1);
%                if ~isempty(index)
%                   fordata = data(index);
%                   uv = nan*zeros(length(index),2);
%                   parfor i=1:1:length(index)
%                       Lyn = fordata{i};uv(i,:) = sscanf(Lyn(3:end),'%f');
%                   end
%                else
%                   uv = [];
%                end  
%                clear list index;  
%              % reading faces    
%                list = strcmp('f',names);index = find(list==1);
%                if ~isempty(index)
%                   F3 = nan*zeros(length(index),3);F3t = F3;
%                   %simply post a warning when Quads are presented
%                   fordata = data(index);
%                   parfor i=1:1:length(index)
%                       Lyn = fordata{i};
%                       str = textscan(Lyn(2:end),'%s'); str = str{1};
%                       nf = length(strfind(str{1},'/'));
%                       tEmpty = length(strfind(str{1},'//'));
%                       [tok,str] = strtok(str,'//');     %#ok<*STTOK> % add vertex
%                       F3(i,:) = str2double(tok)';
%                       if (nf > 0) && tEmpty==0
%                           [tok,~] = strtok(str,'//');   % add texture coordinates
%                           F3t(i,:) = str2double(tok)';
%                       end
%                    end
%                end
%                Findex = index;
%                obj.Faces = F3;
%              % TEXTURE
%                if isempty(uv), return; end
%              % reading in material lib file, and texturemaps
%                try
%                     [usemtl,mtl,maps] = parseMTL(filename,usePrefix);
%                catch
%                   tmp = dir([filename(1:end-4) '.bmp']);
%                   if isempty(tmp)
%                      usemtl = false;
%                   else
%                      usemtl = true;
%                      mtl{1} = 'map';
%                      maps{1} = imread([filename(1:end-4) '.bmp']);
%                   end
%                end
%                if ~usemtl, warning('TextureMaps are being ignored, because .mtl file is missing'); return; end
%                switch length(mtl)>1
%                    case false% SINGLE IMAGE
%                        map = maps{1};
%                    case true % MUTLIPLE IMAGES/FOR NOW I EXPECT THREE SEPERATE ONES FOR FACES!
%                        % look for usemtl definitions in .obj file
%                        nmtl = length(mtl);
%                        usemtlind = +inf*ones(1,nmtl);      
%                        list = strcmp('usemtl',names);index = find(list==1);
%                        listdata = data(index); 
%                        for i=1:1:nmtl
%                           for j=1:1:length(listdata)
%                               if contains(listdata{j},mtl{i})
%                                  usemtlind(i) = index(j);
%                                  break;
%                               end
%                           end
%                         end
%                        % convert to a single texture map, for three image exports
%                         [map,convert] = convertVectra2SingleTexture(maps);
%                         Fmtl = Findex;
%                         counter = 0;
%                         for i=1:1:length(usemtlind)
%                             Fmtl(Findex>=usemtlind(i)) = counter;
%                             counter = counter+1;
%                         end
%                         convert = convert(2:end);
% 
%                         %newuv = uv;
%                         for i=1:1:3
%                            forF3t = F3t(Fmtl==i,:);
%                            index = unique(forF3t);
%                            foruv = uv(index,:);
%                            forconvert = convert{i};
% 
%                            foruv(:,1) = rescaleUV(foruv(:,1),forconvert(1,1),forconvert(1,2));
%                            foruv(:,2) = rescaleUV(foruv(:,2),forconvert(2,1),forconvert(2,2));
%                            %newuv(index,:) = foruv;
%                            uv(index,:) = foruv;
%                         end
%                         %uv = newuv;
%                end
%               % readjusting UV to Vertex order 
%                 [~,I,~] = unique(F3(:),'last');
%                 Vert_Index = F3(I);% nice 1 to last
%                 crop(obj,'VertexIndex',Vert_Index,'Action','crop');% adjust vertices and triangles (remove those vertices not referenced in any triangles)
%            
%                 UV_Index = F3t(I);
%                 UV = nan*zeros(length(Vert_Index),2);
%                 good = find(~isnan(UV_Index));
%                 UV(good,:) = uv(UV_Index(good),:);
%                 UV(isnan(UV_Index),:) = 0;% bad references
%                 UV(:,2) = 1-UV(:,2);% flip image height reference
%                 obj.UV = UV;
%                 obj.TextureMap = map;
%                 obj.ColorMode = 'Texture';% set texture on when rendering
%                % MANIFOLD CLEANING
%                 cleanManifold(obj);
%                                    
% end
% % 
% %                 [w,h,~] = size(map);
% %                 figure;image(map);
% %                 hold on;
% %                 plot(obj.UV(:,1)*h,obj.UV(:,2)*w,'b.');  
% 
%         

function importWavefront(obj,filename,path,usePrefix)
               if nargin<4, usePrefix= []; end
               if nargin<3, path = pwd; end
               if ~strcmp(path(end),'/'), path = [path '/'];end
               if ~strcmp(filename(end-3:end),'.obj'), filename = [filename '.obj'];end
               obj.Tag = filename(1:end-4);
               filename = [path filename];
             % quick inreading of vertices, normals and texturecoordinates;
             try
               fid = fopen(filename,'r');
               names = textscan(fid,'%s%*[^\n]');names = names{1};
               fclose(fid);
               fid = fopen(filename,'r');
               data = textscan(fid,'%[^\n\r]');data = data{1};
               fclose(fid);
             catch
                 error('File does not exist, check path and filename');
             end
             % reading Vertices
               list = strcmp('v',names);index = find(list==1);
               if isempty(index), return; end %nothing in the file
               % small test to preallocate memory
               Lyn = data{index(1)};tmp = sscanf(Lyn(2:end),'%f');
               Location = nan*zeros(length(index),size(tmp,1));
               fordata = data(index);
               for i=1:1:length(index)
                   Lyn = fordata{i};Location(i,:) = sscanf(Lyn(2:end),'%f');
               end
               obj.Vertices = Location(:,1:3);
               if size(Location,2)>3,  obj.VertexRGB = Location(:,4:end)./255; else, obj.VertexRGB = []; end 
               clear list index Location;
             % reading normals 
               % not required because computed in shape3D class on the fly    
             % reading texture coordinates
               list = strcmp('vt',names);index = find(list==1);
               if ~isempty(index)
                  fordata = data(index);
                  uv = nan*zeros(length(index),2);
                  for i=1:1:length(index)
                      Lyn = fordata{i};uv(i,:) = sscanf(Lyn(3:end),'%f');
                  end
               else
                  uv = [];
               end  
               clear list index;  
             % reading faces    
               list = strcmp('f',names);index = find(list==1);
               if ~isempty(index)
                  F3 = nan*zeros(length(index),3);F3t = F3;
                  %simply post a warning when Quads are presented
                  fordata = data(index);
                  for i=1:1:length(index)
                      Lyn = fordata{i};
                      str = textscan(Lyn(2:end),'%s'); str = str{1};
                      nf = length(strfind(str{1},'/'));
                      tEmpty = length(strfind(str{1},'//'));
                      [tok,str] = strtok(str,'//');     %#ok<*STTOK> % add vertex
                      F3(i,:) = str2double(tok)';
                      if (nf > 0) && tEmpty==0
                          [tok,~] = strtok(str,'//');   % add texture coordinates
                          F3t(i,:) = str2double(tok)';
                      end
                   end
               end
               Findex = index;
               obj.Faces = F3;
             % TEXTURE
               if isempty(uv), return; end
             % reading in material lib file, and texturemaps
               try
                    [usemtl,mtl,maps] = parseMTL(filename,usePrefix);
               catch
                  tmp = dir([filename(1:end-4) '.bmp']);
                  if isempty(tmp)
                     usemtl = false;
                  else
                     usemtl = true;
                     mtl{1} = 'map';
                     maps{1} = imread([filename(1:end-4) '.bmp']);
                  end
               end
               if ~usemtl, warning('TextureMaps are being ignored, because .mtl file is missing'); return; end
               switch length(mtl)>1
                   case false% SINGLE IMAGE
                       map = maps{1};
                   case true % MUTLIPLE IMAGES/FOR NOW I EXPECT THREE SEPERATE ONES FOR FACES!
                       % look for usemtl definitions in .obj file
                       nmtl = length(mtl);
                       usemtlind = +inf*ones(1,nmtl);      
                       list = strcmp('usemtl',names);index = find(list==1);
                       listdata = data(index); 
                       for i=1:1:nmtl
                          for j=1:1:length(listdata)
                              if contains(listdata{j},mtl{i})
                                 usemtlind(i) = index(j);
                                 break;
                              end
                          end
                        end
                       % convert to a single texture map, for three image exports
                        [map,convert] = convertVectra2SingleTexture(maps);
                        Fmtl = Findex;
                        counter = 0;
                        for i=1:1:length(usemtlind)
                            Fmtl(Findex>=usemtlind(i)) = counter;
                            counter = counter+1;
                        end
                        convert = convert(2:end);

                        %newuv = uv;
                        for i=1:1:3
                           forF3t = F3t(Fmtl==i,:);
                           index = unique(forF3t);
                           foruv = uv(index,:);
                           forconvert = convert{i};

                           foruv(:,1) = rescaleUV(foruv(:,1),forconvert(1,1),forconvert(1,2));
                           foruv(:,2) = rescaleUV(foruv(:,2),forconvert(2,1),forconvert(2,2));
                           %newuv(index,:) = foruv;
                           uv(index,:) = foruv;
                        end
                        %uv = newuv;
               end
              % readjusting UV to Vertex order 
                [~,I,~] = unique(F3(:),'last');
                Vert_Index = F3(I);% nice 1 to last
                crop(obj,'VertexIndex',Vert_Index,'Action','crop');% adjust vertices and triangles (remove those vertices not referenced in any triangles)
           
                UV_Index = F3t(I);
                UV = nan*zeros(length(Vert_Index),2);
                good = find(~isnan(UV_Index));
                UV(good,:) = uv(UV_Index(good),:);
                UV(isnan(UV_Index),:) = 0;% bad references
                UV(:,2) = 1-UV(:,2);% flip image height reference
                obj.UV = UV;
                obj.TextureMap = map;
                obj.ColorMode = 'Texture';% set texture on when rendering
               % MANIFOLD CLEANING
                cleanManifold(obj);
                                   
end

        

