function out = readLandmarksMeVisLabXML(fn)
    f = fopen(fn);
    l = fgetl(f);
    landmarkCount = 0;
    while ischar(l)
        checkListSize = strfind(l,'<ListSize>');
        if ~isempty(checkListSize)
           chunks = strsplit(l,'>' );
           chunks =strsplit(chunks{2},'<');
           nLandmarks = str2num(chunks{1});
           out = zeros(nLandmarks,3);
        end
        
        
        
        check = strfind(l,'<pos>');
        if ~isempty(check)
            landmarkCount = landmarkCount+1
            chunks = strsplit(l((check+5):end),' ');
            
            for i = 1:3
                out(landmarkCount,:) = str2num(chunks{i});
            end
        end
            
        l = fgetl(f);
            
        
    end
    fclose(f);
end
    
    




