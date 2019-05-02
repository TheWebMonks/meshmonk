function [m, rawcelloutput] = readTextLandmarkFile(filename,delimeter,startcolumn,endcolumn,startrow,endrow)
%reads a delimited text file (filename) where columns are separated by a demlimiter character (delimiter) and attempts to convert the entries
%between startcolumn and endcolumn and strat row and endrow into a matrix
%of type double
numlines = countLines(filename);
rawcelloutput = cell(numlines,1);
if nargin<6 || isempty(endrow)
    endrow = numlines;
end

if nargin<5 || isempty(startrow)
    startrow = 1;
end

if nargin<4 || isempty(endcolumn)
    endcolumn = NaN;
end

if nargin<3 || isempty(startcolumn)
    startcolumn = 1;
end
    
fid= fopen(filename);

currrow = 0;
while currrow<endrow
    currrow = currrow+1;
    line = fgetl(fid);
    rawcelloutput{currrow} = split(line,delimeter);
end

for r = startrow:endrow
    rowvals = rawcelloutput{r};
    if r==startrow
        if isnan(endcolumn)
            endcolumn = numel(rowvals);
        end
            %initialise matrix
            m = zeros(endrow-(startrow-1),endcolumn-(startcolumn-1));
    end
    for c = startcolumn:endcolumn
        m((r-(startrow-1)),(c-(startcolumn-1))) = str2double(rowvals{c});
    end
end

fclose(fid);

    

end

function out = countLines(filename)
    fid = fopen(filename);
    out = 0;
    line = fgetl(fid);
    while ischar(line)
        out= out+1;
        line = fgetl(fid);
    end
    fclose(fid);
end