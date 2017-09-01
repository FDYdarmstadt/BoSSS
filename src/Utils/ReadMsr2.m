function Mtx = ReadMsr2(filename)

fid = fopen(filename);

NoOfRows = fscanf(fid,'%d',1);
NoOfCols = fscanf(fid,'%d',1);
NonZeros = fscanf(fid,'%d',1);


%Mtx = spalloc(NoOfRows,NoOfCols,NonZeros);
cnt = 1;

iCol = zeros(NonZeros,1);
iRow = zeros(NonZeros,1);
entries = zeros(NonZeros,1);
l0 = 0;
str = char(zeros(1,6));

for i = 1:NoOfRows
    NonZerosInRow = fscanf(fid,'%d',1);
    
    if(l0 ~= NonZerosInRow)
        str = char(zeros(1,NonZerosInRow*6));
        for j = 1:NonZerosInRow
            i0 = 1+(j-1)*6;
            str(i0:i0+1) = '%f';
            str(i0+3:i0+4) = '%f';
        end
    end
    
    R = fscanf(fid,str,2*NonZerosInRow);
    R2 = reshape(R',2,NonZerosInRow);
    
    ind = cnt:(cnt+NonZerosInRow-1);
    iCol(ind) = R2(1,:);
    iRow(ind) = i;
    entries(ind) = R2(2,:);
    
    
    
%         %entry = fscanf(fid,'%f',1);
    
    
%     for j = 1:NonZerosInRow
%         
%         %iCol  = fscanf(fid,'%d',1);
%         %entry = fscanf(fid,'%f',1);
%         iCol=1;
%         entry=4;
%         
%         Mtx(i,iCol+1) = entry;
%         
%         
%         if(cnt >= NonZeros)
%             error(0,'Number of nonzeros is wrong');
%         end
%         
%         cnt = cnt + 1;
%     end
    
    cnt = cnt + NonZerosInRow;
    i;
end

fclose(fid);


Mtx = sparse(iRow,iCol+1,entries,NoOfRows,NoOfCols,NonZeros);
