function [AvgRingInter] = CalculateInterceptsMean(T,ParentNode,NumCols,NumRows,XorY)




xy = T.get(ParentNode);
x = xy(:,1);
y = xy(:,2);

if XorY == 1
    ymin = min(y);
    ymax = max(y);
    RingIntercepts = zeros(ymax - ymin,1);
    for k = ymin:ymax
        Slice = zeros(NumCols+2,1);
        [I,~] = find((y==k));
        Slice(x(I)+1) = 1;
        RingIntercepts(k - ymin + 1) = sum(islocalmax(Slice, 'FlatSelection','center'));
    end
    
else % vertical slices
    xmin = min(x);
    xmax = max(x);
    RingIntercepts = zeros(xmax - xmin,1);
    for k = xmin:xmax
        Slice = zeros(NumRows+2,1);
        [I,~] = find((x==k));
        Slice(y(I)+1) = 1;
        RingIntercepts(k - xmin + 1) = sum(islocalmax(Slice, 'FlatSelection','center'));
    end
    
    
    
    
end
AvgRingInter = mean(RingIntercepts);
