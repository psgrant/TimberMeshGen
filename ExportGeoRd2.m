%% So we're gonna like hopefully get the .geo file exported
% load 'RingData.mat'
xc = M;
yc = N;
clc

lc = 50;
% define boundary bpoins
CornerPoints = [0 0;
    xc 0;
    xc yc;
    0 yc];
[BoundaryPoints] = FindBoundaryPoints(MinEndPoints,MaxEndPoints,LateEndPoints,N,M);

% BoundaryList = unique(BoundaryList,'rows');
% Number of nodes per length
NodeDensity = 20;
NumRings = length(IsolatedRingsIDs);
CornerIds = [1 0 0 0];
% Gets the number of polynomial points on the edge including the 4 corners
NumBoundaryPoints = length([BoundaryPoints{1};BoundaryPoints{2};BoundaryPoints{3};BoundaryPoints{4}])+4;
BoundaryList = zeros(NumBoundaryPoints,2);
BoundaryList(1,:) = CornerPoints(1,:);
JumpDist = 2;
% Loop over each edge
for i = ((1:4))
    NumEdgePoints = size(BoundaryPoints{i},1);
    if i < 4
        BoundaryList((JumpDist + NumEdgePoints),:) = CornerPoints(i+1,:);
        CornerIds(i+1) = JumpDist + NumEdgePoints;
    end
    BoundaryList(JumpDist:(JumpDist + NumEdgePoints-1),:) = BoundaryPoints{i};

    JumpDist = JumpDist + NumEdgePoints+1;

end
LastPoint = length(BoundaryList) + 1;
SurfCornerIds = [CornerIds, LastPoint];

% Move nodes that are close to the corners onto the corners.
% 2.5% of the mean of both dimensions
CutOff = 0.025*0.5*(M+N);
OldBoundaryList = BoundaryList;
% loop over each corner
for i = 1:4
    % Find the boundary points within the cutoff distance
    idx = find(sqrt(sum((BoundaryList - CornerPoints(i,:)).^2,2)) < CutOff)
    % remove corner point from selection
    idx(idx == SurfCornerIds(i)) = [];
    % Assign the value of that point to the corner
    if ~isempty(idx)
        BoundaryList(idx,:) = repmat(CornerPoints(i,:),[length(idx),1]);
    end
end


[~,a] = unique(BoundaryList,'rows');
% BoundaryList = BoundaryList(sort(a),:);
NumBoundaryPoints = size(BoundaryList,1);
% hold on
% plot(BoundaryList(:,1),BoundaryList(:,2),'go','markersize',20,'linewidth',4)


% Start generating the string for the points


% Print geo string%
% Initalise boudnary point array
geoStrPoints = cell(NumBoundaryPoints,1);
StartPoint = 0;
StartLine = 0;
PointId = StartPoint + 1;
LineId = StartLine + 1;
% Loop through each Point
for pointi = ((1:NumBoundaryPoints))
    geoStrPoints{pointi} = sprintf('Point(%d) = {%f,%f,0,%g};\r\n',...
        PointId,BoundaryList(pointi,1),BoundaryList(pointi,2),lc);
    PointId = PointId +1;
end
geostrBP = strjoin(geoStrPoints);

% Connect the boundary points together together

LinePointsArray = ([(1:NumBoundaryPoints-1); (2:NumBoundaryPoints)])';
LinePointsArray(end+1,:) = [NumBoundaryPoints, 1];


% Initalise boundary line array
GeoBoundaryLines = cell(NumBoundaryPoints,1);

% Loop through each Point
for Line = ((1:NumBoundaryPoints))
    %     Calculates the distance between each pair of boundary points
    Dist = pdist(BoundaryList(LinePointsArray(Line,:),:));
    NodesOnLine = max(round( Dist / BoxDensity(1) ),3);
    % linestr = 'Line(%d) = {%d,%d};\r\nTransfinite Line{%d} = %d;\r\n';
    linestr = 'Line(%d) = {%d,%d};\r\n';
    GeoBoundaryLines{Line} = sprintf(linestr,LineId, ...
        LinePointsArray(Line,1)+StartPoint,LinePointsArray(Line,2)+StartPoint);
    %         GeoBoundaryLines{Line} = sprintf('Line(%d) = {%d,%d};\r\n',...
    %         Line,LinePointsArray(Line,1),LinePointsArray(Line,2));
    LineId = LineId + 1;
end

geostrBL = strjoin(GeoBoundaryLines);


% Make the curve loop and create the surface objects
CurveCell = cell(3,1);
CurveCell{1} = 'Curve Loop(1) = {';
CurveCell{2} = sprintf('%d:%d,',StartLine+1,NumBoundaryPoints + StartLine);
CurveCell{2}(end) = [];
CurveCell{3} = sprintf('};\r\n Plane Surface(1) = {1};\r\n');
geostrSurf = strjoin(CurveCell);



% Start with embedding the middle points into the surfaces

% get the x and y positions of the middle points for both rings
% MiddlePoints = cell(2,1);
% MiddlePoints{1}  = zeros(NumRings,2);
% MiddlePoints{2}  = zeros(NumRings,2);
% for i = 1:NumRings
%    Data =  OutputSpline{i}{1};
%    MiddlePoints{1}(i,:) = Data(2,:);
%    Data =  OutputSpline{i}{2};
%    MiddlePoints{2}(i,:) = Data(2,:);
%
% end

%  We will lopp over each polynomial and find the point Ids for the
%  boundary points, then we will embed the middle point (noting the point
%  ID). Then we will need to create a line through the points (transfinite
%  node spacing based on line length).


PointCount = 1;
geostrRingLines = cell(2*NumRings,1);
RemoveIndexes = [];
% Loop over each ring
for i = 1:NumRings

    % For min and max lines
    StartPointId = PointId;
    for k = kout
        Index = (k-1)*NumRings + i;
        P = SplineCell{i,k};
        if isempty(P)
            RemoveIndexes = [RemoveIndexes Index]; %#ok
            continue
        end
        StartPoint = SplinePoints{i,k}(1,1);
        EndPoint = SplinePoints{i,k}(2,1);
        if StartPoint > EndPoint
            x = EndPoint:StartPoint;
        else
            x = StartPoint:EndPoint;
        end
        
        y = P(x);
        %         plot(x,y,'kx','markersize',15,'linewidth',2)
        
        % If the ring is in a corner we skip over generating lines for ut

        InterArray = PolyIntersections{Index};


        % Isolate x and y positions for start point and end point
        P1 = [x(1), y(1)];
        P3 = [x(end), y(end)];

        % Get the boundary point for the 1st and 3rd boundary point
        P1Num = find(sum(abs(OldBoundaryList - P1) < 1,2)==2)+0;%StartPointId;
        P3Num = find(sum(abs(OldBoundaryList - P3) < 1,2)==2)+0;%StartPointId;

        WholeRingStr = cell(size(InterArray,1),1);


        for section = 1:size(InterArray,1)
            SecLength = pdist([InterArray(section,2:3);InterArray(section,4:5)]);
            SplineInterpPoints = max(round(SecLength/50),3);
            XSection = round(linspace(InterArray(section,2),InterArray(section,4),SplineInterpPoints));
            YSection = P(XSection);

            % The end point of the xection
            P2 = [XSection(end), YSection(end)];

            % For each interal point (not the two boundary points) put them
            % in the mesh geomtry
            PrevPointId = PointId - 1;
            PointArray = [];
            PointStops = length(XSection);

            if size(InterArray,1) == 1
                PointStops = PointStops - 1;
            end
            PointStr = cell(PointStops,1);
            for point = 2:PointStops
                % Save as string
                PointStr{point-1} = sprintf('Point(%d) = {%f,%f,0,%g};\r\n',...
                    PointId,XSection(point),YSection(point),BoxDensity(InterArray(section,1)));
                PointArray = [PointArray, PointId]; %#ok
                PointId = PointId  + 1;
                PointCount = PointCount + 1;
            end
            % Point Id for the point at the end of the spline
            % Put spline through

            % if the ring only goes through one section

            if size(InterArray,1) == 1
                TempStr = [sprintf('%.0f,',P1Num), sprintf('%.0f,',PointArray),sprintf('%.0f,',P3Num)];

            else

                switch section % Differnt spine ording for different sections

                    case 1 % First section need to get the point id of the 1st boundary point
                        TempStr = [sprintf('%.0f,',P1Num), sprintf('%.0f,',PointArray)];

                    case size(InterArray,1) % The end section
                        TempStr = [sprintf('%.0f,',PrevPointId), sprintf('%.0f,',PointArray),sprintf('%.0f,',P3Num)];

                    otherwise % Middle sections
                        TempStr = [sprintf('%.0f,',PrevPointId), sprintf('%.0f,',PointArray)];

                end
            end

            % Compute the number of nodes for each section

            NodesOnSection = max(round(SecLength / BoxDensity(InterArray(section,1))),3);
            if k == 3
                NodesOnSection = round(NodesOnSection*0.9);
            end
            % Remove extra comma
            TempStr = TempStr(1:end-1);
            % Export String
            PointStr{point} = [sprintf('Spline(%.0f) = {',LineId),TempStr,...
                sprintf('};\r\nTransfinite Curve {%.0f} = %g;\r\nCurve{%d} In Surface{1};\r\n'...
                ,LineId,NodesOnSection,LineId)]; % Add in node nums
                % sprintf('};\r\nCurve{%d} In Surface{1};\r\n',LineId)];
                
            LineId = LineId + 1;


            WholeRingStr{section} = strjoin(PointStr);
        end


        geostrRingLines{Index} = strjoin(WholeRingStr);
    end
end

SurfStart = SurfCornerIds(1:end-1) + 2;
SurfEnd = SurfCornerIds(2:end) + 1;
SurfTagCell = cell(6,1);
j = 0;
for i = 11:16
    j = j + 1;
    switch j
        case 1
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Bottom", %g) = {%g:%g};\r\n',i,SurfStart(j),SurfEnd(j));
        case 2
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Left", %g) = {%g:%g};\r\n',i,SurfStart(j),SurfEnd(j));
        case 3
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Top", %g) = {%g:%g};\r\n',i,SurfStart(j),SurfEnd(j));
        case 4
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Right", %g) = {%g:%g};\r\n',i,SurfStart(j),SurfEnd(j));
        case 5
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Front", %g) = {%g};\r\n',i,1);
        case 6
            SurfTagCell{j} = sprintf('Physical Surface("Mesh Back", %g) = {%g};\r\n',i,SurfEnd(end)+1);
        otherwise
    end
end

SurftagStr = strjoin(SurfTagCell);


geostrRingLines(RemoveIndexes) = [];
B = [];
for idx = 1:numel(geostrRingLines)
    if ~isempty(geostrRingLines{idx})
        B = [B, idx];
    end
end
A = geostrRingLines(B);
% for i = 1:5
% A{i} = geostrRingLines{i}
% end
geostrEmbSplines = strjoin(A);
% geostrEmbSplines = (geostrRingLines{i+12});
ExtrudeStr = 'out[] = Extrude {0, 0, 2702} {Surface{1}; Layers{ {2,3,3,2}, {0.25,0.5,0.75,1}}; Recombine;};\r\nPhysical Volume("Mesh Volume",1)={out[1]};\r\nPhysical Volume("Mesh Surface",1)={out[0]};\r\nPhysical Surface("Mesh Surface", 2) = {1:1000};\r\nMesh 3;';
% ExtrudeStr = '';
GeomStr = 'Geometry.OldNewReg=0;\r\n';
totGeoStr = cat(2,GeomStr,geostrBP,geostrBL,geostrSurf,geostrEmbSplines,ExtrudeStr,SurftagStr);


% 



fileID = fopen('TestMesh.geo','w');
fprintf(fileID, totGeoStr);
fclose(fileID);


