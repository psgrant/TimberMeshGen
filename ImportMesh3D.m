tic
% clear all


% CubeTest21
% TriangleElementsRefined
% BoardTest73K
% SingleBoardTest3
% BoardUniform
% SmallBoardNonUniform
% DoubleLayerBoard
% TimberBoard16k
% SquareVisGov3D
% GovMesh3D1
% CLTTest
% BigCLTBot
% BigCLTBot2
% BigCLTMid
% BigCLTMid2
% TEST13
% CLTMIDTEST1
% TEST15Mid
% DiffusionTestMesh
BCTestMESH
% PaperMeshTop
% PaperMeshBot
% PaperMeshMid


% Papermesh
% UniformCube
% TEST151mid
% TEST151midHiRes
% TEST15Bot
% TEST15Top
% KnotTest18k
% BigCLTTop
% BigCLTTop2
% TestBoard2k
% Cube
% CLTTest2
% Radiata1Coarse
% UniformBlock
% BigCLTTop
% BigCLTMid
% BigCLTBot
% Load in all the mesh data
NodePos = msh.POS;
NodePos1 = NodePos;

% NodePos(NodePos1(:,2) > 2771,:) = NodePos(NodePos1(:,2) > 2771,:) - 2772;
% NodePos(NodePos1(:,2) < 1386,:) = NodePos(NodePos1(:,2) < 1386,:) - 2772;
MIN = msh.MIN;
MAX = msh.MAX;
NodePos(:,1) = NodePos(:,1) .* (0.09/MAX(1));
NodePos(:,2) = NodePos(:,2) .* (0.045/MAX(2));
NodePos(:,3) = NodePos(:,3) .* (0.09/MAX(3));

msh.ScalingFacts(3) = 0.09/MAX(3);
msh.ScalingFacts(2) = 0.045/MAX(2);
msh.ScalingFacts(1) = 0.09/MAX(1);


MIN = msh.MIN .* [(0.09/MAX(1)), (0.045/MAX(2)), (0.09/MAX(3))];
MAX = msh.MAX .* [(0.09/MAX(1)), (0.045/MAX(2)), (0.09/MAX(3))];

Elements = msh.PRISMS(:,1:6);
% Elements = msh.TRIANGLES(:,1:3);
NumNodes = msh.nbNod;
NumElements = size(Elements,1);
Triangles = msh.TRIANGLES(:,1:3);
Quads = msh.QUADS(:,1:4);

SurfTriangles = msh.TRIANGLES;
SurfQuads = msh.QUADS;
% BoundaryNodeIDs = unique([Triangles(:); Quads(:)]);
BoundaryNodeIDs = [];

HalfTriElements = msh.TRIANGLES(1:(size(msh.TRIANGLES,1)/2),1:3);
for i = 1:3


    idxMin = find(NodePos(:,i) == MIN(i));
    idxMax = find(NodePos(:,i) == MAX(i));

    BoundaryNodeIDs = [BoundaryNodeIDs; idxMin; idxMax];

end

BoundaryNodeIDs = unique(BoundaryNodeIDs);

%% Node associtivity
% Return an array where each row and column respond to a node-node pair. A
% 1 indicates thata node is connected (in a control volume way). A zero
% indicates no connection
% NodeMatObj = matfile('NodeConnections.mat');
% NodeMatObj.Properties.Writable = true;
% NodeMatObj.NodeConnections(NumNodes,NumNodes) = 0;
%%


% NodeConnections = sparse(zeros(2));
NodeConnections = spalloc(NumNodes,NumNodes,10*NumNodes);
% NodeConnections(NumNodes,NumNodes) = 0;
LongConnections = zeros(NumNodes,2);
% Loop over each element
for element = 1:NumElements

    ElementNodes = Elements(element,:);
    % For each node in the element find it's connections and set it to 1 in
    % the connectivity matrix.
    for node = 1:3

        % Connect the nodes on the triangular face of the 'front' of the
        % element
        NodeConnections(ElementNodes(node),ElementNodes(1:3)) = 1;


        % Connect the nodes that are directly infront of a node, this
        % ensures that we don't have any nodes being connected 'diagonally'
        % across the element
        NodeConnections(ElementNodes(node),ElementNodes(node+3)) = 1;
        NodeConnections(ElementNodes(node+3),ElementNodes(node)) = 1;

        % Update vector describing the longitudinal connections
        LongConnections(ElementNodes(node),1) = 1;
        LongConnections(ElementNodes(node+3),2) = 1;

    end

    for node = 4:6
        % Connect nodes along the backside of the element
        NodeConnections(ElementNodes(node),ElementNodes(4:6)) = 1;
    end


end



%%
TriTemp = HalfTriElements;
TriNodes = unique(TriTemp);
X = 1:length(TriNodes);
XNodes = NodePos(TriNodes,:);
TriEleTemp = HalfTriElements;
for i = 1:length(TriNodes)

    TriEleTemp(TriEleTemp == TriNodes(i)) = i;

end
msh.TriElements = TriEleTemp;
msh.TriNodes = XNodes;


%% Compute RCM ordering and adjust variables
%
%compute RCM ordering
RCM = symrcm(NodeConnections);

% Reoder Nodes
NodePos = NodePos(RCM,:);
msh.NodePos = NodePos;
BoundaryNodeTemp = BoundaryNodeIDs;
% Reorder Elements
ElementsTemp = Elements;
TriTemp = HalfTriElements;

% For surface quads and tris
TrianglesTemp = SurfTriangles(:,1:3);
QuadsTemp = SurfQuads(:,1:4);

for i = 1:NumNodes

    Index = (Elements == RCM(i));
    ElementsTemp(Index) = i;

    Index = (HalfTriElements == RCM(i));
    TriTemp(Index) = i;

    Index = (Triangles == RCM(i));
    TrianglesTemp(Index) = i;

    Index = (Quads == RCM(i));
    QuadsTemp(Index) = i;
end
Elements = ElementsTemp;
msh.Elements = Elements;
HalfTriElements = TriTemp;

SurfTriangles(:,1:3) = TrianglesTemp;
SurfQuads(:,1:4) = QuadsTemp;
% Reorder triangle nodepos to omit nodes that aren't on the triangle face

% TriNodes = unique(TriTemp);
% X = (1:NumElements)';
% RemoveIndex = false(NumElements,1);
% n = 1;
% i = 1;
% while i < length(TriNodes)
%
%     if X(n) ~= TriNodes(i)
%         RemoveIndex(n) = 1;
%     else
%         i = i + 1;
%     end
% n = n + 1;
% end
% X(RemoveIndex) = [];





%%
HalfTriElements = TriTemp;
% msh.TriElements = TriElements;
% Reorder Quads
QuadsTemp = Quads;
for i = 1:NumNodes

    Index = (Quads == RCM(i));
    QuadsTemp(Index) = i;
end
Quads = QuadsTemp;

% Reorder Triangles
TriTemp = Triangles;
for i = 1:NumNodes

    Index = (Triangles == RCM(i));
    TriTemp(Index) = i;
end
Triangles = TriTemp;
msh.RCM = RCM;
% Produce Plot
subplot(121)

spy(NodeConnections)
title('Before RCM')
subplot(122)
spy(NodeConnections(RCM,RCM))
title('After RCM')
% Reorder Node associtivity matrix
NodeConnections = NodeConnections(RCM,RCM);

% Reorder longitudinal connections vector
LongConnections = LongConnections(RCM,:);

LongEndNodes = sum(LongConnections,2) == 1;
NodeConnections = sparse(NodeConnections - eye(size(NodeConnections)));
drawnow
%% Compute SCV and CV volumes
NodeClasses
% FaceBoundaryAreas
% BoundaryAreasScript;
% [BoundaryFaceAreas, BoundaryFaces] = ComputeBoundaryAreas(NumNodes, NodePos,...
%     Elements,Quads,Triangles,MAX);
[CenterPos] = ComputeCenters(NodePos, Elements);

[ElementMidX, ElementMidY, ElementMidZ] = ComputeMidPoints(NodePos,Elements...
    ,NumElements);

[SubVolumes, FullVolumes,SubAreas,~,NormalVectorsX,NormalVectorsY,...
    NormalVectorsZ,EdgeIndex,SCVEdgeAreas] = ...
    ComputeVolumes(NodePos, Elements, CenterPos, ElementMidX, ElementMidY, ElementMidZ,MIN,MAX,NodeClass,NodeType);
%%

ExpectedArea = 2 * (MAX(1) * MAX(2) + MAX(2) * MAX(3) + MAX(1) * MAX(3));
ExpectedVolume = prod(MAX);
% if abs(ExpectedArea - sum(sum(SCVEdgeAreas(:,[3 5 7])))) > sqrt(eps)
%    error('The returned surface area are does not match what is expected.\n     Expected Area = %g\n     Returned Area = %g\n'...
%        ,ExpectedArea,sum(BoundaryAreas(:)))
% end
%
% if abs(ExpectedVolume - sum(FullVolumes(:))) > sqrt(eps)
%    error('The returned volume (of the CVs) are does not match what is expected.n    Expected Volume = %g\n    Returned Volume = %g\n',ExpectedVolume, sum(FullVolumes))
% end
%
% %
% if abs(ExpectedVolume - sum(SubVolumes(:))) > sqrt(eps)
%    error('The returned volume (of the SCVs) are does not match what is expected.\n    Expected Volume = %g\n    Returned Volume = %g\n',ExpectedVolume, sum(FullVolumes))
% end

clf
% Plot_Mesh3D_Shell(NodePos(:,[1 3 2]),msh,[0.2 0.8 0.5].*ones(NumNodes,1))
% Plot_Mesh3D_Shell(msh,[0.2 0.8 0.5].*ones(NumElements,1),1,1)
% Plot_Mesh3D_Shell(msh,[0.2 0.8 0.5].*ones(NumElements,1),NodePos,1,1)
PlotNodes = NodePos(:,[1 3 2]);
PlotNodes(:,3) = 0.045 - PlotNodes(:,3);
Plot_Mesh3D(PlotNodes,Elements,([92 211 217])/255.*ones(NumElements,1),1,1)
view([-45 25])
% view([0 0])
xlabel('X Direction [Metres]')
ylabel('Z Direction [Metres]')
zlabel('Y Direction [Metres]')


view([0 0])
pbaspect([2 1 1])
zlim([0 0.045])

drawnow
pause(0.5)

%%
% Compute boundary areas
SurfQuads(:,5) = SurfQuads(:,5) - 10;
SurfTriangles(:,4) = SurfTriangles(:,4) - 10;

SurfTriArea = zeros(length(SurfTriangles),1);
SurfQuadArea = zeros(length(SurfQuads),1);


% Compute boundary areas for the triangles
for i = 1:length(SurfTriangles)

    x = NodePos(SurfTriangles(i,1:3),1);
    y = NodePos(SurfTriangles(i,1:3),2);
    SurfTriArea(i) = 1/3 * polyarea(x,y);
end

% Compute boundary areas for the Quads
for i = 1:length(SurfQuads)

    switch SurfQuads(i,5)
        case 1
            x = NodePos(SurfQuads(i,1:4),1);
            z = NodePos(SurfQuads(i,1:4),3);
            SurfQuadArea(i) = 1/4 * polyarea(x,z);
        case 2
            x = NodePos(SurfQuads(i,1:4),2);
            z = NodePos(SurfQuads(i,1:4),3);
            SurfQuadArea(i) = 1/4 * polyarea(x,z);
        case 3
            x = NodePos(SurfQuads(i,1:4),1);
            z = NodePos(SurfQuads(i,1:4),3);
            SurfQuadArea(i) = 1/4 * polyarea(x,z);
        case 4
            x = NodePos(SurfQuads(i,1:4),2);
            z = NodePos(SurfQuads(i,1:4),3);
            SurfQuadArea(i) = 1/4 * polyarea(x,z);
        otherwise
    end

end

% Get the normal vectors for each face

SurfTriNormalVectors = zeros(length(SurfTriangles),3);
SurfQuadNormalVectors = zeros(length(SurfQuads),3);

%Quads
QuadFaceNormals = [0 -1 0;
    -1 0 0;
    0 1 0;
    1 0 0];

% Four Quad Faces
for i = 1:4
    SurfQuadNormalVectors(SurfQuads(:,5) == i,:) = ...
        ones(sum(SurfQuads(:,5) == i),1) .* QuadFaceNormals(i,:);
end

% Two end faces 5 at z = 0, 6 at z = Mz
SurfTriNormalVectors(SurfTriangles(:,4) == 5,:) = ...
    ones(sum(SurfTriangles(:,4) == 5),1) .* [0 0 -1];

SurfTriNormalVectors(SurfTriangles(:,4) == 6,:) = ...
    ones(sum(SurfTriangles(:,4) == 6),1) .* [0 0 1];
%%
MeshGeoProps = struct();

MeshGeoProps.NodePos = NodePos;
MeshGeoProps.Elements = Elements;

MeshGeoProps.Triangles = SurfTriangles(:,1:3);
MeshGeoProps.Quads = SurfQuads(:,1:4);
MeshGeoProps.TrianglesSurfIndex = SurfTriangles(:,4);
MeshGeoProps.QuadsSurfIndex = SurfQuads(:,5);
MeshGeoProps.SurfTriArea = SurfTriArea;
MeshGeoProps.SurfQuadArea = SurfQuadArea;
MeshGeoProps.SurfTriNormalVectors = SurfTriNormalVectors;
MeshGeoProps.SurfQuadNormalVectors = SurfQuadNormalVectors;

MeshGeoProps.ScaleFacts = msh.ScalingFacts;

MeshGeoProps.NormalVecsX = NormalVectorsX;
MeshGeoProps.NormalVecsY = NormalVectorsY;
MeshGeoProps.NormalVecsZ = NormalVectorsZ;
MeshGeoProps.CenterPos = CenterPos;

MeshGeoProps.SubVolumes = SubVolumes;
MeshGeoProps.SubAreas = SubAreas;
MeshGeoProps.FullVolumes = FullVolumes;


%%
index_vector = 1:NumNodes;

% Duplicate the index vector
duplicated_index = repelem(index_vector, 2)';
InterlacedFullVolumes = FullVolumes (duplicated_index);
MeshGeoProps.InterlacedFullVolumes = InterlacedFullVolumes;
toc
%% Getting the Pith location and the grain angle!
Imstr = 'CLTImages\BigCLTMid.png';
[XY0,ElementPhi,NodalPhi] = PithLocator(MeshGeoProps,Imstr,MinEndPoints,MaxEndPoints, ...
    SplineCell,MaxCoeffs,MinCoeffs);

MeshGeoProps.ElementalGrainAngle = ElementPhi;
MeshGeoProps.NodalGrainAngle = NodalPhi;
%%
[MeshingSplines,kout] = ExtractLate2EarlySplines(XY0,SplineCell,MinEndPoints,MaxEndPoints,msh.MAX(1),msh.MAX(2));

[TriElementRingsFromPith,ReorderedSplines] = CalculateGrowthRingElements(msh,MinEndPoints,MaxEndPoints, ...
    SplineCellSave,MaxCoeffs,MinCoeffs,XY0,Triangles,NodePos,kout);
pbaspect([1 2.5 1])

%%
% 
% clf
% hold on
% Smoothness = 1-1e-8;
% xp1 = [ReorderedSplines{4}.p.xp; ReorderedSplines{5}.p.xp];
% yp1 = [ReorderedSplines{4}.p.yp; ReorderedSplines{5}.p.yp];
% 
% plot(ReorderedSplines{4}.p.xp,ReorderedSplines{4}.p.yp,'kx')
% plot(ReorderedSplines{5}.p.xp,ReorderedSplines{5}.p.yp,'ks')
% plot(xp1,yp1)
% 
% [fitresult, ~] = createFit2(xp1, yp1,Smoothness);
% fitresult.p
% StartX = min(xp1);
% EndX = max(xp1);
% 
% xplot = linspace(StartX,EndX,100)';
% yplot = fnval(fitresult.p, xplot);
% 
% ReorderedSplines{4}.p = fitresult.p;
% ReorderedSplines{4}.p.xp = xplot;
% ReorderedSplines{4}.p.yp = yplot;
% 
% ReorderedSplines = ReorderedSplines(1:4);
% 
% 
% plot(xplot,yplot,'b','linewidth',2)



%%


[TriRPL] = ComputeRadialPropLength(msh,XY0,Triangles,NodePos, ...
    ReorderedSplines,TriElementRingsFromPith);
pbaspect([1 2.5 1])
%%
BetaParams = [573.3 262 1.1 5.1 9.8];
[TriRho0] = DensityLogisticFit(TriRPL,BetaParams);



%%



[PrismRho0] = Extrude2DTrito3DTriPrism(Elements,Triangles,NodePos,TriRho0);
clf
map = flipud(copper(256));
map = (map(1:end-30,:));
ElementCol = map(round(rescale(PrismRho0)*(length(map)-1))+1,:);
Plot_Mesh3D(NodePos(:,[1 3 2]),Elements,ElementCol,1,1)
colormap(map)
clim(BetaParams([2 1]))
colorbar
view(3)

MeshGeoProps.ElementRho0 = PrismRho0;


%% Volume Averaging (get nodal quantity)

[NodalRho0] = ComputeVolumeAveraging(Elements, SubVolumes,FullVolumes,PrismRho0);
MeshGeoProps.NodeRho0 = NodalRho0;

%%



toc

% [SubVolumes, FullVolumes,SubAreas,~,NormalVectorsX,NormalVectorsY,...
%     NormalVectorsZ,EdgeIndex,SCVEdgeAreas]

% clearvars -except MeshGeoProps MinEndPoints MaxEndPoints  SplineCellSave SplineCell MaxCoeffs MinCoeffs


