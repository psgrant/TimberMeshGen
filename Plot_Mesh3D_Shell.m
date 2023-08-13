function Plot_Mesh3D_Shell(msh,NodeCol,NodePos,FAlpha,EAlpha)

Triangles = msh.TRIANGLES(:,1:3);
Rectangles = msh.QUADS(:,1:4);
RCM = msh.RCM;
NumNodes = msh.nbNod;
% Nodes3D = msh.POS(msh.RCM,[1 3 2]);
% NodePos = msh.POS(:,[1 3 2]);
% NodePos
% Reoder Nodes
% Nodes3D = NodePos(RCM,:);
Nodes3D = NodePos(:,[1 3 2]);
% Reorder Rectangles
QuadsTemp = Rectangles;
for i = 1:NumNodes
    
    Index = (Rectangles == RCM(i));
    QuadsTemp(Index) = i;
end
Rectangles = QuadsTemp;

% Reorder Triangles
TriTemp = Triangles;
for i = 1:NumNodes
    
    Index = (Triangles == RCM(i));
    TriTemp(Index) = i;
end
Triangles = TriTemp;






% Triangles=[];
% Rectangles=[];
% Triangles=[Triangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,3)]];
% Triangles=[Triangles ; [Prisms(:,4) Prisms(:,5) Prisms(:,6)]];
% Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,5)  Prisms(:,4)]];
% Rectangles=[Rectangles ; [Prisms(:,2) Prisms(:,3) Prisms(:,6)  Prisms(:,5)]];
% Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,3) Prisms(:,6)  Prisms(:,4)]];


% All the triangular faces
p0 = patch('Faces',Triangles,'Vertices',Nodes3D,'FaceVertexCData',NodeCol,...
    'FaceColor','interp','EdgeColor','black','FaceAlpha',FAlpha,'EdgeAlpha',EAlpha);
% All the rectangular faces
p1 = patch('Faces',Rectangles,'Vertices',Nodes3D,'FaceVertexCData',NodeCol,...
    'FaceColor','interp','EdgeColor','black','FaceAlpha',FAlpha,'EdgeAlpha',EAlpha);
view(3);


end