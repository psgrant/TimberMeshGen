function Plot_Mesh3D(Nodes3D,Prisms,ElementCol,FAlpha,EAlpha)
% Plots triangular prism mesh where each element's faces are coloured in
% the same colour. Given an input of ElementCol (size(ElementCol) ==
% size(Prisms)). Note the colours for the input are scaled RGB so each
% RGB value in ElementCol is between 0 and 1.

% Generate
Triangles=[];
Rectangles=[];
Triangles=[Triangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,3)]];
Triangles=[Triangles ; [Prisms(:,4) Prisms(:,5) Prisms(:,6)]];
Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,5)  Prisms(:,4)]];
Rectangles=[Rectangles ; [Prisms(:,2) Prisms(:,3) Prisms(:,6)  Prisms(:,5)]];
Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,3) Prisms(:,6)  Prisms(:,4)]];

% Create Face Colour arrays
TriFaceCol =  [ElementCol; ElementCol];
RectFaceCol = [ElementCol; ElementCol; ElementCol];

% All the triangular faces
patch('Vertices',Nodes3D,'Faces',Triangles,'FaceVertexCData',TriFaceCol,...
    'FaceColor','flat','EdgeColor','black','FaceAlpha',FAlpha,'EdgeAlpha',EAlpha);

% All the rectangular faces
patch('Vertices',Nodes3D,'Faces',Rectangles,'FaceVertexCData',RectFaceCol,...
    'FaceColor','flat','EdgeColor','black','FaceAlpha',FAlpha,'EdgeAlpha',EAlpha);

end