%% Bounding box polynomial sampling
% clf
hold on
BoxProps = [35 61]/200;
BoxDensity = [60 75 85]; % Length Per Node
NumBoxes = length(BoxProps);
BoxLocsX = [0 M M 0];
BoxLocsY = [0 0 N N];

% Get the x and y positions of the corners of the box starting with the
% bottom left and going counterclockwise
for i = 1:NumBoxes
    Mp = BoxProps(i) * (M+N) * 0.5;
    BoxLocsX(i+1,:) = [Mp, (M - Mp), (M - Mp), Mp];
    BoxLocsY(i+1,:) = [BoxProps(i)*N,BoxProps(i)*N, (1-BoxProps(i))*N, (1-BoxProps(i))*N];
end
rectangle('Position',[0 0 M N],'Linewidth',2);
pbaspect([M N 1])
hold on
axis tight

for i = 1:NumBoxes
    plot([BoxLocsX(i,:),BoxLocsX(i,1)],[BoxLocsY(i,:),BoxLocsY(i,1)],'Linewidth',2)
end

% rectangle('Position',[BoxLocsX(1,1), BoxLocsY(1,2), BoxLocsX(1,3) - BoxLocsX(1,1),  BoxLocsY(1,4) - BoxLocsY(1,2)],'Linewidth',2,'edgecolor','r');
% for each ring
for i = 1:NumRings
    for k = 1:2
        P = SplineCell{i,k};
        if isempty(P)
            continue
        end
        StartPoint = SplinePoints{i,k}(1,1);
        EndPoint = SplinePoints{i,k}(2,1);
        % need to do everything for the minimum line and the maximum line
        % Get the polynomial coeffs, start and end points of the polyonials
        P = MinCoeffs(i,:);
        StartPoint = MinEndPoints(i,1);
        EndPoint = MinEndPoints(i,3);

        x = StartPoint:EndPoint;
        y = polyval(P,x);

        plot(x,y,'Linewidth',2)


    end
end

%%
% hold off
RingZoning
ExportGeoRd2