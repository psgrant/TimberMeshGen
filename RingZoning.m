
clf
clc
hold on
ColorArray = [1 0 0;
    0 0 1;
    0.5 0 0.5;
    0.5 0.5 0;];
for i = 1:NumBoxes+1
    plot([BoxLocsX(i,:),BoxLocsX(i,1)],[BoxLocsY(i,:),BoxLocsY(i,1)],'color',ColorArray(i,:),'Linewidth',2)
end
hold on


% Find if and if so, where the polynomials intersect the density boxes
PolyIntersections = cell(NumRings*3,1);
Zones = zeros(N,M,NumBoxes);
[X,Y] = meshgrid(1:M,1:N);
mindist = 1000;
for Ring = 1:NumRings
    for k = kout
        InternalZoneRecorded = 0;
        InterArray = [];
        Index = (k-1) * NumRings + Ring;
%         switch k
%             case 1
%                 P = MinCoeffs(Ring,:);
%                 StartPoint = MinEndPoints(Ring,1);
%                 EndPoint = MinEndPoints(Ring,3);
%             case 2
%                 % The max line
%                 P = MaxCoeffs(Ring,:);
%                 StartPoint = MaxEndPoints(Ring,1);
%                 EndPoint = MaxEndPoints(Ring,3);
%             case 3
%                 % The late line
%                 P = LateCoeffs(Ring,:);
%                 StartPoint = LateEndPoints(Ring,1);
%                 EndPoint = LateEndPoints(Ring,3);
%             otherwise
%         end

        P = SplineCell{Ring,k};
        if isempty(P)
            continue
        end
        StartPoint = SplinePoints{Ring,k}(1,1);
        EndPoint = SplinePoints{Ring,k}(2,1);
        x = StartPoint:EndPoint;
        y = P(x);
        % Ring is in a corner
        if length(x) < 1
            continue
        end
        
        
        for Zone = 1:NumBoxes+1
            
            if Zone < NumBoxes+1
                xv = [BoxLocsX(Zone,:), BoxLocsX(Zone,1), NaN, fliplr([BoxLocsX(Zone+1,:), BoxLocsX(Zone+1,1)])];
                yv = [BoxLocsY(Zone,:), BoxLocsY(Zone,1), NaN, fliplr([BoxLocsY(Zone+1,:), BoxLocsY(Zone+1,1)])];
            else
                xv = [BoxLocsX(Zone,:), BoxLocsX(Zone,1)];
                yv = [BoxLocsY(Zone,:), BoxLocsY(Zone,1)];
            end
            xq = StartPoint:EndPoint;
            yq = P(xq);
            [in, on] = inpolygon(xq,yq,xv,yv);
            in = logical(min(in + on,1));
            
                        plot(xv,yv,'LineWidth',2) % polygon
                        axis equal
            
                        hold on
                        plot(xq(in),yq(in),'r+') % points inside
                        plot(xq(~in),yq(~in),'bo') % points outside
                        hold off
                        axis([0 M 0 N])
            in = [0, in, 0]; %#ok
            
            % find leaving the zone
            XLeaving = strfind(in,[1 0]) - 1;
            XLeaving;
            % find entering the zone
            XEntering = strfind(in,[0 1]);
            
            % The whole ring is contained within the one zone
            if sum(in(2:end-1)) == length(x)
                InterArray = [InterArray; Zone, x(1), y(1), x(end), y(end)]; %#ok
                
            else % If the ring passes through multiple zones
                
                for Section = 1:length(XLeaving)
                    
                    if Zone == NumBoxes+1 && InternalZoneRecorded == 0
                        InterArray = [InterArray; Zone, x(XEntering(1)),y(XEntering(1)),...
                            x(XLeaving(1)),y(XLeaving(1))]; %#ok
                        
                        InternalZoneRecorded = 1;
                        
                    else
                        InterArray = [InterArray; Zone, x(XEntering(Section)),y(XEntering(Section)),...
                            x(XLeaving(Section)),y(XLeaving(Section))]; %#ok
                    end
                end
            end
        end
        % Remove sections that are too small
        
        % Sort by values on the x axis
        [~,I] = sort(InterArray(:,2));
        InterArray = InterArray(I,:);
        i = 2;
        
        for i = 2:size(InterArray,1) - 1
            
            SegDist = pdist([InterArray(i,2:3);InterArray(i,4:5)]);
            xp = InterArray(i,2):InterArray(i,4);
            plot(xp, P(xp),'color',ColorArray(InterArray(i,1),:),'linewidth',2);
            
            if SegDist < (0.1*(0.5*(N+M)))
                %                 InterArray(i-1,1) = InterArray(i-1,1) - 1;
                %                 mindist = min(mindist,SegDist);
                if i ~= size(InterArray,1)+1
                    InterArray(i-1,4:5) = InterArray(i+1,4:5);
                    InterArray(i:i+1,:) = [];
                    
                    if k == 3
                        
                    end
                    break
                end
            end
            %             plot(xp, (polyval(P,xp)),'color',ColorArray(InterArray(i,1),:),'linewidth',2);
        end
        
        
        PolyIntersections{Index} = InterArray;
    end
end



%% Visualise
clc
clf
if 1
    
    
    
    hold on
    ColorArray = [1 0 0;
        0 0 1;
        0.5 0 0.5;
        0.5 0.5 0;];
    for i = 1:NumBoxes+1
        plot([BoxLocsX(i,:),BoxLocsX(i,1)],[BoxLocsY(i,:),BoxLocsY(i,1)],'color',ColorArray(i,:),'Linewidth',2)
    end


    for Ring = 1:NumRings

        for k = kout
            Index = (k-1) * NumRings + Ring;
            P = SplineCell{Ring,k};
            if isempty(P)
                continue
            end
            StartPoint = SplinePoints{Ring,k}(1,1);
            EndPoint = SplinePoints{Ring,k}(2,1);
            x = StartPoint:EndPoint;
            y = P(x);

            InterArray = PolyIntersections{Index};
            for j = 1:size(InterArray,1)

                xp = InterArray(j,2):InterArray(j,4);

                plot(xp, P(xp),'color',ColorArray(InterArray(j,1),:),'linewidth',2);

            end

        end
    end
end

% hold off
axis([0 M 0 N])
pbaspect([M, N, 1])
xlabel('Horizontal Distance [Pixels]')
ylabel('Vertical Distance [Pixels]')
