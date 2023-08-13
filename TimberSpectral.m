

% Initalise and read in x and y points from mask
% figure
clf
hold on
colormap gray
% MaskFinal = B;
MaskFinal = Mask;
[y,x] = find(MaskFinal);
[NumRows, NumCols] = size(MaskFinal);
[T,CurrentNode] = PaTreeck([x,y]);
ScaleFacts = (size(GrayImage) ./ (size(MaskFinal)-1));
sigma_d = 3;
XorY = 0; % 0 for vertical 1 for horizontal
% Decision varibale to determine if we need to apply clustering (meaning
% that there are more than one ring in the current iteration).

% Set up structure storage for Tree-Like data for each of the
% Form W = Gaussian distribution based on distances
NodeIDs = 0;
CompletedNodeIDs = [];
RingsSegmented = 0;
% Cut off value for detecting number of rings
AverRingCutOff = 1.05;

% Set up the distance and weighting matrix
% vidfile = VideoWriter('Clustering5.avi');
% vidfile.FrameRate = 2;
% open(vidfile);
n = 1;
tic
while RingsSegmented == 0

    %     ClusteringNeeded = 0;

    [AvgRingInter] = CalculateInterceptsMean(T,CurrentNode,NumCols,NumRows,XorY);

    % If the average number of scanned rings is less than 1.2 indicate the
    % end of a branch
    if AvgRingInter < AverRingCutOff

        % Update that that node has reach the end of it's branch, updates
        % the Current Node with the sibling node
        [T,CurrentNode] = T.endOfBranch(CurrentNode);

        % We test to see if the sibling node would be completed
        [AvgRingInter] = CalculateInterceptsMean(T,CurrentNode,NumCols,NumRows,XorY);

        % See if the sibling node is only 1 ring
        if AvgRingInter < AverRingCutOff

            [T,CurrentNode] =  T.endOfBranch(CurrentNode);

        else % otherwise cluster
            ClusteringNeeded = 1;
        end

    else % If the children are complete
        ClusteringNeeded = 1;
    end




    if ClusteringNeeded == 1
        GoodCut = 0;
        % We need to check to see if the cut is between rings or through
        % several cuts we want:
        %      ||  C  ||        ||      ||
        %      ||  C  ||  NOT   CCCCCCCCCC
        %      ||  C  ||        ||      ||

        sigma_d = 4;
        while GoodCut == 0

            xy = T.get(CurrentNode);
            x = xy(:,1);
            y = xy(:,2);
            % Perform spacial spectral clsutering
            tic
            [x1,y1,x2,y2,EigenVals,EigenVecs] = SegmentData(x,y,sigma_d,250);
            toc

            %                         clf
            %                         EigenClustersScript
            % imagesc(SmallGrayImage)
            PlotClusters
            drawnow
            frame = getframe;
            n = n+1;


            % writeVideo(vidfile,frame);
            %             pause(0.5)
            %                         drawnow
            %                         pause(0.5)

            % Calcualte the distance between all the boundary nodes of the
            % two clusters. A bad cut is where a cut is made along a ring
            % and thus the distance between the two boundary nodes will be
            % close together.

            % Finds boundary nodes, the 1 is a strictness factor, so that
            % every boundary node will be selected. Do this for both
            % clusters.
            k1 = boundary(x1,y1,1);
            k2 = boundary(x2,y2,1);

            % Isolate the boundary nodes
            kx1 = x1(k1);
            kx2 = x2(k2);
            kx = [kx1; kx2];

            ky1 = y1(k1);
            ky2 = y2(k2);
            ky = [ky1; ky2];

            % get the length of the boundary nodes for the first
            % clsuercluster
            ClusterIndex = length(kx1);

            BoundaryDistance = squareform(pdist([kx ky]));
            BoundaryDistance(BoundaryDistance == 0) = NaN;

            ClusterArray = zeros(length(kx));
            ClusterIndex = ClusterIndex + 1;
            ClusterArray(ClusterIndex:end,:) = ClusterArray(ClusterIndex:end,:) + 1;
            ClusterArray(:,ClusterIndex:end) = ClusterArray(:,ClusterIndex:end) + 1;


            BoundaryDistance(ClusterArray ~= 1) = NaN;


            [U,~] = min(BoundaryDistance(:));



            % Calculate the number of overlaping y data points between the
            % 2 data sets - this is used to determine if a cut seperates
            % the rings.
            [idy] = ismember(unique(y1),unique(y2));
            [idx] = ismember(unique(x1),unique(x2));
            OverlappingRatiox = sum(idx) / length(idx);
            OverlappingRatioy = sum(idy) / length(idy);
            OverlappingRatio = min(OverlappingRatiox,OverlappingRatioy);
            % We will accept the cut if there is more than 80% overlap
            if U >= 2
                GoodCut = 1;
                fprintf("Clustering Successful\n");
                % Reduce the sigma parameter and loop untill we get an acceptable cut
            else
                sigma_d = sigma_d - 0.25;
                fprintf(2,"Sigma updated to: %g\n",sigma_d);
                if sigma_d < 0 % Make sure sigma doesn't go below 1
                    sigma_d = 2;

                end
            end

        end
        % Create a branch in the tree to store the data points for the 2
        % created clusters
        [T] = T.addNode([x1, y1],CurrentNode);
        [T] = T.addNode([x2, y2],CurrentNode);

        CurrentNodes = T.getChildren(CurrentNode);

        % We will test to see if any of the newly segmented data is a
        % single ring
        [AvgRingInter] = CalculateInterceptsMean(T,CurrentNodes(1),NumCols,NumRows,XorY);

        [~, gof] = createFit(x1, y1);
        gof.rsquare;


        if AvgRingInter < AverRingCutOff || gof.rsquare > 0.9
            [T] = T.endOfBranch(CurrentNodes(1));
        end

        [AvgRingInter] = CalculateInterceptsMean(T,CurrentNodes(2),NumCols,NumRows,XorY);

        [~, gof] = createFit(x2, y2);
        gof.rsquare;

        if (AvgRingInter < AverRingCutOff) || gof.rsquare > 0.9
            [T,CurrentNode] = T.endOfBranch(CurrentNodes(2));

        else
            if T.BranchComplete(CurrentNodes(1)) == 1
                CurrentNode = CurrentNodes(2);
            else
                CurrentNode = CurrentNodes(1);
            end

        end



    end
    [xt,yt] = treelayout(T.Parent);
    % We will find all the nodes that are at the end of the tree
    [I,BottomNodes] = find(yt == min(yt));
    % If the completed number of rings match the number of nodes at the
    % end of the tree break out of the loop
    if (T.BranchComplete(1)) == 1
        break
    end

end
toc
JustSegmented = 1;
% frame = getframe;
% n = n+1;
% writeVideo(vidfile,frame);
% EigenClustersScript
% drawnow
%             pause(0.5)
imagesc(SmallGrayImage)
PlotClusters
drawnow
frame = getframe;
n = n+1;
% writeVideo(vidfile,frame);
% %%
% writeVideo(vidfile,frame);
% frame = getframe;
% % n = n+1;
% writeVideo(vidfile,frame);
% close(vidfile)
%%
%
if JustSegmented == 0
    T = T2;
else
    T2 = T;
    JustSegmented = 0;
end
[xt,yt] = treelayout(T.Parent);

[N,M] = size(GrayImage);
clf
IsolatedRingsIDs = find(yt == min(yt));
colvar = [0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 0.5 0;
    1 1 1;
    0.55 0 0.5;
    0 0.25 0.5;
    0.5 0.25 0];
colvar = [colvar; colvar; colvar];
GrayImage1 = imresize(GrayImage,size(MaskFinal));
imagesc((Image2))
map = linspecer;
hold on
ScaleFacts = (size(GrayImage) ./ (size(MaskFinal)-1));
% map = linspecer;
% map = [linspace(1,0,length(map))',zeros(length(map),1), linspace(0,1,length(map))'];
A = logical(mod(IsolatedRingsIDs,2));
% Mapcol = round(linspace(1,round(length(map)),length(IsolatedRingsIDs(A))));
% Mapcol2 = round(linspace(3*length(map)/4,length(map),length(IsolatedRingsIDs(~A))));
% % Mapcol = [Mapcol1, Mapcol2];
% Mapcol = round(linspace(1,round(length(map)),length(IsolatedRingsIDs)));
% Mapcol = ones(size(Mapcol));
% A1 = imresize(A,size(GrayImage));
Mapcol = round(linspace(1,round(length(map)),9));
for i =  1:length(IsolatedRingsIDs)



    xy = T.get(IsolatedRingsIDs(i));
    x = xy(:,1);
    y = xy(:,2);
    x = ((x-1) * ScaleFacts(2));
    y = ((y-1) * ScaleFacts(1));
    randcol = randi([1 length(map)]);
    %     scatter1 = plot(x,y,'.','color',map(Mapcol(i),:),'markersize',8);
    hold on
    %     scatter1 = plot(x,y,'g.','markersize',12);
    %     scatter1.MarkerFaceAlpha = 0.2;
    %     scatter1.MarkerEdgeAlpha = 0.2;
    %     axis([0 105 0 41])
    k = boundary(x,y,1);


    xlabel('Horizontal Direction [Pixels]')
    ylabel('Vertical Direction [Pixels]')
    %     plot(x(k),y(k),'wo','linewidth',2);
    [X,Y] = meshgrid(linspace(0,M,150),linspace(0,N,500));

    XGrid = X(:);
    YGrid = Y(:);

    xk = x(k);
    yk = y(k);

    if ~(isempty(find(xk==0))) && ~(isempty(find(xk==M))) %#ok

        minIL = find(xk==0);
        [TopMinL, MinTopIL] = min(yk(minIL));

        minIR = find(xk==M);
        [TopMinR, MinTopIR] = min(yk(minIR));
        xBot = xk(minIL(MinTopIL):minIR(MinTopIR));
        yBot = yk(minIL(MinTopIL):minIR(MinTopIR));


        if yBot(1) > yBot(2)
            set(gca,'YDir','Normal')

            xBot1 = xBot(2:3);
            yBot1 = yBot(2:3);
            [LM] = LinearFit(xBot1, yBot1);
            m = LM.p1;
            c = LM.p2;

            y(k(minIL(MinTopIL))) = c;
            %             plot(xk(minIL(MinTopIL):minIR(MinTopIR)),yk(minIL(MinTopIL):minIR(MinTopIR)),...
            %                 'r.','markersize',15)
        end

        if yBot(end) > yBot(end-1)
            set(gca,'YDir','Normal')

            xBot1 = xBot(end-2:end-1);
            yBot1 = yBot(end-2:end-1);
            [LM] = LinearFit(xBot1, yBot1);
            m = LM.p1;
            c = LM.p2;

            y(k(minIR(MinTopIR))) = m * M + c;
            %             plot(xk(minIL(MinTopIL):minIR(MinTopIR)),y(k(minIL(MinTopIL):minIR(MinTopIR))),...
            %                 'r.','markersize',15)
        end
    end
    plot(x(k),y(k),'k','linewidth',1.5);
    in = inpolygon(XGrid,YGrid,x(k),y(k));
    xnew = XGrid(in);
    ynew = YGrid(in);
    scatter1 = plot(xnew,ynew,'.','color',map(Mapcol(i),:),'markersize',8);
    [T] = set(T, [xnew, ynew], IsolatedRingsIDs(i));
end
% axis([0, size(GrayImage,2)+1, 0, size(GrayImage,1)+1])
pbaspect([fliplr(size(GrayImage)) 1]);
set(gca,'YDir','Normal')


%%
% CubicRingFitting