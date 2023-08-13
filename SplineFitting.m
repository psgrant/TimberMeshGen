%%
clf
tiledlayout(1,1)
% % Image = imread('CLTImages/BigCLT003PS.png');
imagesc((flipud(Image)));
set(gca, 'YDir','normal')
colormap gray
hold on

NumRings = length(IsolatedRingsIDs);
MinCoeffs = zeros(length(IsolatedRingsIDs),9);
MaxCoeffs = MinCoeffs;
LateCoeffs = MinCoeffs;
EarlyCoeffs = MinCoeffs;
SplineCell = cell(NumRings,2);
SplinePoints = SplineCell;
CornerArray = zeros(length(IsolatedRingsIDs),2);
MinEndPoints = zeros(length(IsolatedRingsIDs),4);
MaxEndPoints = zeros(length(IsolatedRingsIDs),4);
LateEndPoints = zeros(length(IsolatedRingsIDs),4);
% MinEndPoints = zeros(length(IsolatedRingsIDs),4);


Smoothness = 1e-7;


xlimit = [0 M];
ylimit = [0 N];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);

NumRings = length(IsolatedRingsIDs);
for i = 1:NumRings
    [N,M] = size(GrayImage);
    xy = (T.get(IsolatedRingsIDs(i)));
    x = xy(:,1);
    y = xy(:,2);
    %     x = round(((x-1) * ScaleFacts(2)));
    %     y = round(((y-1) * ScaleFacts(1)));


    %     scatter1 = plot(x,y,'.','color',colvar(i,:),'markersize',12);
    hold on
    k = boundary(x,y,1);
    kx = x(k);
    ky = y(k);

    %     plot(x(k),y(k),'color',colvar(i,:),'linewidth',2);
    %     plot(x(k),y(k),'ko','linewidth',2);

    % We want to loop through each column and find the min and max value of the
    % ring
    ukx = unique(kx);
    kymax = zeros(length(ukx),1);
    kymin = zeros(length(ukx),1);

    % Return the min and max of the boundary points
    for j = 1:length(ukx)

        kmin = min(ky(kx == ukx(j)));
        kmax = max(ky(kx == ukx(j)));


        % Ensure that the top y potoion is on the border for the min
        % line
        if (j == 1) && ((kmin ~= N) && (ukx(1) ~= 0))
            %             kmin = N;
        end

        if (j == length(ukx)) && ((kmax ~= 0) && (ukx(end) ~= M))
            %             kmax = 0;
        end
        % If it's the same boundary point we want to see if it is closer to
        % the min or max line
        if kmin == kmax && (j > 1 && j < length(ukx))
            DistToPrevMax = pdist([ukx(j), kmax;ukx(j-1),kymax(j-1)]);
            DistToPrevMin = pdist([ukx(j), kmax;ukx(j-1),kymin(j-1)]);

            % If it is closer to the previous maximum point we assign it to
            % a maximum
            if DistToPrevMax < DistToPrevMin
                kymax(j) = kmax;
                kymin(j) = min(y(x == ukx(j)));
                % We use the minumim value in the column (not a boundary point)

            else % If the last min point is closer
                kymin(j) = kmax;
                kymax(j) = max(y(x == ukx(j)));
            end
        else
            kymin(j) = kmin;
            kymax(j) = kmax;
        end



    end


    % Removing multiple 0 points from the min array
    ukxmin = ukx;
    ukxmax = ukx;
    ZerosIndex1 = find(kymin(:,1)==0,1,'last')-1;
    ZFirst = find(kymin(:,1)==0,1,'first');

    x = find(round(kymin(:,1))== 0);

    % there are 2 spots we need to remove points from on the top border,
    % the start and the end
    if max(diff(x)) > 1

        y = zeros(max(x),1);
        y(x) = 1;
        [~,a] = findpeaks(1-[0;y;0]);
        FirstZero = a(1) - 3;

        [~,a] = findpeaks([0;y;0]);
        LastZero = a(end) ;
        kymin(LastZero:end) = [];
        ukxmin(LastZero:end) = [];

        kymin(1:FirstZero) = [];
        ukxmin(1:FirstZero) = [];

    else

        if ~isempty(ZerosIndex1)
            if ZFirst > 2
                kymin(ZFirst+1:end) = [];
                ukxmin(ZFirst+1:end) = [];
            else
                kymin(1:ZerosIndex1) = [];
                ukxmin(1:ZerosIndex1) = [];
            end
        end
        ZerosIndex2 = find(kymin(:,1)==N,1,'last')-1;
        ZFirst = find(kymin(:,1)==N,1,'first');
        if ~isempty(ZerosIndex2)
            if ZFirst > 2
                kymin(ZFirst+1:end) = [];
                ukxmin(ZFirst+1:end) = [];
            else
                kymin(1:ZFirst-1) = [];
                ukxmin(1:ZFirst-1) = [];
            end
        end
    end
    % Removing max points that are at the top of the image

    ZerosIndex = find(round(kymax(:,1))== N,1,'last');
    x = find(round(kymax(:,1))== N);

    % there are 2 spots we need to remove points from on the top border,
    % the start and the end
    if max(diff(x)) > 1

        y = zeros(max(x),1);
        y(x) = 1;
        [~,a] = findpeaks(1-[0;y;0]);
        FirstZero = a(1) - 3;

        [~,a] = findpeaks([0;y;0]);
        LastZero = a(end) ;
        kymax(LastZero:end) = [];
        ukxmax(LastZero:end) = [];

        kymax(1:FirstZero) = [];
        ukxmax(1:FirstZero) = [];

    else
        ZerosIndex1 = find(kymax(:,1)==0,1,'last');
        if ZerosIndex1 == length(kymax)
            ZerosIndex1 = [];
        end
        x = find(round(kymax(:,1))== 0);

        if ~isempty(ZerosIndex1)
            if ZerosIndex1 > length(kymax)*0.75
                kymax(ZFirst+1:end) = [];
                ukxmax(ZFirst+1:end) = [];
            elseif ZerosIndex1 < length(kymax)*0.25
                kymax(1:ZerosIndex1-1) = [];
                ukxmax(1:ZerosIndex1-1) = [];
            end
        end
        ZerosIndex2 = find(kymax(:,1)==N,1,'last');
        ZFirst = find(kymax(:,1)==N,1,'first');
        if ~isempty(ZerosIndex2)
            if ZFirst > 2
                kymax(ZFirst+1:end) = [];
                ukxmax(ZFirst+1:end) = [];
            else
                kymax(1:ZerosIndex2-1) = [];
                ukxmax(1:ZerosIndex2-1) = [];
            end
        end
    end

    % returns the start and finish points of the ring for the max line
    [m1,I1] = min(ukxmax);
    [m2,I2] = max(ukxmax);
    MaxEndPoints(i,:) = [m1,kymax(I1),m2,kymax(I2)];

    % Returns the start and finish points of the ring for the max line
    [m1,I1] = min(ukxmin);
    [m2,I2] = max(ukxmin);
    MinEndPoints(i,:) = [m1,kymin(I1),m2,kymin(I2)];


    % Compute and store the polynomials coeefs for the max line
    % Fixed end points too!
    xplot = min(ukxmax)-150:max(ukxmax)+150;


    if length(ukxmax) > 1
        [fitresult, ~] = createFit2(ukxmax, kymax,Smoothness);
        v = fnval(fitresult.p,xplot);

        %         [P] = ComputePolyCoeffs(ukxmax,kymax);
        SplineCell{i,2} = fitresult;
        [xi,yi] = polyxpoly(xplot,v,xbox,ybox);
        [~,I] = sort(xi);
        xi = xi(I);
        yi = yi(I);
        SplinePoints{i,2} = [xi,yi];
        MaxEndPoints(i,:) = [xi(1),yi(1),xi(2),yi(2)];
        plot(xplot,v,'Color',[100 100 255]/255,'linewidth',3)
        plot(ukxmax,kymax,'k.','markersize',10,'linewidth',1.25)
        
    end



    xplot = min(ukxmin)-250:max(ukxmin)+250;

    if length(ukxmin) > 1
        % Compute and store the polynomials for the min line

        
        [fitresult, gof] = createFit2(ukxmin, kymin,Smoothness);
        v = fnval(fitresult.p,xplot);

        SplineCell{i,1} = fitresult
        [xi,yi] = polyxpoly(xplot,v,xbox,ybox);
        [~,I] = sort(xi);
        xi = xi(I);
        yi = yi(I);
        SplinePoints{i,1} = [xi,yi];
        MinEndPoints(i,:) = [xi(1),yi(1),xi(2),yi(2)];
        plot(xplot,v,'Color',[100 100 255]/255,'linewidth',3)
        plot(ukxmin,kymin,'kx','markersize',5,'linewidth',1.25)
        
    end

end

% for both the min and max lines
axis([0 M 0 N])
pbaspect([M, N, 1])
xlabel('Horizontal Distance [Pixels]')
ylabel('Vertical Distance [Pixels]')
SplineCellSave = SplineCell;
grainAngleFitting
[MeshingSplines,kout] = ExtractLate2EarlySplines(XY0,SplineCellSave,MinEndPoints,MaxEndPoints,M,N)
%%
SamplePolynomial

