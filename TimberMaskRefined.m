%% Inital mask thresholding
%  this script produces the image mask that can later be used in the image
%  segmentation script to identify and isolate each of the timber rings.
%
%
% The first section of the code reads in the desired image and converts it
% to grayscale. Then a guassian filter is applied to the full scale image
% to intially reduce noise

clf 
tic
ImageStr = 'BigCLTBot.png';
% Image = imread('CLTImages\BigCLTBot.png');
Image = (imread(ImageStr));
% Image = Image(1:end-50,:,:);
% [650 350 500 200]
% Image = Image(375:(350+200),650:(650+500),:);|
GrayImage = (rescale(mean(Image,3)));
% GrayImage = GrayImage;
% GrayImage = (imresize(GrayImage,[1385,2702]));

% GrayImage = rot90(GrayImage);
Image2 = Image;
GrayImage2 = (rescale(mean(Image2,3)));
GrayImage2 = GrayImage2;

RowCutoff = 25;
ReductionFactor = 16;
Threshold = 0.625;
GuassianSigma = 2;
GrayImage = imgaussfilt(GrayImage,GuassianSigma);
GrayImage2 = imgaussfilt(GrayImage2,0.05);
% GrayImage = GrayImage(RowCutoff:(end-RowCutoff),:);

SmallGrayImage = imresize(GrayImage,round(size(GrayImage) / ReductionFactor));
SmallGrayImage2 = imresize(GrayImage2,round(size(GrayImage2) / ReductionFactor));
Mask = double(SmallGrayImage) <= Threshold;


[NumLines, NumCols] = size(Mask);

imagesc((SmallGrayImage2));
colormap gray
hold on
% T = adaptthresh(SmallGrayImage, 1);
% BW = imbinarize(1-SmallGrayImage,T);
% Mask = BW;
spy((Mask),10)
hold off
pbaspect([fliplr(size(Mask)) 1])
xlabel('Horizontal Direction [Pixels]')
ylabel('Vertical Direction [Pixels]')
%% Steepness Vertical
% This stage is to remove all the small artifacts in the mask and hopefully
% only the ring data - (Which will most likely need to be cleaned up in later sections)
BackStep = 5;
WeightArray = linspace(1,1,BackStep);
WeightArray = [WeightArray 1 WeightArray];
VerticalSteepness = zeros(NumLines,NumCols);
HorizontalSteepness = zeros(NumLines,NumCols);

% Loop through each pixel and caulaute the weights of the pixels infront
% and behind in the vertical direction
for j = 1:NumCols
    for i = (BackStep + 1):(NumLines - BackStep)
        VerticalSteepness(i,j) = sum(Mask((i - BackStep):(i + BackStep),j).*(WeightArray'));
%         VerticalSteepness1(i,j) = sum(Mask((i - BackStep):i,j).*(WeightArray'));
    end
end

% Do the same but for the horizontal steepness
for j = 1:NumLines
    for i = (BackStep + 1):(NumCols - BackStep)
        HorizontalSteepness(j,i) = sum(Mask(j,(i - BackStep):(i + BackStep)).*(WeightArray));
%         VerticalSteepness1(i,j) = sum(Mask((i - BackStep):i,j).*(WeightArray'));
    end
end


Steepness = VerticalSteepness + HorizontalSteepness;

% Apply adjustments to the steepness arrays to reduce the effect of the
% banding. 
BorderArray = ones(size(Steepness));
% Generate array where the banded regions are 1
BorderArray(BackStep+1:(end-BackStep),BackStep+1:(end-BackStep)) = 0;
% Multiply the banded areas by 2
Steepness(BorderArray==1) = Steepness(BorderArray==1)*2;

% Do the same however for the corners
CornerArray = ones(size(Steepness)); % generate logical corner array
CornerArray(BackStep+1:(end - BackStep),:) = 0;
CornerArray(:,BackStep+1:(end - BackStep)) = 0;
% if the mask is already defined in the corner set those steepness to a
% high value
Steepness(CornerArray == 1) = max(Steepness(:)) * Mask(CornerArray == 1);

clf
ax(1) = subplot(221);
imagesc(VerticalSteepness)
title('Vertical')
colorbar
pbaspect([fliplr(size(Mask)) 1])
caxis([0 8])

ax(2) = subplot(222);
imagesc(HorizontalSteepness)
title('Horizontal')
colorbar
pbaspect([fliplr(size(Mask)) 1])
caxis([0 8])

ax(3) = subplot(223);
imagesc(Steepness)
title('Both')
colorbar
pbaspect([fliplr(size(Mask)) 1])
colormap jet
caxis([0 15])

ax(4) = subplot(224);
MaskSteep = 1 * (Steepness > BackStep);
imagesc(SmallGrayImage)
title('Updated Mask')
hold on
spy(MaskSteep,3)
colormap(ax(4),gray)
pbaspect([fliplr(size(Mask)) 1])
colorbar
xlabel('')
%% Morpohologial operations
clf
MaskSteep  =Mask;
% Create the morphological structing element and then apply that to fill
% the holes in the input mask
r = 3; % Radius of the square structuing element
SE = strel('disk',r); % create a square structuring element 
MaskClose = imclose(MaskSteep>0,SE); % Applying the close operation
r = 2;
% SE = strel('disk',1);
% MaskClose = imerode(MaskClose>0,SE);
% SE = strel('line',2,0);
% % MaskClose = imerode(MaskClose>0,SE);
% SE = strel('line',3,90);
% MaskClose = imdilate(MaskClose>0,SE);
% 
MaskClose = medfilt2(MaskClose, [3,3]);
MaskClose(1,1) = 1;
% subplot(222)
% MaskClose(1:10,1:10) = 0;
imagesc(flipud(SmallGrayImage2));
hold on
% MaskClose(115,225) =1;
spy(flipud((MaskClose)),8)
% title('Final Mask - Reduced ')
hold off
xlabel('Horizontal Direction [Pixels]')
ylabel('Vertical Direction [Pixels]')
% Mask = MaskClose;
colormap gray
set(gca,'YDir','normal')
pbaspect([fliplr(size(Mask)) 1])

Mask = MaskClose;
% M1 = imgaussfilt(double(Mask),3);
% Mask = double(M1 > 0.4);
% Switch to a disk structuuring element and apply the erode operation to
% remove

sum(Mask(:) > 0)