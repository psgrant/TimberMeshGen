
clf
% ICol = imread('PineScan.png');
ICol = imread('Wood1ChaulkCombinedColour.jpg');

IGray = (im2gray(ICol));
IGray2 = IGray;

imagesc(flipud(ICol))

set(gca,'YDir','normal')

x0 = 500;
xw = 4000 - x0;
y0 = 1000;
yw = 1550-y0;

rectangle('Position',[x0 y0 xw yw],'EdgeColor','y','LineWidth',2)

IGray2 = IGray(y0:(y0+yw),x0:(x0+xw));

xlabel('Horizontal Distance [Pixels]')
ylabel('Vertical Distance [Pixels]')
%%
clf
clc

% IGray2 = modifiedImage;
% IGray2 = IGray2(:,200:3300);
% modifiedImage = modifiedImage(:,200:3300);
T = adaptthresh(IGray2, 0.65);
BW = imbinarize(IGray2,T);
BW2 = BW;
% subplot(311)
% imshow(~BW)
% subplot(312)
% imshow(~modifiedImage)
% subplot(313)
% 
% imshowpair(~BW,~modifiedImage)


%%
clf
subplot(311)
xlocs = zeros();
imagesc(IGray2(:,200:3300))
title('Input Image')
set(gca,'YDir','normal')
% rectangle('Position',[3 70 595 130],'EdgeColor','r','LineWidth',2)
pbaspect([size(IGray2,[2 1]) 1])
xlabel('Horizontal Distance [Pixels]')
ylabel('Vertical Distance [Pixels]')


subplot(312)
I = 120;
% BW = IGray2 > I;
BW =  1 - BW2(:,200:3300);
% BW =  1 - modifiedImage;
imagesc((1-BW))
set(gca,'YDir','normal')
pbaspect([size(IGray2,[2 1]) 1])
colormap gray
ImWidth = size(BW,2);
title(sprintf('Image Mask'))
ImHeight = size(IGray2,1);
BarWidthArray = 50;
xlocsCell = cell(length(BarWidthArray),1);
DensityCell = cell(length(BarWidthArray),1);
hold on
% set(gca,'XTick',[])
% set(gca,'YTick',[])
xlabel('Horizontal Distance [Pixels]')
ylabel('Vertical Distance [Pixels]')
for n = 1:length(BarWidthArray)


    BarWidth = BarWidthArray(n);
    xlocs = zeros(ceil(ImWidth/BarWidth),1);
    NumBars = floor(ImWidth / BarWidth);
    DensityArray = zeros(NumBars,1);
    for i = 1:( NumBars)
        xlocs(i) = (i-1)*BarWidth+1;
        DensityArray(i) = sum(BW(:,(i-1)*BarWidth+1:i*BarWidth),'all');
            xline(xlocs(i),'r','linewidth',1.75)
    end

%     if i == NumBars
%         xlocs(i+1) = ImWidth;
%         DensityArray(i+1) = sum(BW(:,(i+1)*BarWidth+1:end),'all');
%     end
    xlocs = xlocs(1:length(DensityArray));
    xlocsCell{n} = xlocs ./ max(xlocs);
    DensityArray = DensityArray ./ (BarWidth*ImHeight)*1530/1.7;
    DensityCell{n} = DensityArray;

end
xline(xlocs(i)+15,'r','linewidth',1.75)
xlocs1 = cell2mat(xlocsCell);
DenseArray = (cell2mat(DensityCell));
% DensityArray = smooth((DensityArray));
% DensityArray = DensityArray ./ (BarWidth*ImHeight)*1530;
subplot(313)
plot(xlocs,DensityArray,'k.','markersize',10);
pbaspect([3 1 1])
grid on
grid minor
xlim([0 max(xlocs)])
title('Density Over the Image [Kg m^{3}]')
xlabel('Horizontal Distance [Pixels]')
pbaspect([size(IGray2,[2 1]) 1])
ylabel('Density [kg m^{-3}]')
% xlim([0 600])
% ylim([250 700])

%%
clf
plot(rescale(xlocs1),DenseArray,'k.','Linewidth',2);
pbaspect([3 1 1])
ylim([250 700])
grid on
grid minor
title('Density Over the Image [Kg m^{3}]')
xlabel('Nondimensionalised Distance')
ylabel('Density [kg m^{-3}]')

%%

clf
subplot(313)
x3 = rescale(x3);
% y3 = DenseArray;
plot(x3,y3,'k.','MarkerSize',10)
[fitresult, gof] = LogisticFit(x3, y3);
gof
hold on
title('Density Over the Image [Kg m^{3}]')
xlabel('Fractional Radial Position')
ylabel('Density [kg m^{-3}]')
x = 0:0.01:1;
plot(x,fitresult(x),'r','linewidth',2)
pbaspect([size(IGray2,[2 1]) 1])
grid on
grid minor
subplot(312)
% x3 = rescale(xlocs1);
% y3 = DenseArray;
plot(x3,y3,'k.','MarkerSize',10)
hold on
title('Density Over the Image [Kg m^{3}]')
xlabel('Fractional Radial Position')
ylabel('Density [kg m^{-3}]')
pbaspect([size(IGray2,[2 1]) 1])
grid on
grid minor
%%
clf
subplot(133)
x2 = rescale(1:length([22:30,1:19]));
y2 = DenseArray([22:30,1:19]);

plot(rescale(1:length([22:30,1:19])), DenseArray([22:30,1:19]),'k.','MarkerSize',10);
% plot(rescale(A(:,1)), A(:,2)*1000,'k.','MarkerSize',10);
[fitresult, gof] = LogisticFit(rescale(1:length([22:30,1:19])), DenseArray([22:30,1:19]));
hold on
% pbaspect([size(IGray,[2, 1]) 1])


grid on 
grid minor
xlabel('Radial Proportion')
ylabel('Density [kg m^{-3}]')

x = 0:0.01:1;
plot(x,fitresult(x),'r','linewidth',2)
% for i = [-1 1]
% plot(x+i,fitresult(x),'b','linewidth',2)
% end

p11 = predint(fitresult,x,0.95);
% plot(x,p11,'r--','LineWidth',2)
ylim([250 700])
pbaspect([3 1 1])

