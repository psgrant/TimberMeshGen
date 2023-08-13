clf
imagesc(GrayImage2)
hold on

x1p = ((x1-1) * ScaleFacts(2));
y1P = ((y1-1) * ScaleFacts(1));

x2p = ((x2-1) * ScaleFacts(2));
y2p = ((y2-1) * ScaleFacts(1));
plot(x1p,y1P,'b.','markersize',10)
plot(x2p,y2p,'r.','markersize',10)
k = boundary(x1p,y1P,1);
plot(x1p(k),y1P(k),'b','linewidth',1);
plot(x1p(k),y1P(k),'ko','linewidth',1,'markersize',3);

k = boundary(x2p,y2p,1);
h1 = plot(x2p(k),y2p(k),'r','linewidth',1);
h2 = plot(x2p(k),y2p(k),'ko','linewidth',1,'markersize',3);

% hold off
% title('Clustered Data')

axis([0, size(GrayImage,2)+1, 0, size(GrayImage,1)+1])
pbaspect([fliplr(size(GrayImage)) 1])


CompletedClusters = find(T.BranchComplete);

for i = 1:length(CompletedClusters)
    j = CompletedClusters(i);
    xcp = T.Data{j}(:,1);
    ycp = T.Data{j}(:,2);
    xcp = ((xcp-1) * ScaleFacts(2));
    ycp = ((ycp-1) * ScaleFacts(1));
    k = boundary(xcp,ycp,1);
    h3 = plot(xcp,ycp,'.','color',[0.75 0.75 0],'markersize',10);
    plot(xcp(k),ycp(k),'color',[0.75 0.75 0],'linewidth',1.5);
    plot(xcp(k),ycp(k),'ko','linewidth',1.5,'markersize',3);
end
