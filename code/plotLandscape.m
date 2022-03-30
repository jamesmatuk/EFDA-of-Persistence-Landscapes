function plotLandscape(t,landscape,zl)


f = zeros(size(landscape,1)*10,length(t));
for j = 1:size(landscape,1)
    f(j*10,:) = landscape(j,:);
end

[m,n] = size(f);
ty = linspace(0,m/10,m)';
tx = t';
% figure
% figure('units','normalized','outerposition',[0 0 1 1])
s = surf(tx,ty,f);
% s.EdgeColor = 'interp';
% s.linewidth = .01
% s.FaceAlpha = .01;
s.LineStyle = '-';
set(gca, 'XDir','reverse')
zlim(zl)
set(gca,'xtick',[])
set(gca,'fontsize',18)
view(170,5)
box on

end