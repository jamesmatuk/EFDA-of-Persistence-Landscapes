% Load code and mex c file for alignment
addpath('./code') 
mex ./code/DynamicProgrammingQ.c

%% Main Paper - circles

close all;clear all;

load('./data/main_circles.mat')

t = linspace(0,1,length(land_mat));

% Show some example data and landscapes
for i = [10,18]
    figure
    plot(rawData(1,:,i),rawData(2,:,i),'ro','markerFaceColor','r');
    axis square
    xlim([-1,3]);ylim([-1,1])
    set(gca,'fontsize',18)
    figure
    plot(t,land_mat(1,:,i),'r','linewidth',2);hold on
    plot(t,land_mat(2,:,i),'r','linewidth',2);hold off
    ylim([0,.7])
    set(gca,'fontsize',18)
end


% Align landscapes
[muQ,gam] = align_landscapes(t',land_mat);

% show landscapes for all observations
% first landscape function
figure
plot(t,squeeze(land_mat(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.7])
% second landscape function
figure
plot(t,squeeze(land_mat(2,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.31])

% show aligned landscapes
% first landscape function
landAligned = zeros(size(land_mat));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat(:,:,i),gam(:,i));
end
figure
plot(t,squeeze(landAligned(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.7])
% Second landscape function
figure
plot(t,squeeze(landAligned(2,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.31])

% Show reparameterizations used for alignment
figure
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square

% Compare pointwise mean to elastic mean
pMean = mean(land_mat,3);
eMean = q_to_curve(ProjectC(muQ));
% first landscape function
figure
plot(t,pMean(1,:),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,eMean(1,:),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.7])
set(gca,'fontsize',18)
% second landscpae function
figure
plot(t,pMean(2,:),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,eMean(2,:),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.31])
set(gca,'fontsize',18)

% Show noisy diagrams
figure
plot(diagPoints(:,1),diagPoints(:,2),'ro');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply warping and show denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(i,:) = interp1(gam(:,i),t,diagPoints(i,:));
    diagPointsWarped(i+20,:) = interp1(gam(:,i),t,diagPoints(i+20,:));
end
figure
plot(diagPointsWarped(:,1),diagPointsWarped(:,2),'ro')
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)


%% Main Paper - spirals

close all;clear all;

load('./data/main_spirals.mat')

t = linspace(0,1,length(land_mat));
% Plot a few examples of data and landscapes
for i = [18,20]
    figure
    plot(rawData(1,1:1000,i),rawData(2,1:1000,i),'ro','markerFaceColor','r');hold on
    plot(rawData(1,1001:end,i),rawData(2,1001:end,i),'bo','markerFaceColor','b');hold off
    axis square
    xlim([-1.15,1.15]);ylim([-1.15,1.15])
    set(gca,'fontsize',18)

    figure
    plot(t,land_mat(1,:,i),'r','linewidth',2);hold on
    ylim([0,.15])
    set(gca,'fontsize',18)
end

% Align landscapes
[muQ,gam] = align_landscapes(t',land_mat);

% Show landscapes for all observations
figure
plot(t,squeeze(land_mat(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.15])

% Show aligned landscapes
landAligned = zeros(size(land_mat));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat(:,:,i),gam(:,i));
end
figure
plot(t,squeeze(landAligned(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.15])

% Show reparameterizations used for alignment
figure
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square


% compare pointwise mean to elastic mean
pMean = mean(land_mat,3);
eMean = q_to_curve(ProjectC(muQ));
figure
plot(t,pMean(1,:),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,eMean(1,:),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.15])
set(gca,'fontsize',18)

% Show Noisy Diagrams
figure
plot(diagPoints(:,1),diagPoints(:,2),'ro');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply warping and show denoised landscapes
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(i,:) = interp1(gam(:,i),t,diagPoints(i,:));
end
figure
plot(diagPointsWarped(:,1),diagPointsWarped(:,2),'ro')
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

%% Main paper - torus

close all;clear all;

load('./data/main_torus.mat')

% Show example data and corresponding landscapes
t = linspace(0,1,length(land_mat));
for i = [3,19]
    figure
    scatter3(rawData(1,:,i)',rawData(2,:,i)',rawData(3,:,i)','k')
    view(45,70)
    xlim([-3,3]);ylim([-3,3]);zlim([-1,1])
    set(gca,'fontsize',18)
    figure
    plot(t,land_mat(1,:,i),'k','linewidth',2);hold on
    plot(t,land_mat(2,:,i),'k','linewidth',2);hold off
    ylim([0,1.3])
    set(gca,'fontsize',18)
end


% Align landscapes
[muQ,gam] = align_landscapes(t',land_mat);

% Show landscapes for all observations
% first landscape function
figure
plot(t,squeeze(land_mat(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,1.3])
% secon landscape function
figure
plot(t,squeeze(land_mat(2,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.5])

% Show aligned landscapes
landAligned = zeros(size(land_mat));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat(:,:,i),gam(:,i));
end
% first landscape function
figure
plot(t,squeeze(landAligned(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,1.3])
% second landscape function
figure
plot(t,squeeze(landAligned(2,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.5])

% show reparameterizations used for alignment
figure
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square

% compare pointwise mean to elastic mean
pMean = mean(land_mat,3);
eMean = q_to_curve(ProjectC(muQ));
% first landscape function
figure
plot(t,pMean(1,:),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,eMean(1,:),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,1.3])
set(gca,'fontsize',18)
% second landscape function 
figure
plot(t,pMean(2,:),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,eMean(2,:),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.5])
set(gca,'fontsize',18)

% Show noisy diagrams
figure
plot(diagPoints(:,1),diagPoints(:,2),'ko');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply warping and show denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(i,:) = interp1(gam(:,i),t,diagPoints(i,:));
    diagPointsWarped(i+20,:) = interp1(gam(:,i),t,diagPoints(i+20,:));
end
figure
plot(diagPointsWarped(:,1),diagPointsWarped(:,2),'ko')
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

