% Load code and mex c file for alignment
addpath('./code') 
mex ./code/DynamicProgrammingQ.c

%% Main paper - simulation 1
clear all; close all;

load('./data/main_sim_1.mat')

% Align Lanscapes
t = linspace(0,1,size(land_mat,2));
[muQ,gam] = align_landscapes(t',land_mat);

% Show landscapes for all observations
figure
plot(t,squeeze(land_mat(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.7])


% Show aligned landscapes
figure
landAligned = zeros(size(land_mat));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat(:,:,i),gam(:,i));
end
plot(t,squeeze(landAligned),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.7])


% Show reparameterizations used for alignment
figure
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square


% Compare pointwise mean to elastic mean
figure
plot(t,mean(land_mat,3),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,q_to_curve(ProjectC(muQ)),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.7])
set(gca,'fontsize',18)

% Plot noisy diagrams
figure
plot(diagPoints(:,1),diagPoints(:,2),'bo');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply Warping and plot denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(i,:) = interp1(gam(:,i),t,diagPoints(i,:));
end
figure
plot(diagPointsWarped(:,1),diagPointsWarped(:,2),'bo')
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)


%% Main paper - simulation 2

clear all; close all;

load('./data/main_sim_2.mat')

t = linspace(0,1,length(land_mat));

% Show data from each group
% one circle
i = 7;
plot(rawData(1,:,i),rawData(2,:,i),'bo','markerFaceColor','b')
axis square
xlim([-1.5,3]);ylim([-1.5,1.5])
set(gca,'fontsize',18)
figure
plotLandscape(t,land_mat(:,:,i),[0,.5])
% two interconnected circles
i = 19;
figure
plot(rawData(1,:,i),rawData(2,:,i),'ro','markerFaceColor','r')
axis square
xlim([-1.5,3]);ylim([-1.5,1.5])
set(gca,'fontsize',18)
figure
plotLandscape(t,land_mat(:,:,i),[0,.5])

% Align landscapes
[muQ,gam,qn] = align_landscapes(t,land_mat);

% Compute elastic PCA
clear qnVec
for i = 1:20
qnVec(i,:) = [qn(1,:,i),qn(2,:,i)];
end
[basisQ,scoreQ,latentQ,tsquaredQ,explainedQ,muQ] = pca(qnVec);

% Visualize elastic principal directions
b = 1;
basisVis = zeros(2,length(basisQ)/2);
count = 1;
for sd = [-1,0,1]
basisVis(1,:) = muQ(1:length(t))' + sd*sqrt(latentQ(b))*basisQ(1:length(t),b);
basisVis(2,:) = muQ((length(t)+1):(2*length(t)))' + sd*sqrt(latentQ(b))*basisQ((length(t)+1):(2*length(t)),b);
basisVis = q_to_curve(ProjectC(basisVis));
figure
plotLandscape(t,basisVis,[0,.5])
count = count + 1;
end

% Show elastic PC scores
figure
plot(scoreQ(1:10,1),scoreQ(1:10,2),'bo','markerFaceColor','b')
hold on
plot(scoreQ(11:end,1),scoreQ(11:end,2),'ro','markerFaceColor','r')
hold off
set(gca,'fontsize',18)
xlabel('PC 1')
ylabel('PC 2')

% Compute poitwise PCA
clear landVec
for i = 1:20
landVec(i,:) = [land_mat(1,:,i),land_mat(2,:,i)];
end
[basis,score,latent,tsquared,explained,mu] = pca(landVec);


% Visualize pointwise principal directions
b = 1;
basisVis = zeros(2,length(basis)/2);
count = 1;
for sd = [-1,0,1]
basisVis(1,:) = mu(1:length(t))' + sd*sqrt(latent(b))*basis(1:length(t),b);
basisVis(2,:) = mu((length(t)+1):(2*length(t)))' + sd*sqrt(latent(b))*basis((length(t)+1):(2*length(t)),b);
figure
plotLandscape(t,basisVis,[0,.4])
count = count + 1;
end


% Visualize pointwise PC scores
figure
plot(score(1:10,1),score(1:10,2),'bo','markerFaceColor','b')
hold on
plot(score(11:end,1),score(11:end,2),'ro','markerFaceColor','r')
hold off
set(gca,'fontsize',18)
xlabel('PC 1')
ylabel('PC 2')

% Visualize noisy diagrams
figure
diagPoints(21:30,:) = 0;
plot(diagPoints([1:10,21:30],1),diagPoints([1:10,21:30],2),'bo');hold on
plot(diagPoints([11:20,31:40],1),diagPoints([11:20,31:40],2),'ro');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply warping and display denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(i,:) = interp1(gam(:,i),t,diagPoints(i,:));
    diagPointsWarped(i+20,:) = interp1(gam(:,i),t,diagPoints(i+20,:));
end
figure
plot(diagPointsWarped([1:10,21:30],1),diagPointsWarped([1:10,21:30],2),'bo');hold on
plot(diagPointsWarped([11:20,31:40],1),diagPointsWarped([11:20,31:40],2),'ro');hold on
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)


% Show reparameterizations used for alignment
figure
plot(t,squeeze(gam(:,1:10)),'b','linewidth',2);hold on;
plot(t,squeeze(gam(:,11:20)),'r','linewidth',2);hold on;
set(gca,'fontsize',18)
axis square


