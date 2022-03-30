% Load code and mex c file for alignment
addpath('./code') 
mex ./code/DynamicProgrammingQ.c

%% Main paper - simulation 1 with noise
clear all; close all;

load('./data/main_sim_1_noise.mat')

% Show examples of data and landscapes
t = linspace(0,1,length(land_mat));
for i = [1,10]
figure
plot(rawData(1,:,i),rawData(2,:,i),'bo','markerFaceColor','b');
axis square
xlim([-1.5,1.5]);ylim([-1.5,1.5])

figure
plot(t,squeeze(land_mat(1,:,i)),'b','linewidth',2);
ylim([0,.9])
set(gca,'fontsize',1)
end

% Align landscapes
[muQ,gam] = align_landscapes(t',land_mat);

% Show landscapes 
figure
plot(t,squeeze(land_mat(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.9])

% show aligned landscapes
figure
landAligned = zeros(size(land_mat));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat(:,:,i),gam(:,i));
end
plot(t,squeeze(landAligned),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.9])

% Show reparameterizations used for alignment
figure
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square

% compare pointwise mean to elastic mean
figure
plot(t,mean(land_mat,3),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,q_to_curve(ProjectC(muQ)),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.9])
set(gca,'fontsize',18)


% Show noisy diagrams
figure
for i = 1:20
plot(diagPoints(:,1,i),diagPoints(:,2,i),'bo');hold on
end
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)
% saveas(gcf,['./newSims/figuresFinal/noisecircleDiagrams.png'])

% Apply warping and display denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(:,:,i) = interp1(gam(:,i),t,diagPoints(:,:,i));
end
figure
for i = 1:20
plot(diagPointsWarped(:,1,i),diagPointsWarped(:,2,i),'bo');hold on
end
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

%% Main paper - simulation 2 with noise

clear all; close all;

load('./data/main_sim_2_noise.mat')

t = linspace(0,1,length(land_mat));

% Show examples from each group
% one circle with noise
i = 3;
figure
plot(rawData(1,:,i),rawData(2,:,i),'bo','markerFaceColor','b')
axis square
xlim([-1.5,3]);ylim([-1.5,1.5])
set(gca,'fontsize',18)
figure
plotLandscape(t,land_mat(:,:,i),[0,.75])

% two interconnected circles with noise
i = 11;
figure
plot(rawData(1,:,i),rawData(2,:,i),'ro','markerFaceColor','r')
axis square
xlim([-1.5,3]);ylim([-1.5,1.5])
set(gca,'fontsize',18)
figure
plotLandscape(t,land_mat(:,:,i),[0,.75])

% compute elastic mean
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
plotLandscape(t,basisVis,[0,.75])
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

% Compute pointwise PCA
clear landVec
for i = 1:20
landVec(i,:) = [land_mat(1,:,i),land_mat(2,:,i)];
end
[basis,score,latent,tsquared,explained,mu] = pca(landVec);

% Visualize Principal directions
b = 1;
basisVis = zeros(2,length(basis)/2);
count = 1;
for sd = [-1,0,1]
basisVis(1,:) = mu(1:length(t))' + sd*sqrt(latent(b))*basis(1:length(t),b);
basisVis(2,:) = mu((length(t)+1):(2*length(t)))' + sd*sqrt(latent(b))*basis((length(t)+1):(2*length(t)),b);
figure
plotLandscape(t,basisVis,[0,.75])
count = count + 1;
end

% Visualize Pointwise PC scores
figure
plot(score(1:10,1),score(1:10,2),'bo','markerFaceColor','b')
hold on
plot(score(11:end,1),score(11:end,2),'ro','markerFaceColor','r')
hold off
set(gca,'fontsize',18)
xlabel('PC 1')
ylabel('PC 2')

% Show noisy diagrams
figure
for i = 1:10
   plot(diagPoints(1:10,1,i),diagPoints(1:10,2,i),'bo');hold on
   plot(0,0,'bo');hold on
   plot(diagPoints(1:10,1,i+10),diagPoints(1:10,2,i+10),'ro');hold on 
end
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% Apply warping and show denoised diagrams
diagPointsWarped = zeros(size(diagPoints));
for i = 1:20
    diagPointsWarped(:,1,i) = interp1(gam(:,i),t,diagPoints(:,1,i));
    diagPointsWarped(:,2,i) = interp1(gam(:,i),t,diagPoints(:,2,i));
end
figure
for i = 1:10
   plot(diagPointsWarped(1:10,1,i),diagPointsWarped(1:10,2,i),'bo');hold on
   plot(0,0,'bo');hold on
   plot(diagPointsWarped(1:10,1,i+10),diagPointsWarped(1:10,2,i+10),'ro');hold on 
end
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)

% show reparameterizations used for alignment
figure
plot(t,squeeze(gam(:,1:10)),'b','linewidth',2);hold on;
plot(t,squeeze(gam(:,11:20)),'r','linewidth',2);hold on;
set(gca,'fontsize',18)
axis square
