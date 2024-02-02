% Load code and mex c file for alignment
addpath('./code') 
mex "./code/DynamicProgrammingQ.c"

%% 
clear all; close all;

load('./data/increasing_noise_example.mat')

% There are 3 options for noise_ind
    % 1: low noise
    % 2: medium noise
    % 3: high noise
noise_ind = 1;

% plot example point clouds
figure(1)
plot(rawData(1,:,1,noise_ind),rawData(2,:,1,noise_ind),'bo')
set(gca,'fontsize',18)
ylim([-2,2])
xlim([-2,2])
axis square

land_mat_slice = land_mat(:,:,:,noise_ind);
% Align Lanscapes
t = linspace(0,1,size(land_mat_slice,2));
[muQ,gam] = align_landscapes(t',land_mat_slice);

% Show landscapes for all observations
figure(2)
plot(t,squeeze(land_mat_slice(1,:,:)),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.6])

% Show aligned landscapes
figure(3)
landAligned = zeros(size(land_mat_slice));
for i=1:size(landAligned,3)
        landAligned(:,:,i) = Group_Action_by_Gamma_Coord(land_mat_slice(:,:,i),gam(:,i));
end
plot(t,squeeze(landAligned),'linewidth',2)
set(gca,'fontsize',18)
ylim([0,.6])

% Show reparameterizations used for alignment
figure(4)
plot(t,squeeze(gam),'linewidth',2)
set(gca,'fontsize',18)
axis square

% Compare pointwise mean to elastic mean
figure(5)
plot(t,mean(land_mat_slice,3),'r','linewidth',2) % plot pointwise mean
hold on
plot(t,q_to_curve(ProjectC(muQ)),'b','linewidth',2) % plot elastic mean
hold off
legend(["pointwise","elastic"])
ylim([0,.6])
set(gca,'fontsize',18)

% Plot noisy diagrams
figure(6)
plot(squeeze(diagPoints(:,1,:,noise_ind)),squeeze(diagPoints(:,2,:,noise_ind)),'bo');hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)


% Apply Warping and plot denoised diagrams
diagPointsWarped = nan(size(diagPoints(:,:,:,noise_ind)));
for i = 1:20
    diagPointsWarped(:,:,i) = interp1(gam(:,i),t,diagPoints(:,:,i,noise_ind));
end
figure(7)
plot(squeeze(diagPointsWarped(:,1,:)),squeeze(diagPointsWarped(:,2,:)),'bo')
hold on
plot(t,t,'k'); hold off
xlim([0,1]);ylim([0,1])
axis square
set(gca,'fontsize',18)


