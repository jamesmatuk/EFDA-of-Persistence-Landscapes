function [mu,gam,qn] = align_landscapes(t,land)

% Compute SRVFs
q = zeros(size(land));
for i = 1:size(land,3)
q(:,:,i) = curve_to_q(land(:,:,i));
end

% Compute SRVFs
mnq = mean(q,3);
dist = zeros(size(q,3),1);
for i = 1:size(q,3)
dqq = zeros(size(q,2),1);
for j = 1:length(t)
  dqq(j) =   norm((q(:,j,i) - mnq(:,j)));
end
dist(i) = trapz(t,dqq);
end
[~, min_ind] = min(dist);
mu = q(:,:,min_ind);

% initialize algorithm for computing elastic mean
iter=1;
maxIter = 20;
e = 1e-6;
E=1;
qn = zeros(size(q));
gam = zeros(length(t),size(q,3));

% run algorithm for computing elastic mean
while (E(iter) > e)
    iter=iter+1;
    for i=1:size(q,3)
        gam(:,i) = DynamicProgrammingQ(q(:,:,i),mu(:,:,iter-1),0,0);
        gam(:,i) = (gam(:,i)-gam(1,i))/(gam(end,i)-gam(1,i));
        qn(:,:,i) = curve_to_q(Group_Action_by_Gamma_Coord(q_to_curve(q(:,:,i)),gam(:,i)));
    end
    qm =  mean(qn,3);
    E(iter) = sqrt(trapz(t,(sum((mu(:,:,iter-1)  - qm).^2,1)')));
    mu(:,:,iter) = mean(qn,3);
    if (iter == maxIter)
       break 
    end
end

% center mean within orbit
gamI = SqrtMeanInverse(gam');
mu = curve_to_q(Group_Action_by_Gamma_Coord(q_to_curve(mean(qn,3)),gamI));
for i=1:size(q,3)
        gam(:,i) = DynamicProgrammingQ(q(:,:,i),mu,0,0);
        gam(:,i) = (gam(:,i)-gam(1,i))/(gam(end,i)-gam(1,i));
        qn(:,:,i) = curve_to_q(Group_Action_by_Gamma_Coord(q_to_curve(q(:,:,i)),gam(:,i)));
end
mu =  mean(qn,3);

end