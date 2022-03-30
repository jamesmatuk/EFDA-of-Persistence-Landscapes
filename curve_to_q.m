function [q] = curve_to_q(p)

[n,N] = size(p);

for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N));
end

for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro'));
        q(:,i) = v(:,i)/(L(i)+eps);
end


end