function pDist = phaseDist(t,g1,g2)

psi1 = sqrt(gradient(g1,t));
psi2 = sqrt(gradient(g2,t));


pDist = acos(trapz(t,psi1.*psi2));

end