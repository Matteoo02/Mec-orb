function [r,v] = kep2car(a,e,i,OM,w,theta,mu)

p = a*(1-e^2);

r_con = p/(1+e*cos(theta));


r = r_con*[cos(theta),sin(theta),0]';
v = sqrt((mu/p))*[-sin(theta),(e+cos(theta)),0]';


R3_1 = [cos((OM)) sin((OM)) 0;-sin((OM)) cos((OM)) 0;0 0 1];
R1_2 = [1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];
R3_3 = [cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];

r = [R3_3*R1_2*R3_1]'*r;
v = [R3_3*R1_2*R3_1]'*v;
r = r';
v = v';
end



