function [Xp,Yp,Zp,X,Y,Z] = plotOrbit(a,e,in,OM,w,theta)


r = zeros(length(theta),3);
for i=1:length(theta)

p = a*(1-e^2);

r_con = p/(1+e*cos(theta(i)));

r(i,:) = r_con*[cos(theta(i)),sin(theta(i)),0];

end
Xp = r(:,1);
Yp = r(:,2);
Zp = r(:,3);




R3_1 = [cos((OM)) sin((OM)) 0;-sin((OM)) cos((OM)) 0;0 0 1];
R1_2 = [1 0 0;0 cos(in) sin(in);0 -sin(in) cos(in)];
R3_3 = [cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];

for i=1:length(theta)

r(i,:) = [R3_3*R1_2*R3_1]'*r(i,:)';

end

X = r(:,1);
Y = r(:,2);
Z = r(:,3);


end