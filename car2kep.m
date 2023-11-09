function [a,e,i,OM,w,theta] = car2kep(r,v,mu)


K = [0,0,1]';
h = cross(r,v);
i =  rad2deg(acos(h(3)/(norm(h))));
e1 = ((cross(v,h)/mu)-r/norm(r));
e = norm(e1);

a = (-mu/((norm(v)^2/2)-mu/norm(r)))*1/2;


N = cross(K,h)/norm(cross(K,h));

if K == h
 N = [1 0 0];
end

if N(2)>=0
    OM = rad2deg(acos(N(1)));
else
    OM = rad2deg(2*pi- acos(N(1)));
end

if e1(3)>=0
    w = rad2deg(acos((e1*N')/(norm(e1))));
else
    w = rad2deg(2*pi-acos((e1*N')/(norm(e1))));
end

v_rad = (r*v')/norm(r);

if v_rad>=0

    theta = rad2deg((acos((e1*r')/(norm(e1)*norm(r)))));
else
    theta = rad2deg((2*pi-acos((e1*r')/(norm(e1)*norm(r)))));
end

end