function theta = RadOfTime(a,e,t1,mu)


t1 = 30*60;
n = sqrt(mu/a^3);
fun = @(E) E-e*sin(E)-n*t1;
% fun2 = @(E) E-e*sin(E);

Ef = fsolve(fun,pi/4);
theta = 2*atan(sqrt((1+e)/(1-e))*tan(Ef/2));

end
