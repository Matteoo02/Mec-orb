function dT = timeOfFlight(a,e,th0,thf,mu)

E0 = 2*atan(sqrt((1-e)/(1+e))*tan(th0/2));
Ef = 2*atan(sqrt((1-e)/(1+e))*tan(thf/2));
if E0<0
    E0 = E0 + 2*pi;
else if Ef<0
        Ef = Ef + 2*pi;
end
end

dM = (Ef-E0)-e*(sin(Ef)-sin(E0));


 dT = dM*sqrt(a^3/mu); %secondi

end