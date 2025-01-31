function[r,v]= OrbCart(a,emag,i,RAAN,omega,theta,mu)
if emag == 0 
    %need to be given u, insert for theta and omega be 0
    p = a*(1-emag^2);
    r_pqw = p/(1+emag*cosd(theta))*[cosd(theta) sind(theta) 0]
    v_pqw = sqrt(mu/p)*[-sind(theta) emag+cosd(theta) 0];
    %only 2 angle transformation
    a = 0;
    B = i;
    y = RAAN;
    Tpqw_ijk = [cosd(-a)*cosd(-y)-cosd(-B)*sind(-a)*sind(-y) cosd(-y)*sind(-a)+cosd(-a)*cosd(-B)*sind(-y) sind(-y)*sind(-B); -cosd(-a)*sind(-y)-cosd(-y)*cosd(-B)*sind(-a) cosd(-a)*cosd(-y)*cosd(-B)-sind(-a)*sind(-y) cosd(-y)*sind(-B); sind(-a)*sind(-B) -cosd(-a)*sind(-B) cosd(-B)]
    r = Tpqw_ijk*r_pqw'
    v = Tpqw_ijk*v_pqw'
    
elseif i == 0
    %need to be given wbar
    p = a*(1-emag^2);
    pospqw = [cosd(theta) sind(theta) 0]';
    rpqw = p/(1+emag*cosd(theta))*pospqw;
    velpqw = [-sind(theta) emag+cosd(theta) 0]';
    vpqw = sqrt(mu/p)*velpqw;
    %only 1 angle transformation
    a = 0;
    B = 0;
    y = wbar
    Tpqw_ijk = [cosd(-a)*cosd(-y)-cosd(-B)*sind(-a)*sind(-y) cosd(-y)*sind(-a)+cosd(-a)*cosd(-B)*sind(-y) sind(-y)*sind(-B); -cosd(-a)*sind(-y)-cosd(-y)*cosd(-B)*sind(-a) cosd(-a)*cosd(-y)*cosd(-B)-sind(-a)*sind(-y) cosd(-y)*sind(-B); sind(-a)*sind(-B) -cosd(-a)*sind(-B) cosd(-B)]
    r = Tpqw_ijk*rpqw
    v = Tpqw_ijk*vpqw
    
elseif (i == 0) && (emag == 0)
    %no real transformation needs to be done.
    %Be given l
    p = a*(1-emag^2);
    r = p/(1+emag*cos(l))*[cos(l) sin(l) 0]
    v = sqrt(mu/p)*[-sin(l) emag+cos(l) 0]
    
else 
    %no singularities 
p = a*(1-emag^2);
pospqw = [cosd(theta) sind(theta) 0]';
rpqw = p/(1+emag*cosd(theta))*pospqw;
velpqw = [-sind(theta) emag+cosd(theta) 0]';
vpqw = sqrt(mu/p)*velpqw;
a = omega;
B = i;
y = RAAN;

Tpqw_ijk = [cosd(-a)*cosd(-y)-cosd(-B)*sind(-a)*sind(-y) cosd(-y)*sind(-a)+cosd(-a)*cosd(-B)*sind(-y) sind(-y)*sind(-B); -cosd(-a)*sind(-y)-cosd(-y)*cosd(-B)*sind(-a) cosd(-a)*cosd(-y)*cosd(-B)-sind(-a)*sind(-y) cosd(-y)*sind(-B); sind(-a)*sind(-B) -cosd(-a)*sind(-B) cosd(-B)]
r = Tpqw_ijk*rpqw
v = Tpqw_ijk*vpqw
end
end