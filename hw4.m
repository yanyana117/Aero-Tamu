%Find outer radius as a cunction of changing inner radius
%assuming constant skeletal loading and consequently constant
%section modulus Z
ro_init = 1.70; %cm, for males
%ro_init = 1.40; %cm, for females
ri_init = 1.20; %cm, for males
%ri_init = 0.90; %cm, for females
ri_rate = 4*0.004; %cm/year
rho = 1.05; %g/cm^3, solid bone density
Z = 0.25*pi*(ro_init^4 - ri_init^4)/ro_init; %section modulus



bone_results = zeros(70,7); %matrix to hold results for each year
ro = ro_init; %initialize outer radius



for t = 1:70,
    bone_results(t,1) = t + 25; % age in years
    ro_side = 0.25*pi*ro^4 - Z*ro; % one side of equation for Z
    ri = ri_init + ri_rate*t;
    ri_side = 0.25*pi*ri^4; % other side of equation for Z
    
    while (ro_side <= ri_side),
        ro = ro + 0.00005; % make tiny increment in ro if one
        ro_side = 0.25*pi*ro^4 - Z*ro; %side of equation is less than
    
    end %other

    area = pi*(ro^2 - ri^2); %once the two sides converge,
    BMD = area*rho/(2*ro); %store results for this t
    bone_results(t,2) = ro; %in matrix
    bone_results(t,3) = ri;
    bone_results(t,4) = ro_side;
    bone_results(t,5) = ri_side;
    bone_results(t,6) = area;
    bone_results(t,7) = BMD;
    t = t + 1; %move on to the next year

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find outer radius as a cunction of changing inner radius
%assuming variable skeletal loading and consequently variable
%section modulus Z
%ro_init = 1.70; %cm, for males
ro_init = 1.40; %cm, for females
%ri_init = 1.20; %cm, for males
ri_init = 0.90; %cm, for females
ri_rate = 0.004; %cm/year
rho = 1.05; %g/cm^3, solid bone density


Z_init = 0.25*pi*(ro_init^4 - ri_init^4)/ro_init; %section modulus
Z = Z_init;

bone_results = zeros(70,8); %matrix to hold results for each year
ro = ro_init; %initialize outer radius

for t = 1:70,
    bone_results(t,1) = t + 25; %age in years
    if t >= 10, %corresponds to age 35
        if t <= 35, %corresponds to age 60
            Z = Z_init*(-0.0120*(t-10)+1);
                %decreased loading ---> decreased Z
                %this expression for changing Z will cause it
                %decrease by 30% over 25 years
        end
    end

ro_side = 0.25*pi*ro^4 - Z*ro; %one side of equation for Z
ri = ri_init + ri_rate*t;
ri_side = 0.25*pi*ri^4; %other side of equation for Z
    
    while (ro_side <= ri_side),
        ro = ro + 0.00005;
        ro_side = 0.25*pi*ro^4 - Z*ro;
    end

if (ro_side >= ri_side+0.004),

    ro = ro - 0.02; %make tiny increment in ro if one
    ro_side = 0.25*pi*ro^4 - Z*ro; %side is less than other

    while (ro_side <= ri_side),
        ro = ro + 0.00005;
        ro_side = 0.25*pi*ro^4 - Z*ro;
    end

end
    
    area = pi*(ro^2 - ri^2); %once the two sides converge,
    BMD = area*rho/(2*ro); %store results for this t
    bone_results(t,2) = ro; %in matrix
    bone_results(t,3) = ri;
    bone_results(t,4) = ro_side;
    bone_results(t,5) = ri_side;
    bone_results(t,6) = area;
    bone_results(t,7) = BMD;
    bone_results(t,8) = Z;
    t = t + 1; %move on to the next year

end 














































