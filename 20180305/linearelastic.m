E = 1000000;
nu = 0.2;

K = E / 3 / (1 - 2*nu);
G = E / 2 / (1 + nu);

for x=1:6
    for y=1:6
        if x<=3
            if y<=3
                if x==y
                    De(x,y)=K+4/3*G; 
                else
                    De(x,y)=K-2/3*G;
                end
            end
        end
        if x>3
            if y>3
                if x==y
                    De(x,y)=2*G;
                end
            end
        end
    end 
end

dStrain = [0,-2e-08,0,0,0,0]';
stress = De * dStrain

dStrain = [7.40741e-16,-2e-08,0,-5.42101e-25,0,0]';
stress = stress + De * dStrain

dStrain = [2.22222e-15,-2e-08,0,8.13152e-25,0,0]';
stress = stress + De * dStrain
