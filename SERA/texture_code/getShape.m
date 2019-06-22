function [Shape] = getShape(V,S)
% V vol ume  need to cm^3   S surface   cm^2


% V = V/1e3;
compactness1 = V/sqrt(pi)/sqrt(S^3);
compactness2 = 36*pi*V^2/S^3;
% sphericity = pi^(1/3)*sqrt((6*V)^3)/S;
sphericity = (36*pi*V^2)^(1/3)/S;
sphericalDisp = S / (36*pi*V^2)^(1/3);
Asphericity = (S^3 / V^2 / 36 / pi)^(1/3) - 1; 
SVratio = S/V;
Scircle = (3*V/(4*pi))^(2/3)*4*pi;
Irregularity = S/Scircle;


Shape.compactness1 = compactness1;
Shape.compactness2 = compactness2;
Shape.sphericity = sphericity;
Shape.SVratio = SVratio;
Shape.Irregularity = Irregularity;
Shape.sphericalDisp = sphericalDisp;
Shape.Asphericity = Asphericity;

end