function [ weight ] = objective(r)
    Length = 9.14;
    rho = 7860;
    weight = ( 6*r(1)^2 + 4*sqrt(2)*r(2)^2 )*Length*rho*pi;
end