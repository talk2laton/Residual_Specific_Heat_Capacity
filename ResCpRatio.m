function Cp_R = ResCpRatio(Pr, Tr)
    % ResCpRatio.m
    % This script compute the ratio of the residual specific heat capacity to
    % the universal gas constant using Hall and yarborough correlation for
    % compresibility factor. 
    %
    % Pr = pseudoreduced pressure;
    % Tr = pseudoreduced temperature;
    % Written by Lateef Kareem - 10/09/2022
    % All equations programed are from the paper in the link below
    % https://www.researchgate.net/profile/Lateef_Kareem/publication/235961260_Specific_Heat_Capacity_of_Natural_Gas_Expressed_as_a_Function_of_Its_Specific_gravity_and_Temperature/links/57dfabfa08ae5292a37f7fd3/Specific-Heat-Capacity-of-Natural-Gas-Expressed-as-a-Function-of-Its-Specific-gravity-and-Temperature.pdf
    if(Pr == 0)
        Cp_R = 0;
        return;
    end
    t = 1/Tr;
    fun = @(p)d2Zdt2_Pr(p, t);
    Cp_R = -t^2*integral(@(P) arrayfun(@(p)fun(p), P), 0, Pr);
end
function v = d2Zdt2_Pr(Pr, t)
    [A, dAdt, d2Adt2, y, dydt, d2ydt2] = computePara(Pr, t);
    v = ((d2Adt2*y^2 + 2*A*dydt^2) - y*(A*d2ydt2 + 2*dAdt*dydt))/y^3;
end
function [A, dAdt, d2Adt2, y, dydt, d2ydt2] = computePara(Pr, t)
    %functions of t
    Afun = @(t) 0.06125 * t * exp(-1.2*(1-t)^2);
    Bfun = @(t) 14.76*t - 9.76*t^2 + 4.58*t^3;
    Cfun = @(t) 90.7*t - 242.2*t^2 + 42.4*t^3;
    Dfun = @(t) 2.18 + 2.82*t;
    A = Afun(t); [dAdt, d2Adt2] = Der(Afun, t);
    B = Bfun(t); [dBdt, d2Bdt2] = Der(Bfun, t);
    C = Cfun(t); [dCdt, d2Cdt2] = Der(Cfun, t);
    D = Dfun(t); [dDdt, d2Ddt2] = Der(Dfun, t);
    f = @(y) -A*Pr + y*(1 + y + y^2 - y^3)/(1-y)^3 - B*y^2 + C*y^D;
    v = 1;
    if(t > 0.8)
        if(0.8 < Pr && Pr < 4)
            v = 2;
        elseif Pr>13
            v = 0.5;
        end
    end
    y = fzero(f,v*A*Pr);
    Num  = dAdt*Pr + dBdt*y^2 - dCdt*y^D - dDdt*C*y^D*log(y);
    Den  = (1+4*y + 4*y^2 - 4*y^3 + y^4)/(1-y)^4 - 2*B*y + C*D*y^(D-1);
    dydt = Num/Den;
    Num  = d2Adt2*Pr + d2Bdt2*y^2 - d2Cdt2*y^D - d2Ddt2*C*y^D*log(y) - ...
        (dDdt*log(y))^2*C*y^D - 2*(-2*dBdt*y + dCdt*D*y^(D-1) + ...
        dDdt*C*y^(D-1)*(D*log(y) + 1))*dydt - ((8 + 20*y - 4*y^2)/...
        (1-y)^5 - 2*B + C*D*(D - 1)*y^(D-2))*dydt^2;
    d2ydt2 = Num/Den;
end
function [f, s] = Der(fun, t)
    dt = 1e-6*(abs(t)+1);
    fmdt = fun(t-dt); fpdt = fun(t+dt); ft = fun(t);
    f = (fpdt - fmdt)/(2*dt);
    s = (fmdt - 2*ft + fpdt)/(dt^2);
end