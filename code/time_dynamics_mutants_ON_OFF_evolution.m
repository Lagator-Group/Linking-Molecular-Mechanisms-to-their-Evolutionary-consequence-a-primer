%% Time dynamics of a given mutant
% Compute temporal dynamics for ON->OFF
function [T, Y] = time_dynamics_mutants_ON_OFF_evolution(weights, values, plt)

energyRNAP = weights{1};
energyRNAPWT = weights{2};
energyCI = weights{3};
energyCIWT = weights{4};
V1new = weights{5};
V2new = weights{6};


RNAP = values(1);
omega0 = values(2);
omega30 = values(3);

tau1 = values(4);
tau2 = values(5);
tau0 = values(6);

R1 = values(7);
R2 = values(8);
n = values(9);

tmax = values(10);
tspan = 0:1e-0:tmax;

slideCI = values(11);

% Prepare parameter values
CIend = R1*tau1; % CI conc. for fully OFF, i.e. CI = highest
omega = omega0/CIend^2;
omega3 = omega30/CIend^4;
CI0 = 0;


[weightOFFCI, weightOFFBoth, weightON, weightONwt] = computeWeightsExpression(energyRNAP, energyRNAPWT, energyCI, energyCIWT, V1new, V2new, RNAP, omega, CIend);


PonWT = sum(weightONwt,1)./(1+sum(weightONwt,1));
Pon = sum(weightON(:,1),1)./(1+sum(weightON(:,1),1));
Poff = sum(weightON(:,2),1)./(sum(weightON(:,2),1)+sum(weightOFFBoth(:,2),1)+mean(weightOFFCI(:,2),1));

if abs(PonWT(1)-PonWT(2))>1e-10
    error 'WT on value wrong?'
end
PonWT = PonWT(1);


YFPstart = Pon/PonWT;
y0 = [CI0 YFPstart];

t0 = 0;
val = [R1 tau1 R2 tau2 RNAP omega omega3 tau0 n t0];
%options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Go by smaller steps so you don't have to simulate all of it
% First define properties to get OFF and ON values

% On and Off YFP values
Yon = Pon/PonWT;
Yoff = Poff/PonWT;

PonWT = PonWT(1);
    


    
Y = [];
T = [];
c = 0;
a = 0;
while c==0 && a<1000/tmax
    a = a+1;
    [TT,YY] = ode45(@(t,y) model_on_off_coop_evolution(t, y, val, energyRNAP, energyRNAPWT, energyCI, energyCIWT, slideCI, V1new, V2new), tspan, y0);
    if a == 1
        T = TT;
    else
        T = [T; T(end)+TT];
    end
    Y = [Y; YY];
    % Time for half maximum for ON -> OFF
    v = Yoff + (Yon - Yoff)/2;
    ind1 = find(Y(:,2) < v);
    if size(ind1,1) == 0
        y0 = [Y(end,1) Y(end,2)];
        t0 = T(end);
        val = [R1 tau1 R2 tau2 RNAP omega omega3 tau0 n t0];
    elseif size(ind1,1) > 0
        c = 1;
    end
end

y0 = [Y(end,1) Y(end,2)];
t0 = T(end);
val = [R1 tau1 R2 tau2 RNAP omega omega3 tau0 n t0];
[TT,YY] = ode45(@(t,y) model_on_off_coop_evolution(t, y, val, energyRNAP, energyRNAPWT, energyCI, energyCIWT, slideCI, V1new, V2new), tspan, y0);
T = [T; T(end)+TT];
Y = [Y; YY];


if sum(find(Y(:,1)/CIend < -1e-3))>0
    error 'CI < 0'
elseif sum(find(Y(:,2) < -1e-3))>0
    error 'YFP < 0'
end


if plt == 1
    figure
    plot(T,Y(:,2),'LineWidth', 2)
    xlabel 'Time [min]'
    ylabel 'YFP [unit of WT]'
    grid on
    box
    
elseif plt == 0
    
else
    error 'Plot or not to plot?'
end



