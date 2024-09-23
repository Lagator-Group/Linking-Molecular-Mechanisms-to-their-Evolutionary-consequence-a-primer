%% Finding dynamical phenotypes
function [M, m, slope1, slope2, half1, half2, www1] = evolutionSequencesMinFunc(seqInd, matRNAP, matCI, slideRNAP, slideCI, ver00, timeEvol)




%% Load the files

load('filesForEvolution_190607.mat')

% Correct the WT sequence
wtVec=[zeros(1,4);sequenceVectorWT;zeros(12,4)];

wtVec(1,1) = 1;
for j = 1:12
    wtVec(68+j,1) = 1;
end
sequenceVectorWT = wtVec;

%%

ind = ([1:80]*4-4)+seqInd';
seq = zeros(4,80);
seq(round(ind)) = 1;
sequenceVector = seq';

% Compute the energies of the mutant
for j = 1:5 % over different spacers
    matRNAP00 = matRNAP{j};
    matRNAP0 = [zeros(4,4); matRNAP00; zeros(22,4)];
    energyRNAP0{j} = find_energy_RNAP_sequence_givenMatrix(sequenceVector, 0, matRNAP0);
    energyRNAPWT0{j} = find_energy_RNAP_sequence_givenMatrix(sequenceVectorWT, 0, matRNAP0);
end
energyCI0 = find_energy_CI_sequence_givenMatrix(sequenceVector, 0, matCI);
energyCIWT0 = find_energy_CI_sequence_givenMatrix(sequenceVectorWT, 0, matCI); % WT energy
if slideCI == 0
    energyCI0 = [energyCI0(11); energyCI0(35)];
    energyCIWT0 = [energyCIWT0(11); energyCIWT0(35)];
end

deltaSpacer = alpha0*deltaSpacer0;

for j = 1:5
    energyRNAP{j} = alpha0*energyRNAP0{j};
    energyRNAPWT{j} = alpha0*energyRNAPWT0{j};
end
energyCI = iota0*energyCI0;
energyCIWT = iota0*energyCIWT0;

% Compute the mutant properties
% On, Off values get from analytical result
tau1 = 60;
R1 = 1;
CI0 = 0;
CIend = R1*tau1; % CI conc. for fully OFF, i.e. CI = highest
omega = omega0/CIend^2;
omega3 = omega30/CIend^4;

[weightOFFCI, weightOFFBoth, weightON, weightONwt, www1] = computeWeightsExpression(energyRNAP, energyRNAPWT, energyCI, energyCIWT, V1new, V2new, RNAP, omega, CIend);


PonWT = sum(weightONwt,1)./(1+sum(weightONwt,1));
Pon = sum(weightON(:,1),1)./(1+sum(weightON(:,1),1));
Poff = sum(weightON(:,2),1)./(sum(weightON(:,2),1)+sum(weightOFFBoth(:,2),1)+mean(weightOFFCI(:,2),1));

% On and Off YFP values
Yon = Pon/PonWT(1);
Yoff = Poff/PonWT(1);

if abs(PonWT(1)-PonWT(2))>1e-10
    error 'WT ON value wrong?'
end
PonWT = PonWT(1);


if timeEvol == 1
    % Do time evolution
    % Time parameter values
    tau2 = tau1;
    tau0a = 80;
    tau0b = 15;
    R2 = 1/(tau2*PonWT);
    
    n = 4;
    beta = 2;
    tmax = 200;
    plt1 = 0;
    
    weights{1} = energyRNAP;
    weights{2} = energyRNAPWT;
    weights{3} = energyCI;
    weights{4} = energyCIWT;
    weights{5} = V1new;
    weights{6} = V2new;
    
    % Compute results for ON -> OFF
    values = [RNAP omega0 omega30 tau1 tau2 tau0a R1 R2 n tmax slideCI];
    % Time evolve the mutant
    [T1, Y1] = time_dynamics_mutants_ON_OFF_evolution(weights, values, plt1);
    TT1 = T1;
    YY1 = Y1;
    
    % Compute results for OFF -> ON
    values = [RNAP omega0 omega30 tau1 tau2 tau0b R1 R2 beta tmax slideCI];
    % Time evolve the mutant
    [T2, Y2] = time_dynamics_mutants_OFF_ON_evolution(weights, values, plt1);
    TT2 = T2;
    YY2 = Y2;
    
    
    
    % Time for half maximum for ON -> OFF
    v = Yoff + (Yon - Yoff)/2;
    ind1 = find(Y1(:,2) < v,1);
    if abs(Yoff-Yon)<1e-14
        half1 = 0;
        half2 = 0;
        slope1 = 0;
        slope2 = 0;
    else
        if size(ind1,1) == 0
            keyboard
            error 'Time too short!'
        end
        T0 = T1(ind1(1)-1) + (T1(ind1(1))-T1(ind1(1)-1)) * (v-Y1(ind1(1)-1,2)) / (Y1(ind1(1),2)-Y1(ind1(1)-1,2));
        half1 = T0;
        %half1(k,1) = T1(ind1(1));
        
        % Time for half maximum for OFF -> ON
        v = Yoff + (Yon - Yoff)/2;
        ind2 = find(Y2(:,2) > v,1);
        if size(ind2,1) == 0
            error 'Time too short!'
        end
        T0 = T2(ind2(1)-1) + (T2(ind2(1))-T2(ind2(1)-1)) * (v-Y2(ind2(1)-1,2)) / (Y2(ind2(1),2)-Y2(ind2(1)-1,2));
        half2 = T0;
        %half2(k,1) = T2(ind2(1));
        
        % Slope at half maximum for ON -> OFF
        derY = diff(Y1(:,2));
        derT = diff(T1(:,1));
        if derT(ind1(1)) == 0
            ind1(1) = ind1(1)+1;
        end
        %slope1(k,1) = der(ind1(1))/( T1(ind1(1)) - T1(ind1(1)-1) );
        slope1 = derY(ind1(1))/derT(ind1(1));
        
        % Slope at half maximum for OFF -> ON
        derY = diff(Y2(:,2));
        derT = diff(T2(:,1));
        if derT(ind2(1)) == 0
            ind2(1) = ind2(1)+1;
        end
        %slope2(k,1) = der(ind2(1))/(T2(ind2(1)) - T2(ind2(1)-1));
        slope2 = derY(ind2(1))/derT(ind2(1));
    end
else
    slope1 = 0;
    slope2 =0;
    half1 = 0;
    half2 = 0;
end

% Compute the time dynamic properties
% Max and min value
M = Yon;
m = Yoff;




