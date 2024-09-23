%% Finding a promoter with optimal values of ON~1 and OFF~0
% We can it 'evolution' as it evolves through 
function [] = computeEvolutionofPromoter_190420(count, Tinit, smallValue)


%%  Load WT sequence

load('doubleMutants_Pl_short_v2.mat');




%% Load the energy matrices

plt = 0;
% Load RNAP energy matrix
load('35.txt')
load('10.txt')
mat0 = [X_35; zeros(8,4); X_10];

% Load CI energy matrix
load('C1_from_DNN.txt')
matCI0 = -C1_from_DNN;

matRNAP0 = nan(size(mat0));
for j = 1:size(mat0,1)
    if strcmp(EM,'New')
        k = find(sequenceVectorWT(j+23-1,:) == 1);
    elseif strcmp(EM,'Old')
        k = find(sequenceVectorWT(j+24-1,:) == 1);
    end
    if length(k)>1
        error 'error'
    end
    matRNAP0(j,:) = mat0(j,:) - mat0(j,k);
end

m = abs(max(abs(matRNAP0(:))));
matRNAP0 = matRNAP0 / m;
mat1a = matRNAP0(1:12,:);
mat1b = matRNAP0(21:32,:);



for j = 1:5 % over different spacers
    matRNAP{j} = [mat1a; zeros(5+j,4); mat1b];
end




matCI = matCI0;
matCInorm = nan(size(matCI0));
for j = 1:size(matCI0,1)
    k = find(sequenceVectorWT(j+11-1,:) == 1);
    if length(k)>1
        error 'error'
    end
    matCInorm(j,:) = matCI(j,:) - matCI(j,k);
end

m = abs(max(abs(matCI(:))));
matCI = matCI / m;



%% Do the simulated annealing "evolution"
tic

ver00 = 1;
slideRNAP = 1;
slideCI = 1;
shrani = 1;

rng('shuffle')
seqInit = randi(4,[67 1]);

kmax = 1e4;
valVec = nan(kmax,1);
Tcount = nan(kmax,1);

[M, m] = evolutionSequencesMinFunc(seqInit, matRNAP, matCI, slideRNAP, slideCI, ver00);
valInit = (M-1)^2 + m^2;

seq = seqInit;
val = valInit;
T = Tinit;

for k = 1:kmax
    if rem(k,100) == 0
        k
        val
        toc
    end
    
    % Update temperature
    T = T/(1+smallValue);
    % Have T=0 for last few iterations
    if k>kmax-100
        T = 0;
    end
    
    
    % Update to a new state
    l1 = randi(length(seq),1);
    l2 = randi(4,1);
    for hh = 1:5
        if seq(l1,1) == l2
            l2 = randi(4,1);
        end
    end
    seqNew = seq;
    seqNew(l1,1) = l2;
    
    % Compute new phenotypes/energy
    [M, m, www1] = evolutionSequencesMinFunc(seqNew, matRNAP, matCI, slideRNAP, slideCI, ver00);
    valNew = (M-1)^2 + m^2;
    dE = (valNew - val);
    
    % Update accordingly
    if exp(-dE/T) > rand
        seq = seqNew;
        val = valNew;
    end
    valVec(k,1) = val;
    Tcount(k,1) = T;
    
    if val<1e-3
        break
    end
end
toc


% Save the results
if shrani == 1
    filename = (['/nfs/scistore12/calingrp/rgrah/evolOfPromoter/evolutionOfSequencesResults_differentSlide_RNAP' num2str(slideRNAP) '_CI' num2str(slideCI) '_newParam' num2str(ver00) '_rnapMat' EM '_K' num2str(count) '_190420.mat']);
    save(filename)
    %, 'A','M','m','slope1','slope2','half1','half2','Awt','Mwt','mwt','slope1wt','slope2wt','half1wt','half2wt','mutSeq0','mutSeq','matRNAPnew','matRNAPnew0','matCInew','matCInew0');
end


