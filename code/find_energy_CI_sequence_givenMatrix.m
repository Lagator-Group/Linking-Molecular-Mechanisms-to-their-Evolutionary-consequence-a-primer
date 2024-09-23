%% Find values of energy for given sequence
% This function finds CI binding energy for a given sequence and a given
% energy matrix
function energy = find_energy_CI_sequence_givenMatrix(sequenceVector, plt, matCI_0)


if size(sequenceVector,2) ~= 4
    error 'Sequence not right size!'
end



%% Find WT sequence

load('sequences_mutants.mat')

% The WT sequence is the 17th sequence
seqWT = sequences(17,:);

%% Transform sequence: A = 0, C = 1, G = 2, T = 3

% WT
kA = strfind(seqWT{1},'A');
kC = strfind(seqWT{1},'C');
kG = strfind(seqWT{1},'G');
kT = strfind(seqWT{1},'T');

sequenceVectorWT = zeros(length(seqWT{1}),4);
sequenceVectorWT(kA,1) = 1;
sequenceVectorWT(kC,2) = 1;
sequenceVectorWT(kG,3) = 1;
sequenceVectorWT(kT,4) = 1;




%% Fix the energy matrix so WT has energy 0

matCI = matCI_0;
m = abs(max(abs(matCI(:))));
matCI = matCI / m;

    
if plt == 1
    figure(102)
    imagesc(matCI')
    colorbar
    axis equal
    ylim([1/2 9/2])
    
    
    for j = 1:size(matCI,1)
        k1 = find(sequenceVectorWT(j+11-1,:) == 1);
        k2 = find(sequenceVectorWT(j+35-1,:) == 1);
        if length(k1)>1
            error 'error'
        elseif length(k2)>1
            error 'error'
        end
        hold on
        plot(j,k1,'rx','MarkerSize',20,'LineWidth',3)
        plot(j,k2,'ro','MarkerSize',20,'LineWidth',3)
        hold off
    end
elseif plt == 0
    
else
    error 'Plot or not to plot?'
end


%% Compute energy for that mutant

matrix = matCI;

% % Old way by sliding
energy = nan(size(sequenceVector,1)-16,1);

for k = 1:size(sequenceVector,1)-16
    
    vec = sequenceVector(k:k+16,:);
    energy(k,1) = sum(sum(matrix.*vec));
    
end
