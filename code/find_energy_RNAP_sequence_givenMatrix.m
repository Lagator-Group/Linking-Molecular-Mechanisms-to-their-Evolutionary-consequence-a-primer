%% Find values of energy for given sequence
% This function finds RNAP binding energy for a given sequence and a given
% energy matrix
% Energy matrix can vary from the 'original', see Figs. 4-6 in the paper
function energy = find_energy_RNAP_sequence_givenMatrix(sequenceVector, plt, mat0)


if size(sequenceVector,2) ~= 4
    keyboard
    error 'Sequence not right size!'
end


%% Fix the energy matrix so WT has energy 0

mat = mat0;



%% Compute energy for that mutant

matrix = mat;
L = size(matrix,1)-1;

% % Old way by sliding
energy = nan(size(sequenceVector,1) - L,1);

for k = 1:size(sequenceVector,1)- L
    
    vec = sequenceVector(k:k+L,:);
    energy(k,1) = sum(sum(matrix.*vec));
    
end

