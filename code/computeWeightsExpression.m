%% Function to compute Boltzmann weights of the model
function [weightOFFCI, weightOFFBoth, weightON, weightONwt, www1] = computeWeightsExpression(energyRNAP, energyRNAPWT, energyCI, energyCIWT, V1new, V2new, RNAP, omega, CIend)

%%

% Compute weights - ALL states
for j = 1:5
    energ = energyRNAP{j};
    energWT = energyRNAPWT{j};
        
    v1 = V1new{j};
    v1(v1==0) = nan;
    v2 = V2new{j};
    v2(v2==0) = nan;
    
    ww = nan(size(v1,1),4,2);
    
    weightOFFCIm = 0;
    weightOFFBothm = 0;
    weightONm = 0;
    weightONwtm = 0;
    for k = 1:size(v1,1)
        k;
        clear weiCI weiBoth weiON weiONwt
        ind1 = v1(k,~isnan(v1(k,:)));
        
        %if ismember(1,ind1) && (sum(ind1==2)+sum(ind1==3))==0 % To see if only RNAP(s) is/are present - use for expression state
        checkVal = 1;% Check if CI is bound before RNAP
        a1 = find(ind1==1); % index of RNAP bound
        a2 = find(ind1==2); % index of CI bound
        if a2>a1 % If CI is bound after RNAP, expression is NOT allowed
            checkVal = 0;
        end
        
        % Check if RNAP is bound, or bound before CI
        if ismember(1,ind1) && checkVal==1
            weiON = 1;
            weiONwt = 1;
        else 
            weiON = 0;
            weiONwt = 0;
        end
        
        % Check if only CI is bound
        if ~ismember(1,ind1)
            weiCI = 1;
        else
            weiCI = 0;
        end
        
        % Check if CI is bound after RNAP
        weiBoth = 0;
        if ismember(1,ind1) && ismember(2,ind1) && a2>a1
            weiBoth = 1;
        end
            
        for k1 = 1:length(ind1)
            k2 = ind1(k1);
            ind2 = v2(k,k1); 
            ind2;
%             if ind2==53
%                 Pkeyboard
%             end
            switch k2
                case 1
                    %error 'ON expression is when NO CI is present!'
                    weiBoth = RNAP*exp(-energ(ind2)) .* weiBoth;
                    weiON = RNAP*exp(-energ(ind2)) .* weiON;
                    weiONwt = RNAP*exp(-energWT(ind2)) .* weiONwt;
                case 2
                    weiCI = omega*[0 CIend^2]*exp(-energyCI(ind2)) .* weiCI;
                    weiBoth = omega*[0 CIend^2]*exp(-energyCI(ind2)) .* weiBoth;
                    %error 'Do we include weiON to this part?'
                    weiON = omega*[0 CIend^2]*exp(-energyCI(ind2)) .* weiON;
                    weiONwt = omega*[0 CIend^2]*exp(-energyCIWT(ind2)) .* weiONwt;
                case 3
                    %wei = omega3*CIend^4*exp(-energyCI(ind2)-energyCI(ind2+24)) * wei;
            end
        end
        weightOFFCIm = weightOFFCIm + weiCI;
        weightOFFBothm = weightOFFBothm + weiBoth;
        weightONm = weightONm + weiON;
        weightONwtm = weightONwtm + weiONwt;
        ww(k,1,1:2) = weiCI;
        ww(k,2,1:2) = weiBoth;
        ww(k,3,1:2) = weiON;
        ww(k,4,1:2) = weiONwt;
        
    end
    www1{j} = ww;
    
    weightOFFCI(j,1:2) = weightOFFCIm;
    weightOFFBoth(j,1:2) = weightOFFBothm;
    weightON(j,1:2) = weightONm;
    weightONwt(j,1:2) = weightONwtm;
    
end


