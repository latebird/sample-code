function [transitionRate,DecayRate] = ...
    initialization(kRuhop,kRuDecay,kOshop,kOsDecay,kOsTrap,kOsUntrap)

    transitionRate = zeros(2,2);
    transitionRate(1,1) = kRuhop;
    transitionRate(1,2) = kOsTrap;
    transitionRate(2,1) = kOsUntrap;
    transitionRate(2,2) = kOshop;
    
    DecayRate = zeros(1,2);
    DecayRate(1) = kRuDecay; 
    DecayRate(2) = kOsDecay;
    

end 
