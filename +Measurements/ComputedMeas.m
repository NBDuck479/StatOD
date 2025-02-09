function [compMeas] = ComputedMeas(yHistRef, Htilde, xMinus, statNumOb, i)

refRangeMeas     = yHistRef.Range(i,statNumOb); 
refRangeRateMeas = yHistRef.RangeRate(i,statNumOb);

refMeas = [refRangeMeas; refRangeRateMeas];

% measurement associated with deivation state
devMeas = Htilde * xMinus; 

% Add to the ref measurement and see what happens
compMeas = rmmissing(refMeas) + devMeas; 

end