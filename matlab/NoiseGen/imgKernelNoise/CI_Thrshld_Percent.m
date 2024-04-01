function [upperMask, lowerMask, upperMskCI,lowerMskCI] = CI_Thrshld_Percent(ciIn,thrld)

pUpper=thrld;
pLower=100-thrld;

PU = prctile(ciIn(:),pUpper);
PL = prctile(ciIn(:),pLower);

upperMask = ciIn >= PU;
upperMskCI=upperMask.*ciIn;
lowerMask = ciIn <= PL;
lowerMskCI=lowerMask.*ciIn;

end