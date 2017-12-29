function [ y] = SigShifted( t,tau,Deltaf,CA)

global If;

y=CAshifted(t,tau,Deltaf,CA)'.*CRshifted(t,tau,Deltaf);

end

