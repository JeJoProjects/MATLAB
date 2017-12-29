function y=SigShifted(t,tau,Deltaf,CA)

global If

y=CAshifted(t,tau,Deltaf,CA)'.*CRshiftet(t,tau,Deltaf);
