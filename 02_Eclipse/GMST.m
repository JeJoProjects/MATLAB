function [HH,MM,SS]=GMST(yy,mm,dd,hh,min,sec)
    GM=mod(18.697374558+24.06570982441908*...
        (G2JD(yy,mm,dd,hh,min,sec)-G2JD(2000,1,1,12,0,0)),24);
    HH=floor(GM);
    FM=(GM-HH)*60.0;
    MM=floor(FM);
    SS=(FM-MM)*60.0;
end