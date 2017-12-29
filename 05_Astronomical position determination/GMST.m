function sd=GMST(yy,mm,dd,hh,mn,ss)
sd=mod(18.697374558+24.06570982441908*(G2JD(yy,mm,dd,hh,mn,ss)-2451545.0),24);
end