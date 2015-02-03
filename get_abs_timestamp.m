
function [mattime,matstring]=get_abs_timestamp(tstring)
% Daniel Buscombe August 2011
% make a matlab time vector out of a aqa filename
% version 1.0
y=tstring(1:4);
m=tstring(5:6);
d=tstring(7:8);
h=tstring(9:10);
mm=tstring(11:12);
ss=tstring(13:14);
mattime=datenum([y,'.',m,'.',d,' ',h,':',mm,':',ss],'yyyy.mm.dd HH:MM:SS');
matstring=datestr(mattime);
