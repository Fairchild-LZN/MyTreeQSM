function [Tmin,Tsec] = sec2min(T)
% 将只有秒数，转换为几分几秒

% Transforms the given number of seconds into minutes and residual seconds

Tmin = floor(T/60);
Tsec = round((T-Tmin*60)*10)/10;