function f = minjerk(varargin)
% min-jerk ramp from 0 to 1 with length length(t), start time t_start, rise
% time t_rise.

% sz = size(t);
% t = t(:);
% 
% t_up = bsxfun(@and,t>t_start,t<=(t_start+t_rise));
% 
% f = 0*t;
% f(t<=t_start) = 0;
% f(t>(t_start+t_rise)) = 1;
% 
% f(t_up) = 10*((t(t_up)-t_start)/t_rise).^3 - 15*((t(t_up)-t_start)/t_rise).^4 + 6*((t(t_up)-t_start)/t_rise).^5;
% f = reshape(f,sz);

if nargin<1
    error('minjerk needs at least 1 argument');
end

t = varargin{1};
if nargin<2
    t_start = t(1);
    t_rise = t(end)-t(1);
    t_duration = Inf;
end

t_start = varargin{2}(1);
if numel(varargin{2})<3
    t_duration = Inf;
    if numel(varargin{2})<2
        t_rise = t(end)-t(1);
    else
        t_rise = varargin{2}(2);
    end
else
    t_start = varargin{2}(1);
    t_rise = varargin{2}(2);
    t_duration = varargin{2}(3);
end

% t_up = bsxfun(@and,t>=t_start,t<(t_start+t_rise));
% t_mid = bsxfun(@and,t>=(t_start+t_rise),t<(t_start+t_rise+t_duration));
% t_down = bsxfun(@and,t>=(t_start+t_rise+t_duration),t<(t_start+2*t_rise+t_duration));
% 
f = 0*t;
% f(t<t_start) = 0;
% f(t>t_down) = 0;
% f(t_mid) = 1;

% f(t_up) = 10*((t(t_up)-t_start)/t_rise).^3 - 15*((t(t_up)-t_start)/t_rise).^4 + 6*((t(t_up)-t_start)/t_rise).^5;
% f(t_down) = 1-(10*((t(t_down)-(t_start+t_rise+t_duration))/t_rise).^3 - 15*((t(t_down)-(t_start+t_rise+t_duration))/t_rise).^4 + 6*((t(t_down)-(t_start+t_rise+t_duration))/t_rise).^5);

f(t>t_start) = 10*((t(t>t_start)-t_start)/t_rise).^3 - 15*((t(t>t_start)-t_start)/t_rise).^4 + 6*((t(t>t_start)-t_start)/t_rise).^5;
f(t>(t_start+t_rise)) = 1;
f(t>(t_start+t_rise+t_duration)) = 1-(10*((t(t>(t_start+t_rise+t_duration))-(t_start+t_rise+t_duration))/t_rise).^3 - 15*((t(t>(t_start+t_rise+t_duration))-(t_start+t_rise+t_duration))/t_rise).^4 + 6*((t(t>(t_start+t_rise+t_duration))-(t_start+t_rise+t_duration))/t_rise).^5);
f(t>(t_start+2*t_rise+t_duration)) = 0;

end