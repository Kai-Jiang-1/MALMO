function x = truncnormrnd(a, b, mu, sigma)
% Draw from a normal(mu, sigma) that is restricted to [a, b].
% Uses MATLAB's normcdf / norminv (Statistics toolbox).
%
% Idea: the CDF at a and b is the allowed probability range.
% Draw a uniform number in that range, then invert the CDF.

u = rand(size(mu));
Fa = normcdf(a, mu, sigma);   % P(Y <= a)
Fb = normcdf(b, mu, sigma);   % P(Y <= b)
p  = Fa + u .* (Fb - Fa);     % uniform between Fa and Fb
p  = min(max(p, eps), 1-eps); % keep p inside (0, 1) for norminv
x  = norminv(p, mu, sigma);
x  = min(max(x, a), b);       % keep the draw inside [a, b]
end

