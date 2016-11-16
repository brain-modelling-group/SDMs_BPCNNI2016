%% Heun.m
% Heun integrator for SDE in Stratonovich form. 
% MJA version 2011-04-30
% modified by JAR 2012-21-09 (small speedup)
%
% Arguments:
% f, g      Function handles defining the RHS of Stratonovich SODE system: dy = f(t,y).dt + g(t,y).dW
% tspan     Vector of fixed-size time steps for integration: [0:0.01:5] 
% y0        Vector of initial conditions
% display   true or false (progressive plots during integration)
%
% Bugs:
% - Input not yet properly validated
% - Progressive display not yet implemented
%
% References: 
% Mannella, R (2002) Integration of Stochastic Differential Equations on a Computer
% Rumelin, W (1982) Numerical Treatment of Stochastic Differential Equations
%
function sol = Heun(f, g, tspan, y0, display)
    % get dimension of system
    n = numel(y0);
    steps = numel(tspan) - 1;
    % validate some arguments
    if numel(g(tspan(1),y0))~=n || numel(f(tspan(1),y0))~=n
        ex = MException('Heun:InvalidArgs', 'Argument dimensions are inconsistent');
        throw(ex);
    end
    % create array to hold results
    sol.y = [y0, zeros(n, steps)];
    sol.x = tspan;
    sol.solver = 'Heun';
    % integrate
    for k = 2:(steps+1)
        tn = tspan(k-1);
        tm = tspan(k);
        yn = sol.y(:, k-1);
        dt = tm - tn;
        dW = sqrt(dt).*randn(n, 1);
        fn=f(tn, yn);
        gn=g(tn, yn);
        ybar = yn + fn.*dt + gn.*dW;
        sol.y(:, k) = yn + 0.5.*(fn + f(tm, ybar)).*dt + 0.5.*(gn + g(tm, ybar)).*dW;

        %if display
            %TODO progressive display
        %end
    end


end
