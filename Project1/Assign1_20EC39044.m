%% Computational Neuroscience -- Project 1
%   -- Avi Amalanshu
%   -- 20EC39044
%% Parameters
mu = [1, 0.1, 100];
TSPAN = {[0,100], [0,100], [0,500]};
x_0 = [1, 0];
times = zeros(2,3); % instantiation
%% ODE45 Solution
figure
for i = 1:3
    t0 = tic();
    [t, x] = ode45(                             ...
                    @(t, x) deq(t, x, mu(i)),   ...
                    TSPAN{i},                   ...
                    x_0                         ...
             );
    t1 = toc(t0);
    times(1,i) = t1;
    subplot(2,3,i)
    plot(t, x(:,1))
    if i == 2
        title('Response')
    end
    legend(sprintf('mu = %3.1f', mu(i)))
    xlabel('Time')
    ylabel('State Variable x_1 = y')
    subplot(2,3,i+3)
    plot(x(:,1), x(:,2),'-.')
    if i == 2
        title('Phase-Plane')
    end
    legend(sprintf('mu = %3.1f', mu(i)))
    xlabel('State Variable x_1 = y')
    ylabel('State Variable x_2 = dy/dx')
end
%% ODE15s Solution
for i = 1:3
    t0 = tic();
    [t, x] = ode15s(                            ...
                    @(t, x) deq(t, x, mu(i)),   ...
                    TSPAN{i},                   ...
                    x_0                         ...
             );
    t1 = toc(t0);
    times(2,i) = t1;
end
%% Times
times
%% State Function
function dxdt = deq(t, x, mu)
    dxdt = [mu*x(2); mu*(1-x(1)^2)*x(2) - x(1)/mu];
end
