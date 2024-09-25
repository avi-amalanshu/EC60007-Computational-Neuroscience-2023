%% Computational Neuroscience -- Project 2
%   -- Avi Amalanshu
%   -- 20EC39044
%% Parameters
global G_Ca G_K G_L E_Ca E_K E_L phi V_1 V_2 V_3 V_4 V_5 V_6 C; %I_ext; 

G_Ca = 4.4; 
G_K = 8.0; 
G_L = 2; 
E_Ca = 120; 
E_K = -84; 
E_L = -60;
phi = 0.02; 
V_1 = -1.2; 
V_2 = 18;  
V_3 = 2; 
V_4 = 30; 
V_5 = 2; 
V_6 = 30;
C = 20;
I_ext = 0;
%% Question 2
% Null-cline plots and quiver plot
figure;
fimplicit(@(V, w) dVdt(V, w/100, I_ext), [-100 100 0 100]);
hold on
fimplicit(@(V, w) dwdt(V, w/100), [-100 100 0 100]);
xlabel('V');
ylabel('w');
title('Nullclines of V and w');
hold on;
V_range = linspace(-100,100,100);
w_range = linspace(0,1,100);
[V_quiv, w_quiv] = meshgrid(V_range, w_range);
V_dot_vec = dVdt(V_quiv, w_quiv, I_ext);
w_dot_vec = dwdt(V_quiv, w_quiv);
quiver(V_quiv, w_quiv*100, V_dot_vec, w_dot_vec*100, 4)
legend('V-nullcline', 'w-nullcline', 'quiver plot');
grid on;

% Finding the equilibrium point (method 1, numerical)
eq_1 = fsolve(@(X)[dVdt(X(1), X(2), I_ext); dwdt(X(1),X(2))], [-60, 0.01]);

% Finding the equilibrium point (method 2, setting off system)
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,2000], [-60, 0.01]);
eq_2 = X_(end,:);
%% Question 3
% Finding the Jacobian
syms V_ w_
J = jacobian([dVdt(V_,w_, I_ext), dwdt(V_,w_)],[V_, w_]);
J_mat = matlabFunction(J);
J_0 = feval(J_mat, eq_1(1), eq_1(2));

% Finding the Eigenvalues
E = eig(J_0);
%% Question 5
% Effect of Phi: Morris Lecar ODE solution
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [0, eq_1(2)]); % time step 0-300, init voltage 0mV, init w = eqb val
figure()
sgtitle('Morris-Lecar Equations: Phase-Plane Plots')
plot(X_(:,1), X_(:,2), LineWidth=3)
hold on
phi = 0.04;
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [0, eq_1(2)]); % time step 0-300, init voltage 0mV, init w = eqb val
plot(X_(:,1), X_(:,2), LineWidth=3)
hold on
phi = 0.01;
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [0, eq_1(2)]); % time step 0-300, init voltage 0mV, init w = eqb val
plot(X_(:,1), X_(:,2), LineWidth=3)
hold on
fimplicit(@(V, w) dVdt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) dwdt(V, w), [-100 100 0 1]);
hold off
xlabel('V');
ylabel('w');
legend('Sys. Traj. \phi = 0.02', 'Sys. Traj. \phi = 0.04', 'Sys. Traj. \phi = 0.01', 'V Null-Cline', 'w Null-Cline')

phi = 0.02;
%% Question 6
% Effect of initial V on peak voltage (observing thresholding)
Vs = -30:0.01:15;
L = length(Vs);
Xs = cell(L);
for i = 1:L
    [~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [Vs(i), eq_1(2)]); % time step 0-300, init voltage 0mV, init w = eqb val
    Xs{i} = X_;
end
max_Vs = zeros(1, L);
for i = 1:length(Vs)
    max_Vs(i) = max(Xs{i}(:,1));
end
figure()
plot(Vs, max_Vs)
title('Voltage peak vs starting voltage')
xlabel('V (mV)')
ylabel('V_{max} (mV)')
% Trajectories for some values of sub-(and super-)threshold initial V.
figure()
sgtitle('Trajectories for some initial values of V')
subplot(121)
for i = 1:15
    plot(Xs{50*i}(:,1), Xs{50*i}(:,2), LineWidth=2)
    hold on
end
fimplicit(@(V, w) dVdt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) dwdt(V, w), [-100 100 0 1]);
hold off
title('Subthreshold')
xlabel('V (mV)')
ylabel('w')
subplot(122)
for i = 1:15
    plot(Xs{floor(3*L/4) + 50*i}(:,1), Xs{floor(3*L/4) + 50*i}(:,2), LineWidth=2)
    hold on
end
fimplicit(@(V, w) dVdt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) dwdt(V, w), [-100 100 0 1]);
hold off
xlabel('V (mV)')
ylabel('w')
title('Superthreshold')
%% Question 7
% Demonstrating the shift in nullclines (and therefore equilibrium points)
figure;
fimplicit(@(V, w) dVdt(V, w/100, I_ext), [-100 100 0 100]);
hold on
I_ext = 86;
fimplicit(@(V, w) dVdt(V, w/100, I_ext), [-100 100 0 100]);
hold on
fimplicit(@(V, w) dwdt(V, w/100), [-100 100 0 100]);
legend('Old V null-cline', 'New V null-cline', 'w null-cline')
xlabel('V');
ylabel('w');

% New equilibrium point
eq_new = fsolve(@(X)[dVdt(X(1), X(2), I_ext); dwdt(X(1),X(2))], [-30, 0.1]);

figure;
fimplicit(@(V, w) dVdt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) dwdt(V, w), [-100 100 0 1]);
hold on
% Trial 1 for new trajectories (starting sys from old equilibrium point)
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [eq_1(1), eq_1(2)]);
plot(X_(:,1), X_(:,2), LineWidth=2)
hold on
% Trial 2 for new trajectories (starting sys at new equilibrium point)
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [eq_2(1), eq_2(2)]);
plot(X_(:,1), X_(:,2), LineWidth=2)
hold on
% Trial 3 for new trajectories (starting sys off new equilibrium point)
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [-27.9, 0.17]);
plot(X_(:,1), X_(:,2), LineWidth=2)
hold off
xlabel('V (mV)')
ylabel('w')
legend('V null-cline', 'w null-cline', 'Trial 1', 'Trial 2', 'Trial 3')

% Finding the new Jacobian
J_new = jacobian([dVdt(V_,w_, I_ext), dwdt(V_,w_)],[V_, w_]);
J_mat_new = matlabFunction(J_new);
J_new_0 = feval(J_mat_new, eq_new(1), eq_new(2));

% Finding the new Eigenvalues
E_new = eig(J_new_0);
%% Question 8
% Running the system backwards in time
figure;
fimplicit(@(V, w) -dVdt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) -dwdt(V, w), [-100 100 0 1]);
hold on
[~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,-300], [-27.9, 0.17]); % Running [0, 300] backward in time
plot(X_(:,1), X_(:,2), LineWidth=3)
hold off
xlabel('V (mV)')
ylabel('w')
legend('V null-cline', 'w null-cline', 'UPO')
%% Question 9
I_ext = 80;
eq_new_80 = fsolve(@(X)[dVdt(X(1), X(2), I_ext); dwdt(X(1),X(2))], [-30, 0.1]);
J_new_80 = jacobian([dVdt(V_,w_, I_ext), dwdt(V_,w_)],[V_, w_]);
J_mat_new_80 = matlabFunction(J_new_80);
J_new_80_0 = feval(J_mat_new_80, eq_new_80(1), eq_new_80(2));
E_new_80 = eig(J_new_80_0);

I_ext = 90;
eq_new_90 = fsolve(@(X)[dVdt(X(1), X(2), I_ext); dwdt(X(1),X(2))], [-30, 0.1]);
J_new_90 = jacobian([dVdt(V_,w_, I_ext), dwdt(V_,w_)],[V_, w_]);
J_mat_new_90 = matlabFunction(J_new_90);
J_new_90_0 = feval(J_mat_new_90, eq_new_90(1), eq_new_90(2));
E_new_90 = eig(J_new_90_0);

I = 80:1:100;
firing_rates = zeros(1,length(I));
for i = 1:length(I)
    I_ext = I(i);
    eq_new_I = fsolve(@(X)[dVdt(X(1), X(2), I_ext); dwdt(X(1),X(2))], [-30, 0.1]);
    [~, X_] = ode45(@(t, X)mle(t, X, I_ext), [0,300], [1.05*eq_new_I(1), 1.05*eq_new_I(2)]);
    v_trace = X_(:,1);
    % Count number of times the membrane voltage rises above 0
    [peak_values, peak_times] = findpeaks(v_trace, 'MinPeakHeight', 0);
    % Calculate firing rate
    num_spikes = length(peak_times); % Number of spikes
    firing_rates(i) = 1000 * num_spikes / length(v_trace); % Firing rate in Hz
end
I_ext = 0;
figure;
plot(I,firing_rates);
title('Effect of I_{ext} on the MLE firing rate')
xlabel('I_{ext} (\mu A/cm^2)')
ylabel('Firing rate (Hz)')
%% Question 10
% New parameters
global G_Ca_alt G_K_alt G_L_alt E_Ca_alt E_K_alt E_L_alt phi_alt V_1_alt V_2_alt V_3_alt V_4_alt V_5_alt V_6_alt C_alt;

G_Ca_alt=4; 
G_K_alt=8.0; 
G_L_alt=2; 
E_Ca_alt=120; 
E_K_alt=-84; 
E_L_alt=-60;
phi_alt=0.0667; 
V_1_alt=-1.2; 
V_2_alt=18;  
V_3_alt=12; 
V_4_alt=17.4; 
V_5_alt=12; 
V_6_alt=17.4;
C_alt=20;
I_ext=30;

% Finding the three equilibrium points for these parameters
eq_alt_params_1 = fsolve(@(X)[dVdt_alt(X(1), X(2), I_ext); dwdt_alt(X(1),X(2))], [-40, 0]);
eq_alt_params_2 = fsolve(@(X)[dVdt_alt(X(1), X(2), I_ext); dwdt_alt(X(1),X(2))], [-19, 0.02]);
eq_alt_params_3 = fsolve(@(X)[dVdt_alt(X(1), X(2), I_ext); dwdt_alt(X(1),X(2))], [4, 0.2]);

% Evaluating the Jacobians at the new equilibrium points
J_alt = jacobian([dVdt_alt(V_,w_, I_ext), dwdt_alt(V_,w_)],[V_, w_]);
J_mat_alt = matlabFunction(J_alt);
J_alt_eq_1 = feval(J_mat_alt, eq_alt_params_1(1), eq_alt_params_1(2));
J_alt_eq_2 = feval(J_mat_alt, eq_alt_params_2(1), eq_alt_params_2(2));
J_alt_eq_3 = feval(J_mat_alt, eq_alt_params_3(1), eq_alt_params_3(2));

[EV_1, D_1] = eig(J_alt_eq_1);
[EV_2, D_2] = eig(J_alt_eq_2);
[EV_3, D_3] = eig(J_alt_eq_3);

E_alt_params_1 = [D_1(1,1); D_1(2,2)];
E_alt_params_2 = [D_2(1,1); D_2(2,2)];
E_alt_params_3 = [D_3(1,1); D_3(2,2)];

% Plotting the null clines and manifolds
figure;
fimplicit(@(V, w) dVdt_alt(V, w, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, w) dwdt_alt(V, w), [-100 100 0 1]);
hold on
% Unstable manifold perturbation in one direction
[~, X_] = ode45(@(t, X)mle_alt(t, X, I_ext), ...
    [0,2000], ...
    [eq_alt_params_2(1)+0.01,eq_alt_params_2(2)+0.01*EV_2(2,1)]);
plot(X_(:,1), X_(:,2), 'k', LineWidth=2)
hold on
% Reversing time and perturbing for stable manifold
[~, X_] = ode45(@(t, X)mle_alt(t, X, I_ext), ...
    [0,-2000], ...
    [eq_alt_params_2(1)+0.0001,eq_alt_params_2(2)+0.0001*EV_2(2,2)]);
plot(X_(:,1), X_(:,2), 'c', LineWidth=2)
hold on
% Unstable manifold perturbation in the other direction
[~, X_] = ode45(@(t, X)mle_alt(t, X, I_ext), ...
    [0,2000], ...
    [eq_alt_params_2(1)-0.01,eq_alt_params_2(2)-0.01*EV_2(2,1)]);
plot(X_(:,1), X_(:,2), 'k', LineWidth=2)
hold on
% Reversing time and perturbing for unstable manifold
[~, X_] = ode45(@(t, X)mle_alt(t, X, I_ext), ...
    [0,-165], ...
    [eq_alt_params_2(1)-0.0001,eq_alt_params_2(2)-0.0001*EV_2(2,2)]);
plot(X_(:,1), X_(:,2), 'c', LineWidth=2)
hold off
xlabel('V (mV)')
ylabel('w')
title('Null clines for the new parameter set and I_{ext} = 30\mu A/cm^2')
legend('V null-cline', 'w null-cline', 'unstable manifold', 'stable manifold')
%% Question 11
I = [30:39, 39.01:0.01:40, 41:45];
Es_alt_params_1 = zeros(2,length(I));
Es_alt_params_2 = zeros(2,length(I));
Es_alt_params_3 = zeros(2,length(I));
firing_rates = zeros(1, length(I));
for i = 1:length(I)
    % Finding the equilibrium points
    eq_1_ = fsolve(@(X)[dVdt_alt(X(1), X(2), I(i)); dwdt_alt(X(1),X(2))], [-40, 0]);
    eq_2_ = fsolve(@(X)[dVdt_alt(X(1), X(2), I(i)); dwdt_alt(X(1),X(2))], [-19, 0.02]);
    eq_3_ = fsolve(@(X)[dVdt_alt(X(1), X(2), I(i)); dwdt_alt(X(1),X(2))], [4, 0.2]);
    J_ = jacobian([dVdt_alt(V_,w_, I), dwdt_alt(V_,w_)],[V_, w_]);
    J_mat_ = matlabFunction(J_);
    % Finding the Jacobians at all the equilibrium points
    J_1_ = feval(J_mat_alt, eq_1_(1), eq_1_(2));
    J_2_ = feval(J_mat_alt, eq_2_(1), eq_2_(2));
    J_3_ = feval(J_mat_alt, eq_3_(1), eq_3_(2));
    % Finding the Eigenvalues at all the equilibrium points
    Es_alt_params_1(:,i) = eig(J_1_);
    Es_alt_params_2(:,i) = eig(J_2_);
    Es_alt_params_3(:,i) = eig(J_3_);
    
    % Finding the largest valid frequency component (hopefully is the firing rate)
    [t_, X_] = ode15s(@(t, X)mle(t, X, I(i)), [0,300], [1.05*eq_2_(1), 1.05*eq_2_(2)]);  
    v_trace = X_(:,1);
    % Count number of times the membrane voltage rises close to the peak
    [peak_values, peak_times] = findpeaks(v_trace, 'MinPeakHeight', max(v_trace)-5);
    % Calculate firing rate
    num_spikes = length(peak_times); % Number of spikes
    firing_rates(i) = 1000 * num_spikes / length(v_trace); % Firing rate in Hz
end

% Plotting the eigenvalues wrt current
figure;
sgtitle('Equilibrium point(s) as I_{ext} varies')
subplot(211)
hold on
plot(I', real(Es_alt_params_1(1,:))', 'r')
plot(I', real(Es_alt_params_2(1,:))', 'g')
plot(I', real(Es_alt_params_3(1,:))', 'b')
plot(I', real(Es_alt_params_1(2,:))', 'r')
plot(I', real(Es_alt_params_2(2,:))', 'g')
plot(I', real(Es_alt_params_3(2,:))', 'b')
legend('First equilibrium point', 'Second equilibrium point', 'Third equilibrium point')
title('Re(eq)')
hold off
subplot(212)
hold on
plot(I, imag(Es_alt_params_1(1,:)), 'r')
plot(I, imag(Es_alt_params_2(1,:)), 'g')
plot(I, imag(Es_alt_params_3(1,:)), 'b')
plot(I, imag(Es_alt_params_1(2,:)), 'r')
plot(I, imag(Es_alt_params_2(2,:)), 'g')
plot(I, imag(Es_alt_params_3(2,:)), 'b')
legend('First equilibrium point', 'Second equilibrium point', 'Third equilibrium point')
title('Im(eq)')
hold off

% Plotting the firing rate wrt current
figure;
plot(I, firing_rates)
title('Modified model: firing rates modulation with input current')
xlabel('Input current (\mu A/cm^2)')
ylabel('Firing rate (Hz)')
%% Question 12
% Hodgin Huxley parameters
global G_K_HH E_K_HH G_Na_HH E_Na_HH G_L_HH E_L_HH C_HH Q_HH T_HH phi_HH V_r h_inf_r m_inf_r n_inf_r;
V_r = -60;
phi_HH = 1;
G_K_HH = 36;
G_Na_HH = 120; 
G_L_HH = 0.3;
E_K_HH = -72;
E_Na_HH = 55;
Q_HH = 3; 
T_HH = 6.3; 
C_HH = 1; 
I_ext = 0;
%% Question 13
% Finding E_L_HH
m_inf_r = m_inf(V_r);
n_inf_r = n_inf(V_r);
h_inf_r = h_inf(V_r);
E_L_HH = V_r - (I_ext-G_K_HH*n_inf_r^4*(V_r-E_K_HH)-G_Na_HH*m_inf_r^3*h_inf_r*(V_r-E_Na_HH))/(C_HH*G_L_HH);
% Showing response with I_ext = 10
I_ext = 10;
[t_, X_] = ode15s(@(t,X)hh(t, X, I_ext), [0,100], [V_r, n_inf_r, m_inf_r, h_inf_r]);
figure;
plot(t_, X_(:,1), LineWidth=2);
xlabel('Time (s)');
ylabel('Membrane voltage (mV)'); 
title('Hodgkin-Huxley time response for I_{ext} = 10\mu A/cm^2');
%% Question 14
I_ext = 0;
impulse = linspace(0,15,100);
max_V = zeros(1, 100);
for i = 1:100
    [~, X_] = ode15s(@(t,X)hh(t,X, I_ext), [0,100], [V_r+impulse(i)/C_HH, n_inf_r, m_inf_r, h_inf_r]);
    max_V(i) = max(X_(:,1));
end
thresh = (min(max_V(:)) + max(max_V(:)))/2;
thresh_i = 500;
for i = 1:100
    if max_V(i) >= thresh && i < thresh_i
        thresh_i = i;
    end
end
figure;
plot(impulse, max_V);
xlabel('q (A/F)');
ylabel('V (mV)');
title('Thresholding behavior in membrane voltage (V) wrt charge (q) injected');
%% Question 15
I_ext = 0;
eq_HH_0 = fsolve(@(X)hh(0, X, I_ext), [-60,0.05,0.05,0.05]);
for i = 8:12
    I_ext = i;
    % Obtaining equilibrium point
    eq_HH = fsolve(@(X)hh(0, X, I_ext), [eq_HH_0(1),eq_HH_0(2),eq_HH_0(3),eq_HH_0(4)]);
    % Linearizing the system about a small region
    syms V_ n_ m_ h_
    J_HH = jacobian([dVdt_HH(V_,n_,m_,h_,I_ext), dndt_HH(V_,n_), dmdt_HH(V_,m_), dhdt_HH(V_,h_)],[V_, n_, m_, h_]);
    J_HH_mat = matlabFunction(J_HH);
    J_HH_0 = feval(J_HH_mat, eq_HH(1), eq_HH(2), eq_HH(3), eq_HH(4));
    % Evaluating the 4-dimensional Eigenvalues of the linearized system
    E_HH(:,i-7) = eig(J_HH_0);
end
I_ext = 0;
%% Question 16
% Showing reduced-system response with I_ext = 10
I_ext = 10;
[t_, X_] = ode15s(@(t,X)hh_red(t, X, I_ext), [0,100], [V_r, n_inf_r]);
figure;
plot(t_, X_(:,1), LineWidth=2);
xlabel('Time (s)');
ylabel('Membrane voltage (mV)'); 
title('V-n reduced Hodgkin-Huxley time response for I_{ext} = 10\mu A/cm^2');
I_ext = 0;
% Showing charge spike response
for i = 1:100
    [~, X_] = ode15s(@(t,X)hh_red(t,X,I_ext), [0,100], [V_r+impulse(i)/C_HH, n_inf_r]);
    max_V(i) = max(X_(:,1));
end
thresh = (min(max_V(:)) + max(max_V(:)))/2;
thresh_i = 500;
for i = 1:100
    if max_V(i) >= thresh && i < thresh_i
        thresh_i = i;
    end
end
figure;
plot(impulse, max_V);
xlabel('q (A/F)');
ylabel('V (mV)');
title('Thresholding behavior in membrane voltage (V) wrt charge (q) injected (reduced system)');
%% Question 17
% Finding the resting state with the -3muA/cm^2 current clamp
I_ext = -3;
[t_1, X_clamped] = ode15s(@(t,X)hh(t, X, I_ext), [0,20], [V_r, n_inf_r, m_inf_r, h_inf_r]);
% Finding the response 
I_ext = 0;
[t_2, X_unclamped] = ode15s(@(t,X)hh(t, X, I_ext), [0,100], [X_clamped(end,1), X_clamped(end,2), X_clamped(end,3), X_clamped(end,4)]);
figure;
hold on
yyaxis left
plot([t_1; t_2+t_1(end)], [-3*ones(length(t_1),1); zeros(length(t_2),1)], LineWidth=2)
ylabel('Membrane Voltage (mV)')
yyaxis right
plot([t_1; t_2+t_1(end)], [X_clamped(:,1); X_unclamped(:,1)], LineWidth=2);
ylabel('I_{ext} (\mu A/cm^2)')
hold off
title('Voltage')
xlabel('Time (ms)')
%% Question 18
% Finding the resting state with the -3muA/cm^2 current clamp
I_ext = -3;
[t_1, X_clamped_2] = ode15s(@(t,X)hh_Vm(t, X, n_inf_r, h_inf_r, I_ext), [0,20], [V_r, m_inf_r]);
n_inf_clamped = n_inf(X_clamped_2(end,1));
h_inf_clamped = h_inf(X_clamped_2(end,1));
% Plotting the null-clines
figure;
fimplicit(@(V, m) dVdt_HH(V, n_inf_r, m, h_inf_r, I_ext), [-100 100 0 1]);
hold on
% Finding the equilibrium points of the initial system, and the
% corresponding Eigenvalues
eq_clamped_init_1 = fsolve(@(X)hh_Vm(0, X, n_inf_r, h_inf_r, I_ext), [-70,0]);
eq_clamped_init_2 = fsolve(@(X)hh_Vm(0, X, n_inf_r, h_inf_r, I_ext), [-50,0.2]);
eq_clamped_init_3 = fsolve(@(X)hh_Vm(0, X, n_inf_r, h_inf_r, I_ext), [50,0.8]);
J_c = jacobian([dVdt_HH(V_,n_inf_r,m_,h_inf_r,I_ext), dmdt_HH(V_,m_)],[V_, m_]);
J_c_mat = matlabFunction(J_c);
J_c_0_1 = feval(J_c_mat, eq_clamped_init_1(1), eq_clamped_init_1(2));
J_c_0_2 = feval(J_c_mat, eq_clamped_init_2(1), eq_clamped_init_2(2));
J_c_0_3 = feval(J_c_mat, eq_clamped_init_3(1), eq_clamped_init_3(2));
% Evaluating the Eigenvalues
E_c_1 = eig(J_c_0_1);
E_c_2 = eig(J_c_0_2);
E_c_3 = eig(J_c_0_3);
% Finding the response
I_ext = 0;
[t_2, X_unclamped_2] = ode15s(@(t,X)hh_Vm(t, X, n_inf_clamped, h_inf_clamped, I_ext), [0,100], [X_clamped_2(end,1), X_clamped_2(end,2)]);
% Finding the equilibrium point of the system right after the clamp is
% released
eq_unclamped_init = fsolve(@(X)hh_Vm(0, X, n_inf_clamped, h_inf_clamped, I_ext), [50,0.8]);
J_u = jacobian([dVdt_HH(V_,n_inf_clamped,m_,h_inf_clamped,I_ext), dmdt_HH(V_,m_)],[V_, m_]);
J_u_mat = matlabFunction(J_c);
J_u_0 = feval(J_c_mat, eq_unclamped_init(1), eq_unclamped_init(2));
% Evaluating the 4-dimensional Eigenvalues of the linearized system
E_u = eig(J_u_0);
% Plotting the null-clines
fimplicit(@(V, m) dVdt_HH(V, n_inf_clamped, m, h_inf_clamped, I_ext), [-100 100 0 1]);
hold on
fimplicit(@(V, m) dmdt_HH(V, m), [-100 100 0 1]);
legend('V nc, I = -3\mu A/cm^2, gates steady state at V_r', ...
    'V nc, I = 0\mu A/cm^2, gates steady state at V_{clamped}', ...
    'm nc')
hold off
%% Helper functions
% MLE V'
function dVdt = dVdt(V, w, I_ext)
    global G_Ca G_K G_L E_Ca E_K E_L V_1 V_2 C;
    % global I_ext; 
    dVdt = (I_ext-G_Ca*(0.5*(1+tanh((V-V_1)./V_2))).*(V-E_Ca)-G_K*w.*(V-E_K)-G_L*(V-E_L))/C;
end
% MLE w'
function dwdt = dwdt(V, w)
    global phi; 
    global V_3; 
    global V_4;
    dwdt = phi*(0.5*(1+tanh((V-V_3)./V_4)) - w).*cosh((V-V_3)/(2*V_4));
end
% MLE combined state equation
function dXdt = mle(t, X, I_ext)
    dXdt = [dVdt(X(1),X(2), I_ext); dwdt(X(1),X(2))];
end
% MLE alt params V'
function dVdt = dVdt_alt(V, w, I_ext)
    global G_Ca_alt G_K_alt G_L_alt E_Ca_alt E_K_alt E_L_alt V_1_alt V_2_alt C_alt;
    % global I_ext; 
    dVdt = (I_ext-G_Ca_alt*(0.5*(1+tanh((V-V_1_alt)./V_2_alt))).*(V-E_Ca_alt)-G_K_alt*w.*(V-E_K_alt)-G_L_alt*(V-E_L_alt))/C_alt;
end
% MLE alt params w'
function dwdt = dwdt_alt(V, w)
    global phi_alt; 
    global V_3_alt; 
    global V_4_alt;
    dwdt = phi_alt*(0.5*(1+tanh((V-V_3_alt)./V_4_alt)) - w).*cosh((V-V_3_alt)/(2*V_4_alt));
end
% MLE alt params combined state equation
function dXdt = mle_alt(t, X, I_ext)
    dXdt = [dVdt_alt(X(1),X(2), I_ext); dwdt_alt(X(1),X(2))];
end
% HH alpha_nn
function alpha_n = an(V)
    global phi_HH;
    if V ~= -50
        alpha_n = -0.01*phi_HH*(V+50)/(exp(-(V+50)/10)-1);
    else
        alpha_n = -0.01*phi_HH/(-1/10);
    end
end
% HH beta_n
function beta_n = Bn(V)
    global phi_HH;
    beta_n = 0.125*phi_HH*exp(-(V+60)/80);
end
% HH alpha_m
function alpha_m = am(V)
    global phi_HH;
    if V ~= -35
        alpha_m = -0.1*phi_HH*(V+35)/(exp(-(V+35)/10)-1);
    else
        alpha_m = -0.1*phi_HH/(-1/10);
    end
end
% HH beta_n
function beta_m = Bm(V)
    global phi_HH;
    beta_m = 4*phi_HH*exp(-(V+60)/18);
end
function beta_h = Bh(V)
    global phi_HH;
    beta_h = phi_HH/(exp(-(V+30)/10)+1);
end
% HH beta_n
function alpha_h = ah(V)
    global phi_HH;
    alpha_h = 0.07*phi_HH*exp(-(V+60)/20);
end
% HH V'
function dVdt = dVdt_HH(V, n, m, h, I_ext)
    global G_K_HH E_K_HH G_Na_HH E_Na_HH G_L_HH E_L_HH C_HH
    dVdt = (I_ext ...
        - G_K_HH*n^4.*(V-E_K_HH) ...
        - G_Na_HH*m^3.*h.*(V-E_Na_HH) ...
        - G_L_HH*(V-E_L_HH) ...
        )/C_HH;
end
% HH m'
function dmdt = dmdt_HH(V, m)
    dmdt = am(V)*(1-m) - Bm(V)*m;
end
% HH h'
function dhdt = dhdt_HH(V, h)
    dhdt = ah(V)*(1-h) - Bh(V)*h;
end
% HH n'
function dndt = dndt_HH(V, n)
    dndt = an(V)*(1-n) - Bn(V)*n;
end
% HH m_inf
function m_inf = m_inf(V)
    m_inf = am(V)/(am(V)+Bm(V));
end
function n_inf = n_inf(V)
    n_inf = an(V)/(an(V)+Bn(V));
end
function h_inf = h_inf(V)
    h_inf = ah(V)/(ah(V)+Bh(V));
end
% HH combined system ODE
function dXdt = hh(t, X, I_ext)
    dXdt = [dVdt_HH(X(1), X(2), X(3), X(4), I_ext);
        dndt_HH(X(1), X(2));
        dmdt_HH(X(1), X(3));
        dhdt_HH(X(1), X(4))
        ];
end
% HH combined reduced system ODE
function dXdt = hh_red(t, X, I_ext)
global h_inf_r
m_inf_v = m_inf(X(1));
    dXdt = [dVdt_HH(X(1), X(2), m_inf_v, h_inf_r, I_ext); 
        dndt_HH(X(1), X(2))];
end
% HH combined V-m reduced system ODE
function dXdt = hh_Vm(t, X, n_inf_v, h_inf_v, I_ext)
    dXdt = [dVdt_HH(X(1), n_inf_v, X(2), h_inf_v, I_ext); 
        dmdt_HH(X(1), X(2))];
end
