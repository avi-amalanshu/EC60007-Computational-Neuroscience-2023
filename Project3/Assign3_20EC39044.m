%% Computational Neuroscience -- Project 3
%   -- Avi Amalanshu
%   -- 20EC39044
%% Setup
clc;
clear;
load('data_cn_project_iii_a23.mat');
%% Question 1
[R_xx, lags] = autocorr(Stimulus, 50); % this will return 0-50 autocorr

plot([-lags(end:-1:1) lags],[R_xx(end:-1:1) R_xx]); % symmetric, -50 to -1 concat w/ 0 to 50
xlabel('Lag (ms)');
ylabel('Autocorrelation');
title('Autocorrelation of Stimulus');
%% Question 2
psth = zeros(4, 20000);

for j = 1:4
    for rep = 1:50
        psth(j,:) = psth(j,:) + histcounts(All_Spike_Times{j,rep}*1000, 0:20000)*1000/50;
    end
end

figure;

for i = 1:4
    subplot(2,2,i)
    plot(psth(i,:))
    xlabel('time (ms)');
    ylabel('spike rate');
    title(['PSTH (neuron ' num2str(i) ')']);
end

%% Question 3
bins = [10, 20, 50, 100, 200, 500];
avgs = zeros(6,4,50);
vars = zeros(6,4,50);

for i = 1:6
    bin = bins(i);
    for j = 1:4
        for rep = 1:50
            spike_count = histcounts(All_Spike_Times{j,rep}*1000, 0:bin:20000);
            avgs(i,j,rep) = mean(spike_count);
            vars(i,j,rep) = var(spike_count);
        end
    end
end
count = 1;
figure;
for i = 1:6
    for j = 1:4
        subplot(6, 4, count)
        avg_values = reshape(avgs(i, j, :), 1, []);
        var_values = reshape(vars(i, j, :), 1, []);

        sz = 25;
        scatter(avg_values, var_values, sz, 'filled')
        hold on
        plot([min(avg_values), max(avg_values)],[min(avg_values), max(avg_values)], 'red')
        xlabel('mean')
        ylabel('variance')
        count = count + 1;

        title(sprintf('bin width = %d, neuron = %d', bins(i), j))
    end
end
%% Question 4

sta = zeros(4,100);

for i = 1:4
    spikes = 0;
    for j = 1:50
        wind = All_Spike_Times{i,j} <= 15;
        spikes = spikes + nnz(wind);
        spike_train = All_Spike_Times{i,j}(wind);
        [m,n] = size(spike_train);
        for k = 1:n
            spike_time = max(1,ceil(1000*(spike_train(k))));
            stim = Stimulus(max(spike_time-99, 1):spike_time);
            sta(i,:) = sta(i,:) + [zeros(1,100-length(stim)), stim]; 
        end
    end
    sta(i,:) = sta(i,:)/spikes;
end

sta(:,:) = sta/max(R_xx); 

figure;
sgtitle('Spike Triggered Average without correction');

for i = 1:4
    subplot(2, 2, i);
    plot(0:99, sta(i, :));
    ylim([-0.8, 0.8]);
    xlabel('\tau');
    ylabel('R_{px}(\tau)');
    title(['Neuron ' num2str(i)]);
end
sgtitle('STA'); 

[R_xx_new, lags] = autocorr(Stimulus, 99);

R_xx_new_mat = toeplitz(R_xx_new);


h = ((R_xx_new_mat\sta')');

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(0:99, h(i, :));
    xlabel('\tau');
    ylabel('R_{px}(\tau)');
    title(['Neuron ' num2str(i)]);
end
sgtitle('STA (suppressed C^{-1}_{xx})')
%% Question 5
y = zeros(4, 15000);
pred = conv(Stimulus(1:15000), h(1, :)); y(1, :) = pred(1:15000);
pred = conv(Stimulus(1:15000), h(2, :)); y(2, :) = pred(1:15000);
pred = conv(Stimulus(1:15000), h(3, :)); y(3, :) = pred(1:15000);
pred = conv(Stimulus(1:15000), h(4, :)); y(4, :) = pred(1:15000);
    
psth_15 = psth(:, 1:15000);   
bin_size = 30;
y_ = zeros(4,500);
lam = zeros(4,500);
for i = 1:ceil(15000/bin_size)
    bin_end = i*bin_size;
    if bin_end>15000
        bin_end = 15000;
   end
        y_(1,i) = mean(y(1,(1+(i-1)*bin_size):bin_end));
        y_(2,i) = mean(y(2,(1+(i-1)*bin_size):bin_end));
        y_(3,i) = mean(y(3,(1+(i-1)*bin_size):bin_end));
        y_(4,i) = mean(y(4,(1+(i-1)*bin_size):bin_end));
    
        lam(1,i) = mean(psth_15(1,(1+(i-1)*bin_size):bin_end));
        lam(2,i) = mean(psth_15(2,(1+(i-1)*bin_size):bin_end));
        lam(3,i) = mean(psth_15(3,(1+(i-1)*bin_size):bin_end));
        lam(4,i) = mean(psth_15(4,(1+(i-1)*bin_size):bin_end));
end
ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [20, 20, 0.1];

[xdata, ydata] = prepareCurveData(y_(1,:), lam(1,:));
mdl1 = fit(xdata, ydata, ft, opts);
[xdata, ydata] = prepareCurveData(y_(2,:), lam(2,:));
mdl2 = fit(xdata, ydata, ft, opts);
[xdata, ydata] = prepareCurveData(y_(3,:), lam(3,:));
mdl3 = fit(xdata, ydata, ft, opts);
[xdata, ydata] = prepareCurveData(y_(4,:), lam(4,:));
mdl4 = fit(xdata, ydata, ft, opts);

figure;

subplot(2,2,1)
xlim([-0.1,0.1])
plot(mdl1)
grid on
hold on
scatter(y_(1,:),lam(1,:), 'filled')
hold off

subplot(2,2,2)
xlim([-1,1])
plot(mdl2)
grid on
hold on
scatter(y_(2,:),lam(2,:), 'filled')
hold off

subplot(2,2,3)
xlim([-1,1])
plot(mdl3)
grid on
hold on
scatter(y_(3,:),lam(3,:), 'filled')
hold off

subplot(2,2,4)
xlim([-0.1,0.1])
plot(mdl4)
grid on
hold on
scatter(y_(4,:),lam(4,:), 'filled')

pred1 = conv(Stimulus(15001:20000), h(1, :)); 
pred1 = pred1(1:5000);

pred2 = conv(Stimulus(15001:20000), h(2, :)); 
pred2 = pred2(1:5000);

pred3 = conv(Stimulus(15001:20000), h(3, :)); 
pred3 = pred3(1:5000);

pred4 = conv(Stimulus(15001:20000), h(4, :)); 
pred4 = pred4(1:5000);

pred1 = mdl1(pred1) + 0.0005*randn(5000, 1);
pred2 = mdl2(pred2) + 0.0005*randn(5000, 1);
pred3 = mdl3(pred3) + 0.0005*randn(5000, 1);
pred4 = mdl4(pred4) + 0.0005*randn(5000, 1);

ground_truth1 = psth(1, 15001:20000);
ground_truth2 = psth(2, 15001:20000);
ground_truth3 = psth(3, 15001:20000);
ground_truth4 = psth(4, 15001:20000);

%% Question 6
R1 = corrcoef(psth(1, 15001:20000),pred1);
R2_1 = R1(2)^2

R2 = corrcoef(psth(2, 15001:20000),pred2);
R2_2 = R2(2)^2

R3 = corrcoef(psth(3, 15001:20000),pred3);
R2_3 = R3(2)^2


R4 = corrcoef(psth(4, 15001:20000),pred4);
R2_4 = R4(2)^2

figure;
subplot(2,2,1)
scatter(ground_truth1,pred1')
title('PSTH vs pred for neuron 1')


subplot(2,2,2)
scatter(ground_truth2,pred2')
title('PSTH vs pred for neuron 2')


subplot(2,2,3)
scatter(ground_truth3,pred3')
title('PSTH vs pred for neuron 3')


subplot(2,2,4)
scatter(ground_truth4,pred4')
title('PSTH vs pred for neuron 4')

A = zeros(100,1);
B = zeros(100,1);

initial_r_sq = R2_2;
old_r_sq = R2_2;
new_r_sq = old_r_sq;
pos = zeros(1,1);
count = 0;
h_new = h(:,:);

while (old_r_sq - new_r_sq) < 1e-6 && count < 1e6
    old_r_sq=new_r_sq;
    min_=Inf;
    for i = 1:100 
        if abs(h_new(2,i))<min_ && abs(h_new(2,i))~=0
            min_=abs(h(2,i));
            pos=i;
        end
    end
    h_new(2,pos)=0;
    pos(1,1)=0;
    pred2 = conv(Stimulus(15001:20000), h_new(2, :)); pred2 = pred2(1:5000);
    pred2 = mdl2(pred2) + 0.0005*randn(5000, 1);
    new_r=corrcoef(psth(2, 15001:20000), pred2); new_r_sq=new_r(2)^2;
    count=count+1;
    A(count,1)=count;
    B(count,1)=new_r_sq;
end


figure
scatter(A,B)
xlabel('iterations')
ylabel('r^2'); 
title('Pruning parameters for neuron 2')

max_correlation_coeff2 = max(B);

f_t2 = fft(h_new(2,:));
f_t1 = fft(h(1,:));

vect = -50:49;
vect = vect';

figure
sgtitle('Linear filters h(t)');

for i = 1:4
    subplot(2, 2, i);
    stem(h(i, :));
    title(['Neuron ' num2str(i)]);
end

figure
sgtitle('Linear filters h(t) in frequency domain');

for i = 1:4
    subplot(2, 2, i);
    if i == 2
        plot(vect, f_t2);
    else
        plot(vect, fft(h(i, :)));
    end
    title(['Neuron ' num2str(i)]);
end
%% Part 7
q = [0 0.001 0.01 0.1 1 10 100];
numIterations = 100;
numSegments = 8;
segmentLength = 100;
MI = zeros(4, numIterations, length(q));

for iter = 1:numIterations
    t = randperm(19901, numSegments);
    t = t / 1000; % Convert to seconds
    spike_seg = cell(numSegments, 50);
    
    for rep = 1:50
        for n = 1:4
            spikes = All_Spike_Times{n, rep};
            for i = 1:numSegments
                t_start = t(i);
                t_end = t(i) + segmentLength / 1000;
                spike_seg{i, rep} = spikes(spikes >= t_start & spikes < t_end);
            end
        end
    end
    
    for m = 1:length(q)
        confusion_matrix = zeros(numSegments, numSegments);
        for i = 1:numSegments
            for rep = 1:50
                mean_dist = meandist(i, rep, spike_seg, q(m));
                mean_dist(i) = mean_dist(i) * 50 / 49;
                [~, k] = min(mean_dist);
                confusion_matrix(i, k) = confusion_matrix(i, k) + 1;
            end
        end
        confusion_matrix = confusion_matrix / 50;
        MI(:, iter, m) = MI(:, iter, m) + this_MI(confusion_matrix) / 100;
    end
end

conf_90 = confidence_interval(0.9, MI);

a = MI(:, 1, :);
b = conf_90(:, 1, :);

for n = 1:4
    subplot(2,2,n);
    plot(log10(q), a(n,:));
    hold on
    plot(log10(q), a(n,:)-b(n,:), 'k--');
    plot(log10(q), a(n,:)+b(n,:), 'k--');
    [~,p] = max(a(n,:));
    p = plot(log10(q(p)), a(n,1,p), '*');
    set(p, 'linewidth', 2)
    hold off
    title(['Discrimination - Neuron ' num2str(n)])
    xlabel('log_1_0(q)')
    ylabel('MI(q)')
end
sgtitle('Mutual information vs lag weight (90% conf. interv.)')
%% Helper Functions
function spikeDistance = spike_dist(time_list_I, time_list_J, q)
    num_spikes_I = length(time_list_I);
    num_spikes_J = length(time_list_J);

    if q == 0
       spikeDistance = abs(num_spikes_I - num_spikes_J);
       return
    elseif q == Inf
       spikeDistance = num_spikes_I + num_spikes_J;
       return
    end
    
    score_matrix = zeros(num_spikes_I + 1, num_spikes_J + 1);
    score_matrix(:, 1) = (0:num_spikes_I)';
    score_matrix(1, :) = (0:num_spikes_J);
    
    if(num_spikes_I > 0 && num_spikes_J > 0)
       for i = 2:num_spikes_I + 1
          for j = 2:num_spikes_J + 1
             score_matrix(i, j) = min([score_matrix(i-1, j) + 1, score_matrix(i, j-1) + 1, score_matrix(i-1, j-1) + q * abs(time_list_I(i-1) - time_list_J(j-1))]);
          end
       end
    end
    spikeDistance = score_matrix(num_spikes_I + 1, num_spikes_J + 1);
end

function mean_dist = meandist(i,rep,spike_seg,q)
mean_dist = zeros(1,8);
for i1 = 1:8
    for rep1 = 1:50
        if (rep1 == rep && i1 == i)
            continue
        end
        mean_dist(i1) = mean_dist(i1) + spike_dist(spike_seg{i,rep},spike_seg{i1,rep1},q);
    end
end
mean_dist = mean_dist/50;
end

function conf = confidence_interval(interval, dist)
    scale = std(abs(dist),0,2);
    T_multiplier = tinv(interval-1/2, 99);
    conf = T_multiplier*scale/sqrt(99);
end


function MI = this_MI(dist)
    MI = 0;
    for i=1:size(dist,1)
        for j=1:size(dist,2)
            if(dist(i,j)~=0)
                MI=MI+dist(i,j)/size(dist,1)*log2(dist(i,j)/sum(dist(:,j)));
            end
        end
    end
end
