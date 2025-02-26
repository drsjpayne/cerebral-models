% function hypertensionexercise4
% SJP 20.ix.24
function hypertensionexercise4

% Exercise
[Dmax, Ton, Toff] = meshgrid(0:0.05:0.2, 30, 35:1:70);
cbf = zeros(size(Dmax));
for i = 1:size(Dmax,1)
    for j = 1:size(Dmax,2)
        for k = 1:size(Dmax,3)
            if Ton(i,j,k) < Toff(i,j,k)
                [q] = hypertensionexercise3(Dmax(i,j,k), Ton(i,j,k), Toff(i,j,k), 2);
                cbf(i,j,k) = q(end);
            end
        end
    end
end
yr = 5:40; scbf = squeeze(cbf);
figure(2)
plot(yr, scbf(5,:), 'b', yr, scbf(4,:), 'g', ...
    yr, scbf(3,:), 'y', yr, scbf(2,:), 'r', yr, scbf(1,:), 'k')
legend('0.2', '0.15', '0.1', '0.05', '0', 'location', 'northwest')
xlabel('Exercise duration (years)'), ylabel('CBF at age 80 years')

% Hypertension
[Dmax, Ton, Toff] = meshgrid(0:0.1:0.4, 30, 35:1:70);
cbf = zeros(size(Dmax));
for i = 1:size(Dmax,1)
    for j = 1:size(Dmax,2)
        for k = 1:size(Dmax,3)
            if Ton(i,j,k) < Toff(i,j,k)
                [q] = hypertensionexercise3(Dmax(i,j,k), Ton(i,j,k), Toff(i,j,k), 1);
                cbf(i,j,k) = q(end);
            end
        end
    end
end
yr = 5:40; scbf = squeeze(cbf);
figure(3)
plot(yr, scbf(5,:), 'b', yr, scbf(4,:), 'g', ...
    yr, scbf(3,:), 'y', yr, scbf(2,:), 'r', yr, scbf(1,:), 'k')
legend('40%', '30%', '20%', '10%', '0', 'location', 'southwest')
xlabel('Hypertension duration (years)'), ylabel('CBF at age 80 years')

% Comparison
cbf0 = hypertensionexercise3([0.4 0], [30 50], [100 70], 3);
cbf1 = hypertensionexercise3([0.4 0], [30 50], [50 70], 3);
cbf2 = hypertensionexercise3([0.4 0.2], [30 50], [100 100], 3);
cbf3 = hypertensionexercise3([0.4 0.2], [30 50], [50 100], 3);
t = 20:0.1:80;
figure(4)
plot(t, cbf0, 'k', t, cbf1, 'b', t, cbf2, 'r', t, cbf3, 'g')
xlabel('Age (years)'), ylabel('Blood flow')
legend('-AHT -Exercise', '+AHT -Exercise', '-AHT +Exercise', '+AHT +Exercise', 'location', 'southwest')
