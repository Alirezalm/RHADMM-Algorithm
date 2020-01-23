%% tank system

clc, clear, close all

%% time structure

t0 = 0;         tf = 20;        dt = 0.1;
time = t0: dt: tf - dt; Nt = numel(time);


%% System Structure

A = [0.875 0.1250; 0.1250 0.8047];   B = [0.3 3]';


x = [4 3]';

k  = 0;
for t = t0: dt: tf - dt
    k = k + 1;
   x = A * x;
   
   X(:,k) = x;
end

plot(time,X)