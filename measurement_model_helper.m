%% This is a helper script to find the measurement model

clear all;
close all;

% load('dataset3.mat')
syms x y z fu cu fv cv b

ul = fu * x / z + cu;
vl = fv * y / z + cv;
ur = fu * (x - b) / z + cu;
vr = fv * y / z + cv;

g = [ul; vl; ur; vr];
X = [x;y;z];

G = vpa(jacobian(g, X), 5)

% pretty(G)