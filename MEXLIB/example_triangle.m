clc
clear 
close all

load('mesh_data_new');
disp(p)
%[pout, tout] = tricall(p');
simpplot(p,t)
