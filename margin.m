clearvars
close all
clc

s = tf('s');
F = 1/(s^4+80*s^3+11200*s^2+200000*s);

margin(F);