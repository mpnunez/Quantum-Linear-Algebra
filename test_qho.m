%test

clear; clc;

x = QHO();
x = x.solve();
x.plot_density();