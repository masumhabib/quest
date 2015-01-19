clear all;

VGG = 0;
VDD = 0;
Vo = 0;
outPath = sprintf('./pn_TE/');
fileName = sprintf('dbg_pot.dat');
out = importTransResult([outPath, fileName]);

disp(['File: ', [outPath, fileName]]);
[X,Y,Z,V] = importPotential([outPath, fileName]);

figure;
plot(X, Y, '.');
daspect([1,1,1]);