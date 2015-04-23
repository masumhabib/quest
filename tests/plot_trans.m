clear all;
% clc;
% close all;

VGG = 0;
VDD = 0;
Vo = 0;
outPath = sprintf('./pn_TE/');
fileName = sprintf('TR_VGG%.3f_Vo%.3f_VDD%.3f.dat',...
                   VGG, Vo, VDD);
out = importTransResult([outPath, fileName]);

EE = out.E;
for iE = 1:length(EE)
    TE(iE) = real(trace(out.TE.M{iE}));
end

figure;
plot(EE, TE, 'k', 'LineWidth', 2);
xlabel('E (eV)');
ylabel('T (E)');
grid on;

