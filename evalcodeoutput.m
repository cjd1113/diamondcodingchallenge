

fval_c = importdata('fval.txt');
alph_c = importdata('falph.txt');

iters = 1:1:length(fval_c);

figure(1)
plot(iters,fval_c)

figure(2)
plot(iters,alph_c)