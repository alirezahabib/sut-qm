%{
    Alireza Habibzadeh
    Student No. 99109393
    Quantum Mechanics 1
    Spring 2023 - Prof. Vaezi
%}


close all
clear
clc
tic

dx = 0.02;
x = (-10:dx:10)';
m = 1; % mass
hbar = 1; % Planck constant

n_x = length(x);
D1_x = (diag(ones(1, n_x - 1), 1) - diag(ones(1, n_x), 0)) / dx;
D1_x(n_x, 1) = 1;
P1_x = (hbar/1i) * D1_x; % momentum operator

% potential operator 
V_x = diag(abs(x));
%V_x = diag(-exp(-abs(x)/4) .* cos(x));
%V_x = diag(x.^4/4 - x.^2/2);
%V_x = diag(exp(-x.^2/2));
%V_x = diag(9999999999 * (abs(x) > 2));


H = (P1_x' * P1_x + P1_x * P1_x')/2/(2*m) + V_x; % Hamiltonian operator
[psi, Hd] = eig(full(H)); % diagonalizing Hamiltonian operator
E = diag(Hd); % energy eigenvalues

figure(1);
hold on;

psi_1 = psi(:, 1); % ground-state energy wavefunction
plot(x, abs(psi_1).^2, 'Linewidth', 1) % plotting probability distribution function (pdf)

psi_2 = psi(:, 2); % 1st excited state wave function
plot(x, abs(psi_2).^2 , 'Linewidth', 1) % plotting probability distribution function (pdf)

psi_3 = psi(:, 3); % 2nd excited state wave function
plot(x, abs(psi_3).^2 , 'Linewidth', 1) % plotting probability distribution function (pdf)

xlabel('x');
ylabel('|\psi(x)|^2');
legend('ground state', '1st excited-state', '2nd excited-state');
title('probability distribution function');

subtitle('$V(x) = |x|$','interpreter', 'latex');
%subtitle('$V(x) = -e^{-|x|/4}\cos{x}$','interpreter', 'latex');
%subtitle('$V(x) = x^4/4 - x^2/2$','interpreter', 'latex');
%subtitle('$V(x) = e^{-x^2/2}$','interpreter', 'latex');
%subtitle('$V(x) = \infty \mathrm{if} \, |x| > 2,\, 0 \, \mathrm{o.w.}$','interpreter', 'latex');

% plotting energy eigenvalues (theory: E(n) = (n+1/2) \hbar \omega),
% n=0, 1, 2, ...
figure(2);
nE = 10;
a = zeros(nE, 1);
b = E(1:nE);
plot(a, b, '-s','MarkerSize', 10, ...
    'MarkerEdgeColor', 'red' , ...
    'MarkerFaceColor', [1 .6 .6])
text(a + 0.1, b, num2str(b), 'Fontsize', 8);
ylabel('E(n)');
title('Lowest 10 energy eigenvalues');
toc

