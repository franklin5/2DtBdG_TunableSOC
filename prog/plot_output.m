%%%%% plot the data.
%% ground state solution
% clear
% clc
% load superfluid.dat
% h = superfluid(:,1);
% Delta = superfluid(:,2);
% Mu = superfluid(:,3);
% Eg = superfluid(:,4);
% set(gca,'fontsize',16);
% figure(1)
% plot(h, Delta, 'r', h, Mu, '--', h, Eg, 'k','linewidth',2)
% xlabel('h/E_F')
% legend('\Delta','\mu','E_g')
%%
clear
%close all
 idata = 1;
% Omega1=0.5;omega=50;besselj(0,Omega1/omega)*1.2
filename = {
    'hi_0.5Omega1_0.5omega_50.dat',...
    'hi_0.5Omega1_0.5omega_10.dat',...
    'hi_0.5Omega1_0.5omega_5.dat',...
    'hi_0.5Omega1_0.5omega_2.dat',...
    'hi_0.1omega_0.2.dat',...
    'hi_0.9omega_10.5.dat',...
    'hi_0.9omega_8.dat',...
    'hi_0.9omega_6.dat',...
    'hi_0.9omega_4.dat',...
    'hi_0.9omega_2.dat',...
    'hi_0.9omega_1.dat',...
    'hi_0.9omega_0.75.dat',...
    'hi_0.9omega_0.5.dat',...
    'hi_0.9omega_0.25.dat',...
     'hi_0.9omega_0.1.dat',...
     'hi_0.9omega_0.05.dat'};
%for idata = 1:length(filename)
data = load(filename{idata});
figure(idata)
t = data(:,1);
Delta = data(:,2) + 1i* data(:,3);
ht = data(:,4);
plot(t,abs(Delta),'r')%,t,ht,'b')
%end