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
%% Artifact bug: difference in null quench at t=dt. when dt->0, the bug is gone.
%(0.6621147649485886,0)
%(0.6396046674142885,-0.008351946232544946)
clear
%close all
 idata = 8;
% 
Omega1=0.5;omega=50;besselj(0,Omega1/omega)*1.2
% 0.5532056086826694
Omega1=0.5;omega=10;besselj(0,Omega1/omega)*1.2
% 0.5530341614556746
Omega1=0.5;omega=5;besselj(0,Omega1/omega)*1.2
% 0.5525006458956896
Omega1=0.5;omega=2;besselj(0,Omega1/omega)*1.2
% 0.5488615656462038
Omega1=0.5;omega=1;besselj(0,Omega1/omega)*1.2
% 0.5371926921601495
Omega1=0.5;omega=0.5;besselj(0,Omega1/omega)*1.2
% 0.505885932881813
Omega1=0.5;omega=0.1;besselj(0,Omega1/omega)*1.2
% 0.07909928736239571
Omega1 = 0.5;
omega_grid = [50 10 5 2 1 0.5 0.1];
DeltaGRST = [
    0.5532056086826694
    0.5530341614556746
    0.5525006458956896
    0.5488615656462038
    0.5371926921601495
    0.505885932881813
    0.07909928736239571
    ];
filename = {
    'hi_0.5Omega1_0.5omega_50.dat',...
    'hi_0.5Omega1_0.5omega_10.dat',...
    'hi_0.5Omega1_0.5omega_5.dat',...
    'hi_0.5Omega1_0.5omega_2.dat',...
    'hi_0.5Omega1_0.5omega_1.dat',...
    'hi_0.5Omega1_0.5omega_0.5.dat',...
    'hi_0.5Omega1_0.5omega_0.1.dat',...
    'hi_1Omega1_0.6omega_0.1.dat',...
};
%for idata = 1:length(filename)
for idata = 1:7
data = load(filename{idata});
figure(idata)
t = data(:,1);
Delta = data(:,2) + 1i* data(:,3);
ht = data(:,4);
N0 = data(:,5);
N1 = data(:,6);
plot(t,abs(Delta),'r',t,DeltaGRST(idata)*ones(1,length(t)),'b',...
    t, N0, 'm',t,N1,'k')
xlabel('t/(E_F^{-1})')
ylabel('|\Delta(t)|/E_F')
title(['h_{eff}=0.5,\alpha_0=1.2,E_b=0.2,\Omega_1 = ',num2str(Omega1), ...
    ', \omega = ' , num2str(omega_grid(idata)),...
    ', \Delta_{grst}=',num2str(DeltaGRST(idata))])
set(gca,'fontsize',16)
legend('|\Delta(t)|','\Delta_{grst}','n_s/n','n_t/n')
saveas(figure(idata),['fig',num2str(idata),'.eps'],'epsc')
end