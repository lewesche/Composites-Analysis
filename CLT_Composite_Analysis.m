%Leif Wesche
%Classical Laminate Theory
%Composite Laminate Stress, Strain, Deformation, Failure Analysis

close all
clear all
clc

% Enter Material Properties
E1=181*10^3; %MPa
E2=10.3*10^3; %MPa
G12=7.17*10^3; %MPa
v12=0.28;
h=0.2/1000; %m

% Enter Play Orientation, starting at the bottom
theta=[0, 45, -45, 90, 0]; %deg

% Enter Ply Strengths
s1p=1500; %MPa
s1m=1500; %MPa
s2p=40; %MPa
s2m=246; %MPa
s12=68; %MPa

% Enter Thermal and Moisture Properties 
alpha1=1000;
alpha2=1000;
alpha12=1000;
beta1=1000;
beta2=1000;
beta12=1000;


% Enter Loading conditions
Nx=0.4; %MPa*m
Ny=-0.3; %MPa*m
Nxy=0; %MPa*m
Mx=0; %
My=0; %
Mxy=0; %

%Temperature and Moisture properties differential
dTt= 50;
dTb= 20;
dCt=0;
dCb=0;

dT0=(dTt+dTb)/2;
dT1=dTt-dT0;
dT2=dTb-dT0;
dC0=(dCt+dCb)/2;
dC1=dCt-dC0;
dC2=dCb-dC0;
alpha=[alpha1; alpha2; alpha12];
beta=[beta1; beta2; beta12];

% Compile Linear Z vector
z=linspace(-h*length(theta)/2, h*length(theta)/2, length(theta)+1);

D=(1-(E2./E1)*(v12.^2));

Q=[E1/D, E2*v12/D, 0;
   E2*v12/D, E2/D, 0;
   0, 0, G12];

%Calculate Q_bar matricies
for i=1:length(theta)
    
T_sig{i}=[cosd(theta(i)).^2, sind(theta(i)).^2, 2.*cosd(theta(i)).*sind(theta(i));
      sind(theta(i)).^2, cosd(theta(i)).^2, -2.*cosd(theta(i)).*sind(theta(i));
      -1.*cosd(theta(i)).*sind(theta(i)), 1.*cosd(theta(i)).*sind(theta(i)), (cosd(theta(i)).^2 - sind(theta(i)).^2)];
  
T_eps{i}=[cosd(theta(i)).^2, sind(theta(i)).^2, 1.*cosd(theta(i)).*sind(theta(i));
       sind(theta(i)).^2, cosd(theta(i)).^2, -1.*cosd(theta(i)).*sind(theta(i));
       -2.*cosd(theta(i)).*sind(theta(i)), 2.*cosd(theta(i)).*sind(theta(i)), (cosd(theta(i)).^2 - sind(theta(i)).^2)];

Q_bar{i}=inv(T_sig{i})*Q*T_eps{i};
end

%Calculate Thermal Load 
Nht=[0;0;0];
Mht=[0;0;0];
for i=1:length(theta)
%Uniform Part
Nht=Nht+dT0*(z(i+1)-z(i)).*Q_bar{i}*(T_eps{i}*alpha);
Nht=Nht+dT0*(z(i+1)-z(i)).*Q_bar{i}*(T_eps{i}*beta);
Mht=Mht+dT0*(1/2)*((z(i+1)).^2-(z(i)).^2).*Q_bar{i}*(T_eps{i}*alpha);
Mht=Mht+dT0*(1/2)*((z(i+1)).^2-(z(i)).^2).*Q_bar{i}*(T_eps{i}*beta);
%Linear Part
Nht=Nht+dT1*(1/2)*((z(i+1)).^2-(z(i)).^2).*Q_bar{i}*(T_eps{i}*alpha);
Nht=Nht+dT1*(1/2)*((z(i+1)).^2-(z(i)).^2).*Q_bar{i}*(T_eps{i}*beta);
Mht=Mht+dT1*(1/3)*((z(i+1)).^3-(z(i)).^3).*Q_bar{i}*(T_eps{i}*alpha);
Mht=Mht+dT1*(1/3)*((z(i+1)).^3-(z(i)).^3).*Q_bar{i}*(T_eps{i}*beta);
end





%Calculate and sum A, B, and D matricies
A=[0,0,0;0,0,0;0,0,0];
B=A;
D=A;
for i=1:length(theta)
A= A + Q_bar{i}*(z(i+1)-z(i));
B= B + (1/2)*Q_bar{i}*((z(i+1)).^2-(z(i)).^2);
D= D + (1/3)*Q_bar{i}*((z(i+1)).^3-(z(i)).^3);
end

%Compile ABD
ABD1=cat(2, A, B);
ABD2=cat(2, B, D);
ABD=cat(1, ABD1, ABD2);

%Calculate midpoint strains and curvature
Loading=[Nx; Ny; Nxy; Mx; My; Mxy];
Loading_ht=[Nht; Mht];
Prod_ht=inv(ABD)*Loading_ht;
Prod=inv(ABD)*Loading;
Strain_mid=Prod(1:3);
Curvature=Prod(4:6);
Strain_mid_plot=Prod(1:3)+Prod_ht(1:3);
Curvature_plot=Prod(4:6)+Prod_ht(4:6);


%Calculate Strain at Ply Top and Bottom
for i=1:length(z)
Strain_x(i)=Strain_mid(1)+z(i)*Curvature(1);
Strain_y(i)=Strain_mid(2)+z(i)*Curvature(2);
Strain_xy(i)=Strain_mid(3)+z(i)*Curvature(3);

Strain_x_plot(i)=Strain_mid_plot(1)+z(i)*Curvature_plot(1);
Strain_y_plot(i)=Strain_mid_plot(2)+z(i)*Curvature_plot(2);
Strain_xy_plot(i)=Strain_mid_plot(3)+z(i)*Curvature_plot(3);
end

%Plot Strains Throughout Composite
figure
plot(Strain_x_plot, z, Strain_y_plot, z, Strain_xy_plot, z, 'linewidth', 2)
title('Strains Throughout The Composite'); ylabel('Position (m)'); xlabel('Strain')
legend('\epsilon_{x}', '\epsilon_{y}', '\gamma_{xy}', 'Location', 'Best')

%Calculate Strain at Center-Ply
z_ply=linspace((-h*length(theta)+h)/2, (h*length(theta)-h)/2, length(theta));
for i=1:length(z_ply)
Strain_x(i)=Strain_mid(1)+z(i)*Curvature(1);
Strain_y(i)=Strain_mid(2)+z(i)*Curvature(2);
Strain_xy(i)=Strain_mid(3)+z(i)*Curvature(3);
end

%Calculate Stress at Center-Ply
for i=1:length(z_ply) 
Stress_refxy{i}=Q_bar{i}*[Strain_x(i); Strain_y(i); Strain_xy(i)];    
Stress_dummy_refxy=Stress_refxy{i};
Stress_x_dummy(i)=Stress_dummy_refxy(1);
Stress_y_dummy(i)=Stress_dummy_refxy(2);
Stress_xy_dummy(i)=Stress_dummy_refxy(3);


Stress_ref12{i}=T_sig{i}*Stress_refxy{i};
Stress_dummy_ref12=Stress_ref12{i};
Stress_1_dummy(i)=Stress_dummy_ref12(1);
Stress_2_dummy(i)=Stress_dummy_ref12(2);
Stress_12_dummy(i)=Stress_dummy_ref12(3);
end

%Expand Stress and z vectors to show ply stress discontinuity
Stress_x=[];
Stress_y=[];
Stress_xy=[];
for i=1:length(z_ply) 
Stress_x=[Stress_x, Stress_x_dummy(i), Stress_x_dummy(i)];
Stress_y=[Stress_y, Stress_y_dummy(i), Stress_y_dummy(i)];
Stress_xy=[Stress_xy, Stress_xy_dummy(i), Stress_xy_dummy(i)];
end
z_stress=[z(1)];
for i=1:(length(z_ply)-1) 
z_stress=[z_stress, z(i+1), z(i+1)];
end
z_stress=[z_stress, z(end)];

%Plot Stress Throughout Composite
figure
plot(Stress_x, z_stress, Stress_y, z_stress, Stress_xy, z_stress, 'linewidth', 2)
title('Stress At The Center of Each Ply'); ylabel('Position (m)'); xlabel('Stress (MPa)')
legend('\sigma_{x}', '\sigma_{y}', '\tau_{xy}', 'Location', 'Best')

%Evaluate Failure

%Stress Failure Theory
 for i=1:length(z_ply) 
     %Direction 1
     if Stress_1_dummy(i) > s1p         
        X = ['Maximum Stress Theory Failure in Direction 1 at ply ',num2str(i)];
        disp(X)
     end
     if Stress_1_dummy(i) < -1*s1m
        X = ['Maximum Stress Theory Failure in Direction 1 at ply ',num2str(i)];
        disp(X)         
     end
     %Direction 2
     if Stress_2_dummy(i) > s2p         
        X = ['Maximum Stress Theory Failure in Direction 2 at ply ',num2str(i)];
        disp(X)
     end
     if Stress_2_dummy(i) < -1*s2m
        X = ['Maximum Stress Theory Failure in Direction 2 at ply ',num2str(i)];
        disp(X)         
     end     
 end

 
%Strain Failure Theory
e1p=s1p/E1;
e1m=s1m/E1;
e2p=s2p/E2;
e2m=s2m/E2;
sig1_3 = @(sig2) E1/(v12*E2)*(sig2+s2m);
sig1_4 = @(sig2) E1/(v12*E2)*(sig2-s2p);
sig2_13 = (s2m/(v12*v12) - E2*s1p/(E1*(v12)))/(1-1/(v12*v12));
sig2_14 = -1*(s2p/(v12*v12) + E2*s1p/(E1*(v12)))/(1-1/(v12*v12));
sig2_24 = (E2*s1m/(v12*E1) - s2p/(v12*v12))/(1-1/(v12*v12));
sig2_23 = (E2*s1m/(v12*E1) + s2m/(v12*v12))/(1-1/(v12*v12));
strain_envelope_1 = [sig1_3(sig2_13), sig1_4(sig2_14) sig1_4(sig2_24), sig1_3(sig2_23), sig1_3(sig2_13)];
strain_envelope_2 = [sig2_13, sig2_14, sig2_24, sig2_23, sig2_13];

for i=1:length(z_ply)
     %Direction 1
     if (Stress_1_dummy(i)./E1 - v12.*Stress_2_dummy(i)./E2) > e1p         
        X = ['Maximum Strain Theory Failure in Direction 1 at ply ',num2str(i)];
        disp(X)
     end
     if (Stress_1_dummy(i)./E1 - v12.*Stress_2_dummy(i)./E2) < -1*e1m
        X = ['Maximum Strain Theory Failure in Direction 1 at ply ',num2str(i)];
        disp(X)         
     end
     %Direction 1
     if (Stress_2_dummy(i)./E2 - v12.*Stress_1_dummy(i)./E1) > e2p         
        X = ['Maximum Strain Theory Failure in Direction 2 at ply ',num2str(i)];
        disp(X)
     end
     if (Stress_2_dummy(i)./E2 - v12.*Stress_1_dummy(i)./E1) < -1*e2m
        X = ['Maximum Strain Theory Failure in Direction 2 at ply ',num2str(i)];
        disp(X)         
     end     
end




%Tsai-Wu Theory
F1 = (1/s1p) - (1/s1m);
F11 = 1/(s1p*s1m);
F2 = (1/s2p) - (1/s2m);
F22 = 1/(s2p*s2m);
F12 = -0.5 * sqrt(F11*F22);
F66 = 1/(s12*s12);
for i=1:length(z_ply);
     %Tsai-Wu
     if (Stress_1_dummy(i)*F1 + Stress_2_dummy(i)*F2 + F11*(Stress_1_dummy(i)).^2 + F22*(Stress_2_dummy(i)).^2 + F66*(Stress_12_dummy(i)).^2 + 2*F12*(Stress_1_dummy(i))*(Stress_2_dummy(i)) ) > 1        
        X = ['Tsai-Wu Theory Failure at ply ',num2str(i)];
        disp(X)
     end 
end


%Plot Stress Failure Envelopes (Neglect s12)
figure
plot([-s1m, -s1m, s1p, s1p, -s1m], [-s2m, s2p, s2p, -s2m, -s2m], strain_envelope_1, strain_envelope_2, 'Linewidth', 3); grid on
title('Stress 1 and Stress 2 Faliure Analysis'); ylabel('s2 (Pa)'); xlabel('s1 (Pa)')
hold on
for i=1:length(z_ply)
fimplicit(@(s1,s2) F1*s1 + F2*s2 + F11*(s1)^2 + F22*(s2)^2 + F66*(Stress_12_dummy(i))^2 + 2*F12*s1*s2 -1, [-s1m*2.5 s1p*2.5 -s2m*2.5 s2p*2.5], 'linewidth', 3)
end
hold on
plot(Stress_1_dummy, Stress_2_dummy, 'go', 'Linewidth', 3)
legend('Max Stress Envelop', 'Max Strain Envelope', '\gamma_{xy}', 'Location', 'Best')










