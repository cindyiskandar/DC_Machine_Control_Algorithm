clear all; close all; clc;
sim('tc_simpower.slx');
figure()
subplot(5,1,1)
plot(t,DutyCycleIGBT14,'linewidth',2),grid;
xlabel('Temps (s)');
ylabel('Tension (V)');
title('Tension aux bornes des IGBT 1 et 4 en fonction du temps');
subplot(5,1,2)
plot(t,DutyCycleIGBT23,'linewidth',2),grid;
xlabel('Temps (s)');
ylabel('Tension (V)');
title('Tension aux bornes des IGBT 2 et 3 en fonction du temps');
subplot(3,1,1)
plot(t,ArmatureVoltage,'linewidth',2),grid;
xlabel('Temps (s)');
ylabel('Tension (V)');
title('Tension aux bornes de l''armature de la MCC en fonction du temps');
subplot(3,1,2)
plot(t,ArmatureCurrent,'linewidth',2),grid;
xlabel('Temps (s)');
ylabel('Courant (A)');
title('Courant aux bornes de l''armature de la MCC en fonction du temps');
subplot(3,1,3)
plot(t,MotorSpeed,'linewidth',2),grid;
xlabel('Temps (s)');
ylabel('Vitesse (rad/s)');
title('Vitesse de la MCC en fonction du temps');
