%Amogh THOUTU%

clear all
close all
%ELECTRON MODELLING

m0 = 9.109*10^-31; %rest mass of electron (kg)
meff = 0.26* m0; %effective mass
w = 200*10^-9; l = 100*10^-9; % dimensions (m);

 %add barriers to 2_D to plot
 %upper box
 u1 = [0.8*10^-7, 0.6*10^-7];
 u2 = [1.2*10^-7, 1*10^-7];
 u_h = linspace(u1(1),u2(1),100);% horizontal line
 u_v = linspace(u1(2),u2(2),100);%vertical line
 %lower box
 L1 = [0.8*10^-7, 0];
 L2 = [1.2*10^-7, 0.4*10^-7];
 L_h = linspace(L1(1),L2(1),100);% horizontal line
 L_v = linspace(L1(2),L2(2),100);%vertical line
 figure(2)
 hold on
 %upper box plot
 plot(u_h,u1(2)*ones(1,length(u_h)),'k','LineWidth',2.0);
 plot(u_h, u2(2)*ones(1,length(u_h)),'k','LineWidth',2.0);
 plot(u1(1)*ones(1,length(u_v)),u_v, 'k','LineWidth',2.0);
 plot(u2(1)*ones(1,length(u_v)),u_v, 'k','LineWidth',2.0);
 %lower box plot
 plot(L_h,L1(2)*ones(1,length(L_h)),'k','LineWidth',2.0);
 plot(L_h, L2(2)*ones(1,length(L_h)),'k','LineWidth',2.0);
 plot(L1(1)*ones(1,length(L_v)),L_v, 'k','LineWidth',2.0);
 plot(L2(1)*ones(1,length(L_v)),L_v, 'k','LineWidth',2.0);

%assign random location for particles
npart = 1000; %number of particles
partloc = [w*rand(npart,1),l*rand(npart,1)]; % particle locations
m = 116.768*10^3; %mean
s = 1000;% std

%remove particles from boxes
%upper box
 u1_log = (partloc(:,1)> u1(1)) &( partloc(:,1)<u2(1)); % upper box verical lines
 u2_log = (partloc(:,2)> u1(2)) & (partloc(:,2) <u2(2));% upper box horizontal lines
 u_remove  = u1_log & u2_log; 
 partloc(u_remove,:) = nan;

 %lower box
 L1_log = (partloc(:,1)> L1(1)) &( partloc(:,1)<L2(1)); % lower box verical lines
 L2_log = (partloc(:,2)> L1(2)) & (partloc(:,2) <L2(2));% lower box horizontal lines
 L_remove  = L1_log & L2_log; 
 partloc(L_remove,:) = nan;
% plot(partloc(:,1),partloc(:,2), "*"); %electron density plot
 

% assign random velocity
 vth = s.*randn(npart,1)+ m; 
 %plot histogram for vth
 figure(1)
 histogram(vth,'Normalization','probability');
 xlabel('vth (m/s)')
 ylabel('probability')
 title(['vth normal distribution (mean: ', num2str(m),' m/s) (std: ',num2str(s),' )'])

%assign random direction
 rand_ang = 360* rand(npart,1); % random angles in degrees
 part_vel = [vth.*cos(rand_ang),vth.*sin(rand_ang)]; % particle[vx,vy] m/s
    
 %update particle location
 n_timesteps = 1000;
 npart_plot = 10;% number of particles to plot
 part_plot = zeros(npart,1); 
 part_plot(round(linspace(1,npart,npart_plot)),1)= 1;

 dt = (l./100)/(116768);% value of dt in seconds

 
 %Probability of scattering
 mt  = 0.2*10^-12; %mean time between collisions
 Pscat = 1 - exp(-dt/mt);
 
 temp  = zeros(1, n_timesteps);

 for i =1:n_timesteps
    
     partloc(:,1) = partloc(:,1) + part_vel(:,1)*dt; %x-pos
     partloc(:,2) = partloc(:,2) + part_vel(:,2)*dt; %y-pos
     
     %Electron scattering
    
     c1_rand = rand(size(part_vel,1),1);% col 1
     c2_rand = rand(size(part_vel,1),1);% col 2
      
     c1_log = Pscat > c1_rand;%logical array
     c2_log = Pscat > c2_rand;
     
     part_vel(c1_log, 1) = (s.*randn()+ m).*cos(360* rand());%assign new velocities and directions
     part_vel(c2_log,2)  = (s.*randn()+ m).*sin(360* rand());
     
     %reflect particles off boxes
     % top box
     ru1_log = (partloc(:,1)> u1(1)) & (partloc(:,1) < u2(1)) & (partloc(:,2) > u1(2)-0.01*u1(2)) & (partloc(:,2) < u1(2)+0.01*u1(2));
     ru2_log = (partloc(:,2)> u1(2)) & (partloc(:,2) < u2(2)) & (partloc(:,1) > u1(1)-0.01*u1(1)) & (partloc(:,1) < u1(1)+0.01*u1(1)); 
     ru3_log = (partloc(:,2)> u1(2)) & (partloc(:,2) < u2(2)) & (partloc(:,1) > u2(1)-0.01*u1(1)) & (partloc(:,1) < u2(1)+0.01*u1(1));
     
     part_vel(ru1_log,2) = -part_vel(ru1_log,2);
     part_vel(ru2_log,1) = -part_vel(ru2_log,1);
     part_vel(ru3_log,1) = -part_vel(ru3_log,1);
     
     %bottom box
     rL1_log = (partloc(:,1)> L1(1)) & (partloc(:,1) < L2(1)) & (partloc(:,2) > L2(2)-0.015*L2(2)) & (partloc(:,2) < L2(2)+0.015*L2(2));
     rL2_log = (partloc(:,2)> L1(2)) & (partloc(:,2) < L2(2)) & (partloc(:,1) > L1(1)-0.01*L1(1)) & (partloc(:,1) < L1(1)+0.01*L1(1));
     rL3_log = (partloc(:,2)> L1(2)) & (partloc(:,2) < L2(2)) & (partloc(:,1) > L2(1)-0.01*u1(1)) & (partloc(:,1) < L2(1)+0.01*u1(1));
     
      part_vel(rL1_log,2) = -part_vel(rL1_log,2);
      part_vel(rL2_log,1) = -part_vel(rL2_log,1);
      part_vel(rL3_log,1) = -part_vel(rL3_log,1);
      

     %Boundary conditions
     
     b1_log = partloc(:,1)<0;% logical arrays for boundary conditions
     b2_log = partloc(:,1)>w;
     b3_log = partloc(:,2)<0 | partloc(:,2) >l;
    
     
     partloc(b1_log,1) = w;  % set new particle positions
     partloc(b2_log,1) = 0;
     part_vel(b3_log,2) = -part_vel(b3_log,2);
 
     
     %Temperature
     part_vel_square = part_vel.^2;
     avg_vel = sum(sqrt(sum(part_vel_square,2)))/npart;% average velocity 
     temp(i) = (0.5*m0*(avg_vel).^2)/(1.380*10^-23);% temperature (KE/k) or (1/2 mv^2 / k)
     avg_temp = sum(temp)/i;
    
     
     %plotting
     figure(2);
     
     %subplot(2,1,1)
%      plot(part_plot.*partloc(:,1),part_plot.*partloc(:,2),'b.','LineWidth',2.0)
%      hold on
%      xlim([0 w]);
%      ylim([0 l]);
%      title('Electron trajectories')
%      xlabel('width (m)')
%      ylabel('length (m)')
     
     %subplot(2,1,2)
     plot(i,temp(i),'g.')
     hold on
     xlim([0 n_timesteps])
     title(['Temperature (K): ',num2str(temp(i)),'  Average: ', num2str(avg_temp),' K'])
     xlabel('number of time steps')
     ylabel('semiconductor temperature (K)')
      pause(0.01)
     
 end    
 hold off