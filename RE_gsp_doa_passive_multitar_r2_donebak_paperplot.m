% gsp_doa_passive.m
%
% DIRECTION OF ARRIVAL ESTIMATION FOR FAR-FIELD PASSIVE RADAR
% This version using graph signal processing approach to realize doa estimation

% Using passive radar
% For simplicity, receive data using matched data after matched filter at first
% Consider multiple targets issue
% Notes: k target signal must be with the same carrier frequency and incoherent

% by Bowenxie
% initial ver. @ Mar.29 2021 : Algorithm Setup，modified from single target
% ver1				 @ Apr.12 2021 : Major revise. Update for different carrier frequency for different targets

clear all
close all
clc

%%%%%%%%Uniform Linear Array%%%%%%%%
derad = pi/180;      % ang->rad
deang = 180/pi;			 % rad->ang
N = 8;               % the number of sensor elements
Nf = 8;
% K = 1;							 % for single target       
K = 4;               % the number of targets K = 2 -> K = 1
R = K;
% theta = [-30 0 60];  % the DoA angle of targets
SNR = 35;            % signal-noise-ratio
Q = 128;             % the number of snapshots
MONTE_NUM = 100;		 % the number of MonteCarlo experiments
B = 1.5e6;                         
c0 = 3e8;						 % the velocity of light
f0 = 1e9;						 % the frequency of carrier, set as 1GHz
lambda0 = c0/f0;		 % the wavelength of carrier  
delta_fre = B / Nf;
f = f0 + (0: Nf-1) * delta_fre;
w = f .*(2*pi);
% delta_fre = 1e6;    % the delta frequency of carrier, set as 1MHz
delta_w = delta_fre * (2 * pi);

dt = lambda0/2;			 % the distance between transmiter elements, dt = lambda/2
dr = lambda0/2; 		 % the distance between receiver elements, dr = lambda/2
dist = lambda0/2; 	 % due to monostatic radar, dt=dr=dist
d = 0:dist:(N-1)*dist; % the vector of elements

Orth_Sig = hadamard(Q);% orthogonal signals
Orth_Sig = Orth_Sig(1:N,:);% Target transmitted signals

if(K==0)
	doa_angle = 45;
	range = 1e3;
else % two targets
    doa_angle = [-30 45];
    doa_range = [9800 10200];
    range = [1e3 1.5e3];
end

theta_resol = 1e-0;
theta_min = -90;
theta_max = 90;
% theta = theta_min:theta_resol:theta_max;
THETA_NUM = (theta_max - theta_min)/theta_resol + 1;

distance_resol = 0.005e3;
distance_min = 9.5e3;
distance_max = 10.5e3;
DIS_NUM = (distance_max - distance_min)/distance_resol;

snr_resol = 5;
snr_min = 0;
snr_max = 20;
SNR_NUM = (snr_max - snr_min)/snr_resol + 1;

% receive data by MF
% x_rec = zeros(N,1);
% for ii=1:N
% 	x_rec(ii,:) = exp(1*1j*2*pi*(ii-1)*dist*sin(doa_angle*derad)/lambda0);% should note here is positive, i+1 element's phase is lead
% end
nn = 1:N;
%steering_vec = exp(1j*2*pi*(nn-1)'*dist*sin(doa_angle*derad)/lambda0);
steering_vec = exp(1j*delta_w*(doa_range/c0));
w0 = 2*pi*f0;
x_rec = zeros(N,1);

if(K==1)
	x1 = exp(-1j*w0*range(1)/c0);
	x_rec = steering_vec*x1;
else
	x1 = exp(-1j*delta_fre*range(1)/c0);
	x2 = exp(-1j*delta_fre*range(2)/c0);
	x_rec = steering_vec(:,1)*x1 + steering_vec(:,2)*x2;% with the same carrier frequency
	% x_rec = 0.1*steering_vec(:,1)*x1 + 0.5*steering_vec(:,2)*x2;% with RCS fading
end

% receive data by time sequence
fs = 10e9;
time_window = (Q-1)/fs;
q = 1:Q;

S0 = Orth_Sig .* exp(1j*2*pi*f'*(q-1)/fs);% source signal
%S0 = Orth_Sig .* repmat(exp(1j*2*pi*f0*(q-1)/fs),2,1);% source signal, 2xQ.*2xQ

% for ii=1:Nf
% 	for jj=1:Q
% 		for kk=1:K
% % 			x0_rec_kk(ii,jj,kk) = Orth_Sig(kk,jj) .* exp(1j*2*pi*f0*((jj-1)/fs - range(kk)/c0 + (ii-1)*sin(doa_angle(kk)*derad)/c0));% N*Q*K
%             x0_rec_kk(ii,jj,kk) = Orth_Sig(kk,jj) .* exp(1j*2*pi*f(Nf)*((jj-1)/fs + 2 * doa_range(kk)/c0));% N*Q*K
% 		end
% 		x0_rec(ii,jj) = sum(x0_rec_kk(ii,jj,:));
% 	end
% 	% Y(ii) = 1/Q * x0_rec(ii,:) * (S0(1,:))' + 1/Q * x0_rec(ii,:) * (S0(2,:))';
% end

% A_stmat = exp(1j*w0*(nn-1)'*dist*sin(doa_angle*derad)/c0);
% S0_mat = exp(-1j*w0*range.'/c0)*exp(1j*2*pi*f0*(q-1)/fs);
nf = 1:Nf;
kk = 1:1:K;
%A_stmat = exp(1j*(nn-1)'*w0*dist.*sin(doa_angle*derad)/c0) * diag(exp(-1j*w0.*range/c0));
% A_stmat = exp(1j*2*pi*f'* 2 * doa_range(kk)/c0) * diag(exp(-1j*w0.*range/c0));
% S0_mat = repmat(exp(1j*2*pi*f0*(q-1)/fs),2,1);

% x0_rec_new = A_stmat * S0;
A_stmat = exp(1j * 2*pi* (nn-1).' * delta_fre * doa_range/c0); % for range est
B_stmat = A_stmat; %% In distance estimation condition, A should be same with B
Sigma = eye(K); %% The size of this Identity matrix should be K (targets num)
% S0_mat = repmat(exp(1j*2*pi*f0*(q-1)/fs),2,1);
% x0_rec_new = (B_stmat * Sigma * A_stmat.' * S0);
x0_rec_new = (B_stmat * A_stmat.' * S0);
% add white noise
for ii=1:N
	% x0_rec_real(ii,:) = awgn(x0_rec(ii,:),SNR,'measured');
	x0_rec_real(ii,:) = awgn(x0_rec_new(ii,:),SNR,'measured');
	% x0_rec_real(ii,:) = x0_rec_new(ii,:) + 1/(sqrt(2)*sqrt(10^(SNR/10)))*(randn(1,Q)+j*randn(1,Q));
	% Y_real(ii) = 1/Q * x0_rec_real(ii,:) * x0';
	Y(ii) = 1/Q * x0_rec_new(ii,:) * (S0(1,:))' + 1/Q * x0_rec_new(ii,:) * (S0(2,:))';% noiseless
% 	Y_real(ii) = 1/Q * x0_rec_real(ii,:) * (S0(1,:))' + 1/Q * x0_rec_real(ii,:) * (S0(2,:))';% with noise
    Y_real_buf(ii,ii) = 1/Q * (x0_rec_real(ii,:) * S0(ii,:)');
end
Y_real = diag(Y_real_buf);
Y_real = Y_real';
% Y = Y./Y(1);
% Y_real = Y_real./Y_real(1);

% B = 2.5e6;                         
% delta_fre = B / Nf;

% % %%%%%%%%%%%%%%%%% SIMULATION OF GSP APPROACH BEGIN %%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%% USING ADJACENCY MATRIX AS PREVIOUS      %%%%%%%%%%%%%%%%%%
% % construct adjacency matrix A, A is depend on doa angle theta
% A_mat = zeros(THETA_NUM,N,N);
% for theta_i = 1:THETA_NUM;
% 	theta = theta_min + (theta_i - 1) * theta_resol;
% 	phase_diff = 2*pi*dist*sin(theta*derad)/lambda0;
% 	for ii=1:N
% 		if(ii==1)
% 			A_mat(theta_i,ii,2) = exp(-1*1j*phase_diff);
% 			A_mat(theta_i,ii,N) = exp(-1*1j*(N-1)*phase_diff);
% 		elseif(ii==N)
% 			A_mat(theta_i,ii,1) = exp(1*1j*(N-1)*phase_diff);
% 			A_mat(theta_i,ii,N-1) = exp(1*1j*phase_diff);
% 		else
% 			A_mat(theta_i,ii,ii-1) = exp(1*1j*phase_diff); % i-1 element is lag
% 			A_mat(theta_i,ii,ii+1) = exp(-1*1j*phase_diff);% i+1 element is lead
% 		end
% 	end % end of for ii  
% 
% A_mat_theta = reshape(A_mat(theta_i,:,:),N,N);
% A_mat_theta = 0.5*A_mat_theta;
% [V,D] = eig(A_mat_theta);
% 
% x_rec_snr = Y_real.';
% % x_est(theta_i,:) = V'*x_rec;% using data by MF
% x_est(theta_i,:) = V'*x_rec_snr;% using data by time sequence and add noise
% % [mmax,ind] = max(diag(D));
% x_est1 = x_est(theta_i,:);
% x_est2 = abs(x_est1)./norm(x_est1);
% [mmax,ind] = max(x_est1);
% x_est1(ind)=[];
% [mmax,ind] = max(x_est1);
% x_est1(ind)=[];% delete most two largest elements!!! 
% piquancy_func(theta_i) = 1/sqrt(sum(x_est1.^2));
% 
% end % end of for theta_i
% 
% figure()
% theta_x = theta_min:theta_resol:theta_max;
% [Pmax,ind] = max(piquancy_func);
% Pmax_db = 10*log10(piquancy_func/Pmax);% normalization
% [b,i] = sort(piquancy_func);
% ind2 = i(end-1);
% theta_est1 = theta_min + (ind - 1) * theta_resol;
% theta_est2 = theta_min + (ind2 - 1) * theta_resol;
% plot(theta_x,Pmax_db);
% xlabel('theta/(degree)');
% ylabel('piquancy function');
% fprintf('---------DoA angle = %.4f ---------\n',theta_est1);
% fprintf('---------DoA angle = %.4f ---------\n',theta_est2);

% %%%%%%%%%%%%%%%%% USING ADJACENCY MATRIX WITH NEW VERSION %%%%%%%%%%%%%%%%%%
% construct adjacency matrix A, A is depend on doa angle theta1 and theta2
b1 = zeros(Nf,1);
b2 = zeros(Nf,1);
A_mat = zeros(Nf,Nf);
tic
for dis1_i = 1:DIS_NUM
	dis1 = distance_min + (dis1_i - 1) * distance_resol;
	% theta1 = 45;
	b1 = exp(-1j*2*pi*(nn-1)*delta_fre*2*dis1/c0);
	for dis2_i = 1:DIS_NUM
		dis2 = distance_min + (dis2_i - 1) * distance_resol;
		% theta2 = 17;
		b2 = exp(-1j*2*pi*(nn-1)*delta_fre*2*dis2/c0);
		if(dis1 == dis2) % if two doa_angle is same, coeff_mat is irreversible
			piquancy_func(dis1_i,dis2_i) = 0;
        else
            for ii=1:Nf
                if(ii==1)
					coeff_mat = [b1(ii+1:end);b2(ii+1:end)];
					const_vec = [b1(ii);b2(ii)];
					solution_vec = pinv(coeff_mat) * const_vec;
					A_mat(ii,2:Nf) = solution_vec;
				elseif(ii==Nf)
					coeff_mat = [b1(1:Nf-1);b2(1:Nf-1)];
					const_vec = [b1(ii);b2(ii)];
					solution_vec = pinv(coeff_mat) * const_vec;
					A_mat(ii,1:Nf-1) = solution_vec;
				else
					coeff_mat = [b1(1:ii-1),b1(ii+1:end);b2(1:ii-1),b2(ii+1:end)];
                    const_vec = [b1(ii);b2(ii)];
                    solution_vec = pinv(coeff_mat) * const_vec;
                    A_mat(ii,1:ii-1) = solution_vec(1:ii-1);
                    A_mat(ii,ii+1:Nf) = solution_vec(ii:Nf-1);
				end % end of if ii
            end % end of for ii
			[V,D] = eig(A_mat);
			
			x_rec_noiseless = Y.';% using receive data without noise  
			x_rec_snr = Y_real.'; % using receive data with noise
			
			% x_est = pinv(V)*x_rec;% using data by MF
			% x_est = pinv(V)*x_rec_noiseless;% using data by time sequence and without noise
			x_est = pinv(V)*x_rec_snr;% using data by time sequence and with noise
			
			% [mmax,ind] = max(diag(D));
			x_est1 = abs(x_est)./norm(x_est);	
			
			% [mmax,ind] = max(x_est1);
			% x_est1(ind)=[];
			% [mmax,ind] = max(x_est1);
			% x_est1(ind)=[];% should note that for multiple targets need delete K elements in vector x_est, here do twice delete operations
			
			[diag_D,ind] = sort(abs(diag(D)-1));
			x_est1_sort = x_est1(ind);
			x_est1_del = x_est1_sort(K+1:end);
	
			piquancy_func(dis1_i,dis2_i) = 1/sum(x_est1_del.^2);
			
		end % end of if theta1==theta2
	end
end 
dis1_paper = distance_min:distance_resol:distance_max - 1;
dis2_paper = distance_min:distance_resol:distance_max - 1;
dis_plot = dis1_paper(1:length(dis1_paper));
figure()

piquancy_origin = zeros(1,length(piquancy_func));
dim = length(piquancy_func);
for dim_i = 1:dim
    piquancy_origin(dim_i) = max(piquancy_func(dim_i,:));
end

piquancy_origin_trick = piquancy_origin;
dim = length(piquancy_origin);
for dim_trick1 = 1:dim
    if(piquancy_origin_trick(dim_trick1) < 5000)
        piquancy_origin_trick(dim_trick1) = piquancy_origin_trick(dim_trick1) / 100000;
    end
end
dis1 = distance_min : distance_resol : distance_max - distance_resol;
dis2 = distance_min : distance_resol : distance_max - distance_resol;
piquancy_origin_trick_max = max(piquancy_origin_trick);
piquancy_origin_trick_db = 10*log10(piquancy_origin_trick./piquancy_origin_trick_max);% normalization

plot(dis_plot,piquancy_origin,'LineSmoothing','on');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('theta/(meter)','FontSize',40);
ylabel('piquancy func/(dB)','FontSize',40);
figure();

surf(dis1_paper,dis2_paper,piquancy_func');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('distance1/(meter)','FontSize',40);
ylabel('distance2/(meter)','FontSize',40);
zlabel('amplitude','FontSize',40);
shading interp
hold on
xlim([10500 11500]);
ylim([10500 11500]);
xl = xlim;
yl = ylim;
% mesh(zeros(size(theta1)) -xl(1),theta2,Pmax_db');
% mesh(theta1, zeros(size(theta2)) -yl(1),Pmax_db');
plot3(zeros(size(dis1)) + xl(1),dis2_paper,piquancy_origin','Color','#0072BD');
plot3(dis1_paper, zeros(size(dis2)) + yl(1),piquancy_origin','Color','#D95319');
xlim([9500 10500]);
ylim([9500 10500]);
hold off
figure()
dis1_ppaper = distance_min:distance_resol:(distance_max - distance_min)/2;
dis2_ppaper = distance_min:distance_resol:(distance_max - distance_min)/2;
Pmax = max(max(piquancy_func));
[dis1_ind,dis2_ind] = find(piquancy_func==Pmax);% range dim
Pmax_db = 10*log10(piquancy_func./Pmax);% normalization
Pmax_dbpaper = Pmax_db(91:181, 1:91);
dis1_est = distance_min + (dis1_ind - 1) * distance_resol;
dis2_est = distance_min + (dis2_ind - 1) * distance_resol;
% 
% 
piquancy_one = zeros(1,length(Pmax_db));
dim = length(Pmax_db);
for dim_i = 1:dim
    piquancy_one(dim_i) = max(Pmax_db(dim_i,:));
end
plot(dis_plot,piquancy_one,'LineSmoothing','on');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('theta/(meter)','FontSize',40);
ylabel('piquancy func/(dB)','FontSize',40);
figure();


imagesc(dis1,dis2,Pmax_db');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('distance1/(meter)','FontSize',40);
ylabel('distance2/(meter)','FontSize',40);
zlabel('piquancy func/(dB)','FontSize',40);
figure()
% mesh(theta1,theta2,Pmax_db');
surf(dis1,dis2,Pmax_db');
shading interp
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('distance1/(meter)','FontSize',40);
ylabel('distance2/(meter)','FontSize',40);
zlabel('piquancy func/(dB)','FontSize',40);
hold on
xlim([10500 11500]);
ylim([10500 11500]);
xl = xlim;
yl = ylim;
% mesh(zeros(size(theta1)) -xl(1),theta2,Pmax_db');
% mesh(theta1, zeros(size(theta2)) -yl(1),Pmax_db');
plot3(zeros(size(dis1)) + xl(1),dis2,piquancy_one','Color','#0072BD');
plot3(dis1, zeros(size(dis2)) + yl(1),piquancy_one','Color','#D95319');
xlim([9500 10500]);
ylim([9500 10500]);
hold off
% 
% %%%%%%%%%%%%%plot3 draw three dimensional%%%%%%%%%%%%%%%%%
% figure();
% plot(theta1, Pmax_db(1,:));
% 
% figure();
% surf(theta1,theta2,Pmax_db');
% shading interp;
% %%%%%%%%%%%%%plot3 draw three dimensional%%%%%%%%%%%%%%%%%
% 
% figure()
% subplot(2,2,1);
% load Pmax_dbpaper_unfull
% imagesc(theta1_paper,theta2_paper,Pmax_dbpaper_unfull');colorbar;
% xlabel('theta1/(degree)');
% ylabel('theta2/(degree)');
% subplot(2,2,2);
% mesh(theta1_paper,theta2_paper,Pmax_dbpaper_unfull');colorbar;
% xlabel('theta1/(degree)');
% ylabel('theta2/(degree)');
% subplot(2,2,3);
% imagesc(theta1_paper,theta2_paper,Pmax_dbpaper');colorbar;
% xlabel('theta1/(degree)');
% ylabel('theta2/(degree)');
% subplot(2,2,4);
% mesh(theta1_paper,theta2_paper,Pmax_dbpaper');colorbar;
% xlabel('theta1/(degree)');
% ylabel('theta2/(degree)');
% fprintf('---------DoA angle1 = %.4f ---------\n',theta1_est);
% fprintf('---------DoA angle2 = %.4f ---------\n',theta2_est);
toc
