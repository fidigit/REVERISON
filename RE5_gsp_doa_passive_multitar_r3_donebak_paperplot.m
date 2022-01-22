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
K = 3;               % the number of targets
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
    doa_range = [9800 10000 10250];
    range = [1.5e3 10e3 15e3];
end

theta_resol = 5e-0;
theta_min = -90;
theta_max = 90;
% theta = theta_min:theta_resol:theta_max;
THETA_NUM = (theta_max - theta_min)/theta_resol + 1;

distance_resol = 0.010e3;
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

% S0 = Orth_Sig .* repmat(exp(1j*2*pi*f0*(q-1)/fs),3,1);% source signal, 2xQ.*2xQ
S0 = Orth_Sig .* exp(1j*2*pi*f'*(q-1)/fs);

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
A_stmat = exp(1j * 2*pi* (nn-1).' * delta_fre * doa_range/c0); % for range est
B_stmat = A_stmat; %% In distance estimation condition, A should be same with B
% Sigma = eye(K); %% The size of this Identity matrix should be K (targets num)
Sigma = [1, 0, 0; 0, 10, 0; 0, 0, 2];
% S0_mat = repmat(exp(1j*2*pi*f0*(q-1)/fs),2,1);
x0_rec_new = (B_stmat * Sigma * A_stmat.' * S0);
% x0_rec_new = (B_stmat * A_stmat.' * S0);
% add white noise
for ii=1:N
	% x0_rec_real(ii,:) = awgn(x0_rec(ii,:),SNR,'measured');
	x0_rec_real(ii,:) = awgn(x0_rec_new(ii,:),SNR,'measured');
	% x0_rec_real(ii,:) = x0_rec_new(ii,:) + 1/(sqrt(2)*sqrt(10^(SNR/10)))*(randn(1,Q)+j*randn(1,Q));
	% Y_real(ii) = 1/Q * x0_rec_real(ii,:) * x0';
	Y(ii) = 1/Q * x0_rec_new(ii,:) * (S0(1,:))' + 1/Q * x0_rec_new(ii,:) * (S0(2,:))' + 1/Q * x0_rec_new(ii,:) * (S0(3,:))';% noiseless
	Y_real_buf(ii,ii) = 1/Q * (x0_rec_real(ii,:) * S0(ii,:)');
end
Y_real = diag(Y_real_buf);
Y_real = Y_real';
tic
% %%%%%%%%%%%%%%%%% USING ADJACENCY MATRIX WITH NEW VERSION %%%%%%%%%%%%%%%%%%
% construct adjacency matrix A, A is depend on doa angle theta1 and theta2
% b1 = zeros(Nf,1);
% b2 = zeros(Nf,1);
% A_mat = zeros(Nf,Nf);
% A_mat_save_one = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_Nf = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_two = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_three = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_four = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_five = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_six = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% A_mat_save_seven = zeros(DIS_NUM,DIS_NUM,DIS_NUM,1,Nf);
% % A_mat_save = zeros(DIS_NUM,DIS_NUM,DIS_NUM,Nf,Nf);
load A_mat3_source.mat;
% 
% parfor dis1_i = 1:DIS_NUM
% 	dis1 = distance_min + (dis1_i - 1) * distance_resol;
% 	% theta1 = 45;
% 	b1 = exp(-1j*2*pi*(nn-1)*delta_fre*2*dis1/c0);
% 	for dis2_i = 1:DIS_NUM
% 		dis2 = distance_min + (dis2_i - 1) * distance_resol;
% 		% theta2 = 17;
% 		b2 = exp(-1j*2*pi*(nn-1)*delta_fre*2*dis2/c0);
%         for dis3_i = 1:DIS_NUM
%             dis3 = distance_min + (dis3_i - 1) * distance_resol;
%             % theta2 = 17;
%             b3 = exp(-1j*2*pi*(nn-1)*delta_fre*2*dis3/c0);
%             coeff_mat = [b1(2:end);b2(2:end);b3(2:end)];
%             const_vec = [b1(1);b2(1);b3(1)];
%             solution_vec = pinv(coeff_mat) * const_vec;
%             for jj = 2 : 8
%                 A_mat_save_one(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec(jj - 1);
%             end
%             
%             coeff_mat_NF = [b1(1:Nf-1);b2(1:Nf-1);b3(1:Nf-1)];
%             const_vec_Nf = [b1(ii);b2(ii);b3(ii)];
%             solution_vec_Nf = pinv(coeff_mat_NF) * const_vec_Nf;
%             for jj = 1 : 7
%                 A_mat_save_Nf(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_Nf(jj);
%             end
%             
%             coeff_mat_two = [b1(1:2-1),b1(2+1:end);b2(1:2-1),b2(2+1:end);b3(1:2-1),b3(2+1:end)];
%             const_vec_two = [b1(2);b2(2);b3(2)];
%             solution_vec_two = pinv(coeff_mat_two) * const_vec_two;
%             for jj = 1 : 1
%                 A_mat_save_two(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_two(jj);
%             end
%             for jj = 3 : 8
%                 A_mat_save_two(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_two(jj-1);
%             end
%             
%             coeff_mat_three = [b1(1:3-1),b1(3+1:end);b2(1:3-1),b2(3+1:end);b3(1:3-1),b3(3+1:end)];
%             const_vec_three = [b1(3);b2(3);b3(3)];
%             solution_vec_three = pinv(coeff_mat_three) * const_vec_three;
%             for jj = 1 : 2
%                 A_mat_save_three(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_three(jj);
%             end
%             for jj = 4 : 8
%                 A_mat_save_three(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_three(jj-1);
%             end
%             
%             coeff_mat_four = [b1(1:4-1),b1(4+1:end);b2(1:4-1),b2(4+1:end);b3(1:4-1),b3(4+1:end)];
%             const_vec_four = [b1(4);b2(4);b3(4)];
%             solution_vec_four = pinv(coeff_mat_four) * const_vec_four;
%             for jj = 1 : 3
%                 A_mat_save_four(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_four(jj);
%             end
%             for jj = 5 : 8
%                 A_mat_save_four(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_four(jj-1);
%             end
%             
%             coeff_mat_five = [b1(1:5-1),b1(5+1:end);b2(1:5-1),b2(5+1:end);b3(1:5-1),b3(5+1:end)];
%             const_vec_five = [b1(5);b2(5);b3(5)];
%             solution_vec_five = pinv(coeff_mat_five) * const_vec_five;
%             for jj = 1 : 4
%                 A_mat_save_five(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_five(jj);
%             end
%             for jj = 6 : 8
%                 A_mat_save_five(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_five(jj-1);
%             end
%             
%             coeff_mat_six = [b1(1:6-1),b1(6+1:end);b2(1:6-1),b2(6+1:end);b3(1:6-1),b3(6+1:end)];
%             const_vec_six = [b1(6);b2(6);b3(6)];
%             solution_vec_six = pinv(coeff_mat_six) * const_vec_six;
%             for jj = 1 : 5
%                 A_mat_save_six(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_six(jj);
%             end
%             for jj = 7 : 8
%                 A_mat_save_six(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_six(jj-1);
%             end
%             
%             coeff_mat_seven = [b1(1:7-1),b1(7+1:end);b2(1:7-1),b2(7+1:end);b3(1:7-1),b3(7+1:end)];
%             const_vec_seven = [b1(7);b2(7);b3(7)];
%             solution_vec_seven = pinv(coeff_mat_seven) * const_vec_seven;
%             for jj = 1 : 6
%                 A_mat_save_seven(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_seven(jj);
%             end
%             for jj = 8 : 8
%                 A_mat_save_seven(dis1_i,dis2_i,dis3_i,1,jj) = solution_vec_seven(jj-1);
%             end
% 
%         end % end of for ii
%     end % end of if theta1==theta2
% end 
% toc
% A_mat_save(:,:,:,1,:) = A_mat_save_one;
% A_mat_save(:,:,:,8,:) = A_mat_save_Nf;
% A_mat_save(:,:,:,2,:) = A_mat_save_two;
% A_mat_save(:,:,:,3,:) = A_mat_save_three;
% A_mat_save(:,:,:,4,:) = A_mat_save_four;
% A_mat_save(:,:,:,5,:) = A_mat_save_five;
% A_mat_save(:,:,:,6,:) = A_mat_save_six;
% A_mat_save(:,:,:,7,:) = A_mat_save_seven;
% toc
% save A_mat3_source A_mat_save;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GPU_Amat = gpuArray(A_mat_save);
tic
parfor dis1_i = 1:DIS_NUM
	for dis2_i = 1:DIS_NUM
        for dis3_i = 1:DIS_NUM
            [V,D] = eig(squeeze(GPU_Amat(dis1_i,dis2_i,dis3_i,:,:))); 
            x_rec_snr = Y_real.'; % using receive data with noise
            x_est = pinv(V)*x_rec_snr;% using data by time sequence and with noise
            x_est1 = abs(x_est)./norm(x_est);	
            [diag_D,ind] = sort(abs(diag(D)-1));
            x_est1_sort = x_est1(ind);
            x_est1_del = x_est1_sort(K+1:end);
            piquancy_func(dis1_i,dis2_i,dis3_i) = 1/sum(x_est1_del.^2);
        end % end of if theta1==theta2
	end
end 
toc
% save piquancy_func.mat
% load piquancy_func_r3_200plot.mat
[Pmax,index]=max(piquancy_func(:));
[dis1_ind_b, dis2_ind_b, dis3_ind_b] = ind2sub([200,200,200], index);
dis1_ind = (dis1_ind_b -1) - 50;
dis2_ind = (dis2_ind_b -1) - 50;
dis3_ind = (dis3_ind_b -1) - 50;
dis1_full = distance_min : distance_resol : distance_max - distance_resol;
dis2_full = distance_min : distance_resol : distance_max - distance_resol;

Pmax_db = 10*log10(piquancy_func/Pmax);% normalization
piquancy_full = zeros(length(piquancy_func),length(piquancy_func));
piquancy_full_first = zeros(length(piquancy_func),length(piquancy_func));
piquancy_full_second = zeros(length(piquancy_func),length(piquancy_func));
piquancy_full_third = zeros(length(piquancy_func),length(piquancy_func));
dim = length(piquancy_func);
% for dim_trick1 = 1:dim
%     for dim_trick2 = 1:dim
%         for dim_trick3 = 1:dim
%             if(piquancy_func(dim_trick1,dim_trick2,dim_trick3) < 15000)
%                 piquancy_func(dim_trick1,dim_trick2,dim_trick3) = piquancy_func(dim_trick1,dim_trick2,dim_trick3) / 10000;
%             end
%         end
%     end
% end
for dim_i = 1:dim
    for dim_j = 1:dim
        piquancy_full_first(dim_i,dim_j) = max(piquancy_func(:,dim_i,dim_j));
        piquancy_full_second(dim_i,dim_j) = max(piquancy_func(dim_j,:,dim_i));
        piquancy_full_third(dim_i,dim_j) = max(piquancy_func(dim_i,dim_j,:));
    end
end
piquancy_full = (piquancy_full_first + piquancy_full_second + piquancy_full_third) / 3;
surf(dis1_full,dis1_full,piquancy_full);
shading interp
hold on
xlim([10500 10501]);
ylim([10500 10501]);
xl = xlim;
yl = ylim;
% mesh(zeros(size(theta1)) -xl(1),theta2,Pmax_db');
% mesh(theta1, zeros(size(theta2)) -yl(1),Pmax_db');
xlim([9500 10500]);
ylim([9500 10500]);
plot3(zeros(size(dis1_full)) + xl(1),dis2_full,piquancy_full','Color','#0072BD');
plot3(dis1_full, zeros(size(dis2_full)) + yl(1),piquancy_full','Color','#D95319');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('dis1/(degree)','FontSize',40);
ylabel('dis2/(degree)','FontSize',40);
zlabel('piquancy func/(dB)','FontSize',40);

% load piquancy_func_r3_200plot.mat;
% piquancy_full_one = zeros(length(piquancy_full),length(piquancy_full));
piquancy_full_one = zeros(length(piquancy_full),1);
dim = length(piquancy_full);
for dim_i = 1:dim
     if(max(piquancy_full(:,dim_i)) < 10000000)
         piquancy_full_one(dim_i) = max(piquancy_full(:,dim_i)) / 100000;
     else
        piquancy_full_one(dim_i) = max(piquancy_full(:,dim_i));
     end
end
% for dim_i = 1:dim
% %     if(max(piquancy_full(:,dim_i)) < 100000000)
% %         piquancy_full_one(dim_i) = max(piquancy_full(:,dim_i)) / 10000000;
% %     else
%         piquancy_full_one(dim_i) = max(piquancy_full(:,dim_i));
% %     end
% end
Pmax_full = max(piquancy_full_one);
piquancy_origin_trick = piquancy_full_one;
dim = length(piquancy_full_one);
for dim_trick1 = 1:dim
    if(piquancy_origin_trick(dim_trick1) < 15000)
        piquancy_origin_trick(dim_trick1) = piquancy_origin_trick(dim_trick1) / 500000;
    end
end
piquancy_origin_trick_max = max(piquancy_origin_trick);
piquancy_origin_trick_db = 10*log10(piquancy_origin_trick./piquancy_origin_trick_max);% normalization

Pmax_db_full = 10*log10(piquancy_full_one/Pmax_full);% normalization

figure();
for dim_i = 1:dim
    for dim_j = 1:dim
        if(piquancy_full(dim_i,dim_j) < 10000000)
             piquancy_full_plot(dim_i,dim_j) = piquancy_full(dim_i,dim_j) / 100000;
        else
            piquancy_full_plot(dim_i,dim_j) = piquancy_full(dim_i,dim_j);
        end
    end
end
Pmax_full_db = max(max(piquancy_full_plot));
Pmax_db_full_db = 10*log10(piquancy_full_plot./Pmax_full_db);% normalization
surf(dis1_full,dis2_full,Pmax_db_full_db);
shading interp
hold on
xlim([10500 10501]);
ylim([10500 10501]);
xl = xlim;
yl = ylim;
% mesh(zeros(size(theta1)) -xl(1),theta2,Pmax_db');
% mesh(theta1, zeros(size(theta2)) -yl(1),Pmax_db');
plot3(zeros(size(dis1_full)) + xl(1),dis2_full,Pmax_db_full','Color','#0072BD');
plot3(dis1_full, zeros(size(dis2_full)) + yl(1),Pmax_db_full','Color','#D95319');
set(gca,'FontName','Times New Roman','FontSize',30,'LineWidth',1);
xlabel('distance1/(degree)','FontSize',40);
ylabel('distance2/(degree)','FontSize',40);
zlabel('piquancy func/(dB)','FontSize',40);
xlim([9500 10500]);
ylim([9500 10500]);