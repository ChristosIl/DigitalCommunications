% ������������ �� ������ ������� ������-����� MPSK, �� ��������
% ������� ��������� ���� �������. ������� ������������
% Gray, ��� ����������� ������� ����������� ��� ����������
% MPSK ���� ��� ������� �������������� ���� Eb/No.
% ������� ���������� �������� ������� ������������ ������� �����
% (���������� � Nyquist), ��� �������� �� BER ���� �� ��������
% ������������� �� �� ���������� ������������� ������.
% M (=2^k, k �����)����� �� ������� ��� ��������� ����������,
% mapping ����� �� �������� �������� ��� � PSK ��������,
% �� ������� ������������� Gray
% Nsymb ����� �� ����� ��� ���������� MPSK,
% nsamp ����� � ����������� ���������������, (# ���������/�)
% EbNo ����� � ����� Eb/No, in db
%
close all; clear all;
k=2; 
M=2^k;
Nsymb=30000;
pulse_type=1; % ������� rtrapezium �� ������� ������������
pulse_type=0; % ������� ������������ ����������� ������
nsamp=16; % ����������� ���������������, # ���������/�
fc=6; % ��������� ��������, ����������� ��� Baud Rate (1/T=2.4)
EbNo=15; % Eb/No, �� db
SNR=EbNo-10*log10(nsamp/k/2); % SNR ��� ������ �������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph1=[pi/4];
theta=[ph1; -ph1; pi-ph1; -pi+ph1];
mapping=exp(1j*theta); % ���������� ������������, M=4
if(k>2)
for j=3:k
theta=theta/2;
mapping= exp(1j*theta);
mapping=[mapping;-conj(mapping)];
theta=angle(mapping);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ������ %%%%%%%%%%%%
x=floor(2*rand(k*Nsymb,1)); % ������ ������� ���������
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb')';
y=[];
for n=1:length(xsym)
y=[y mapping(xsym(n)+1)];
end
%% ������� ������� ������������
if (pulse_type==1) % ������ Nyquist -- rtrapezium
delay = 8; % Group delay (# �������� �)
filtorder = delay*nsamp*2;
rolloff = 0.5; % ����������� ��������� �������
% � rtrapezium.m ������ �� ��������� ��� current directory
shaping_filter = rtrapezium(nsamp,rolloff,delay);
else % ����������� ������
delay=0.5;
shaping_filter=ones(1,nsamp);%/sqrt(nsamp);% �� ��������������!
end
%% ����������� ����
ytx=upsample(y,nsamp);
ytx = conv(ytx,shaping_filter);
 % ����������� & ���������� ��������
figure(1); pwelch(real(ytx),[],[],[],nsamp); % �� ������� 1/T
% quadrature modulation
m=(1:length(ytx));
s=real(ytx.*exp(1j*2*pi*fc*m/nsamp));
figure(2); pwelch(s); % �� ������� 1/T
% ---------------------'
% ��������� ��������, ���� ���. ������.
% if (pulse_type==1) eyediagram(ytx(1:2000),nsamp*2);
% �������� ������ ����������� �������
Ps=10*log10(s*s'/length(s)); % ����� �������, �� db
Pn=Ps-SNR; % ���������� ����� �������, �� db
n=sqrt(10^(Pn/10))*randn(1,length(ytx));
snoisy=s+n; % ��������� ���������� ����
%clear ytx xsym s n; % ��� ������������ ������
%%%%%%%%% ������ %%%%%%%%%%%%%
% �������������
yrx=2*snoisy.*exp(-1j*2*pi*fc*m/nsamp); %clear s;
% �������� ��������� ���������� ������.
% ���� ���� ��������� ������������ Nyquist ��� ����������,
% ���� �� ������������� ������ (Nyquist) ����� ����������
% ��� ��� ���� ���������� ������.
yrx = conv(yrx,shaping_filter);
figure(3); pwelch(real(yrx),[],[],[],nsamp); % ������� 1/T
yrx = downsample(yrx,nsamp); % ������������� ��� ������ nT.
yrx = yrx(2*delay+(1:length(y))); % �������� ����� ���������.
% ----------------------
%yi=real(yrx); yq=imag(yrx); % ��������� ��� �������� ���������
xrx=[]; % �������� �������� ���������� ������ -- ������ ����
q=[0:1:M-1];
for n=1:length(yrx) % ������� ������������ �������
[m,j]=min(abs(angle(mapping)-angle(yrx(n))));
yrx(n)=q(j);
xrx=[xrx; de2bi(q(j),k,'left-msb')'];
end
%scatterplot(yrx); % 
%��������� �������������
% �� ����������� ��� ��������� �������������, ���� ���������
% ������������ �� ���. �����; ���� � �������� ��� ber
% ��� ��� �����������;
%% �������� ������ �����
% �� �������� ��� �������� ���������� (x-xrx)
ber2=sum(not(xrx==x))/length(x);

