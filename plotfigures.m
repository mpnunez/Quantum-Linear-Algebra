function plotfigures

clc
clear

% Accuracy Figure 4
N = [2,4,6,8,10,16,20,26,30,36,40];
N = N + ones(1,length(N));
error = [5.458846195967956,1.758719147601725,0.663230855057197,0.244361437242763,0.080033681185350,8.126202752307687e-04,1.143446196083531e-05,3.197926590203792e-09,4.213518423057394e-12,3.552713678800501e-15,1.532107773982716e-14];

figure
plot(N,log(error),'.')
xlabel('N')
ylabel('log(error)')

% Speed Figure 1
N = [1000
5000
100
500
250
750
2500
];
N = N + ones(length(N),1);
FFT_time = [0.6396041
15.6313002
0.0468003
0.1716011
0.0780005
0.3588023
4.1652267
];
fulleigs = [0.8268053
1.91E+02
0.0312002
0.1404009
0.0624004
0.3276021
21.5905384
];
bottom_six_eigs = [0.2496016
9.4692607
0.0156001
0.0312002
0.0156001
0.1092007
1.8408118
];

N = sort(N);
FFT_time = sort(FFT_time);
fulleigs = sort(fulleigs);
bottom_six_eigs = sort(bottom_six_eigs);

figure
%loglog(N,FFT_time,'.',N,fulleigs,'.',N,bottom_six_eigs,'.')
loglog(N,FFT_time,'-',N,fulleigs,'-',N,bottom_six_eigs,'-')
xlabel('log(N)')
ylabel('log(CPU time)')
legend('FFT','eig','eigs')

% Speed Figure 2

k = [1,2,10,25,50,100,300,500];
CPU = [0.327602099999993,0.327602099999993,0.577203699999998,0.982806300000021,1.778411400000010,2.932818799999993,8.361653599999983,16.941708599999998];

figure
loglog(k,CPU,'.')
xlabel('log(k)')
ylabel('log(CPU time)')

% Accuracy 1



end

