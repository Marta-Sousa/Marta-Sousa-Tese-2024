% Treat SVM data 
clear all
dados = struct();
addpath('Dados')
datapath=uigetdir('Database'); % SVD
%datapath=uigetdir('C/Users/Marta Sousa/Desktop/tese/Tese/github/SVD/todas'); %D1s
%datapath=uigetdir('C/Users/Marta Sousa/Desktop/tese/Tese/github/SVD/todas'); %D2s

d=dir(fullfile(datapath,'*.wav'));

for i=1:numel(d)
  dados(i).name = d(i).name;
  dados(i).record = importdata(fullfile(datapath,d(i).name));
  %dados(i).diagnostic = 1; healthy
  %dados(i).diagnostic = 0; pathological
end

%%
%Calcular o pearson r
%clear all
%load('dados_svd.mat');
feat = readtable(uigetfile('Dados/.csv'));

%Delete rows were f0 is 0
feat(~feat.meanF0Hz,:) = [];


n = height(feat);
r = zeros(n, 1);

for j = 1:n
    
    fs = dados(j).record.fs;
    ts = 1/fs;
    x = dados(j).record.data;
    f0 = feat.meanF0Hz(j);
    delay = round((1/f0)/ts);
    xi = [x; zeros(delay, 1)];
    xk = [zeros(delay, 1); x];     %same recording but delayed
    xi_bar = median(xi); 
    xk_bar = median(xk);

    a = 0;
    b = 0;
    c = 0;

    for i=1:length(xi)
        a = a +(xi(i) - xi_bar)*(xk(i) - xk_bar);
        b = b + (xi(i) - xi_bar)^2;
        c = c + (xi(i) - xi_bar)^2;
    end
    r(j) = a / (sqrt(b*c));
end

r = array2table(r);   

%% Using the Butterworth filter: [2500 3500] and <2500 Hz

clearvars -except dados feat r n

Butter_high = zeros(n, 1); 
Butter_band = zeros(n, 1); 


for j = 1:n  

    fs = dados(j).record.fs;
    ts = 1/fs;
    x = dados(j).record.data;
    f0 = feat.meanF0Hz(j);
    [d_band, n_band] = butter(2, [1000 3000]/8000);
    x_band = filter(d_band, n_band, x);
    
    [d_high, n_high] = butter(2, 2500/8000, 'high');
    x_high = filter(d_high, n_high, x);
    
    delay = round((1/f0)/ts);
    xi_band = [x_band; zeros(delay, 1)]; xk_band = [zeros(delay, 1); x_band];
    xi_high = [x_high; zeros(delay, 1)]; xk_high = [zeros(delay, 1); x_high];     %same recording but delayed
    xi_band_bar = median(xi_band); xk_band_bar = median(xk_band);
    xi_high_bar = median(xi_high); xk_high_bar = median(xk_high);

    [a_band, a_high, b_band, b_high, c_high, c_band] = deal(0,0,0,0,0,0);

    for i=1:length(xi_band)
        a_band = a_band + (xi_band(i) - xi_band_bar)*(xk_band(i) - xk_band_bar);
        b_band = b_band + (xi_band(i) - xi_band_bar)^2;
        c_band = c_band + (xi_band(i) - xi_band_bar)^2;
        
        a_high = a_high + (xi_high(i) - xi_high_bar)*(xk_high(i) - xk_high_bar);
        b_high = b_high + (xi_high(i) - xi_high_bar)^2;
        c_high = c_high + (xi_high(i) - xi_high_bar)^2;
    end

    Butter_high(j) = a_high / (sqrt(b_high*c_high));
    Butter_band(j) = a_band / (sqrt(b_band*c_band));

end

Butter_high = array2table(Butter_high);
Butter_band = array2table(Butter_band); 
r = [r Butter_high Butter_band];

%% filtro eliptico - sugestão dos profs

clearvars -except dados feat r n

Ellip_high = zeros(n, 1); 
Ellip_band = zeros(n, 1);
x_band = 0; x_high = 0; 

for j = 1: n
    
    fs = dados(j).record.fs;
    ts = 1/fs;
    x = dados(j).record.data;
    f0 = feat.meanF0Hz(j);
%[b,a] = ellip(n,Rp,Rs,Wp) returns the transfer function coefficients of an nth-order lowpass digital 
%elliptic filter with normalized passband edge frequency Wp. The resulting filter has Rp decibels of 
%peak-to-peak passband ripple and Rs decibels of stopband attenuation down from the peak passband value.
    [d_band, n_band] = ellip(4,3,10, [1000 3000]/8000);
    x_band = filter(d_band, n_band, x);
    
    [d_high, n_high] = ellip(4,3, 10, 2500/8000, 'high');
    x_high = filter(d_high, n_high, x);
    
    delay = round((1/f0)/ts);
    xi_band = [x_band; zeros(delay, 1)]; xk_band = [zeros(delay, 1); x_band];
    xi_high = [x_high; zeros(delay, 1)]; xk_high = [zeros(delay, 1); x_high];     %same recording but delayed
    xi_band_bar = median(xi_band); xk_band_bar = median(xk_band);
    xi_high_bar = median(xi_high); xk_high_bar = median(xk_high);

    [a_band, a_high, b_band, b_high, c_high, c_band] = deal(0,0,0,0,0,0);

    for i=1:length(xi_band)
        a_band = a_band + (xi_band(i) - xi_band_bar)*(xk_band(i) - xk_band_bar);
        b_band = b_band + (xi_band(i) - xi_band_bar)^2;
        c_band = c_band + (xi_band(i) - xi_band_bar)^2;
        
        a_high = a_high + (xi_high(i) - xi_high_bar)*(xk_high(i) - xk_high_bar);
        b_high = b_high + (xi_high(i) - xi_high_bar)^2;
        c_high = c_high + (xi_high(i) - xi_high_bar)^2;
    end

    Ellip_high(j,1) = a_high / (sqrt(b_high*c_high));
    Ellip_band(j,1) = a_band / (sqrt(b_band*c_band));

end

Ellip_high = array2table(Ellip_high);
Ellip_band = array2table(Ellip_band); 

r = [r Ellip_high Ellip_band];


%% Por tudo na variavel com tudo
dados_svd = readtable('dados_svd.csv');
dados_svd.Properties.VariableNames{1} = 'voiceID';
feat = [feat r];
tabel = join(dados_svd(:,{'voiceID','Diagnostics', 'Age','Gender_f', 'Gender_m'}),feat, 'Keys', 'voiceID');
writetable(tabel, 'Dados/dados_svd_mat.txt')






