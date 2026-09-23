% 准备数据
x = [58; 63; 68; 72; 77; 81; 86; 91; 95; 100; 105; 109; 114; 118; 123; 128];
y = [23.62;31.63;42.53;48.81;62.14;71.36;98.07;139.11;204.22;369.88;762.16;1487.49;3485.89;7359.20;21271.08;64823.46];
% 构建模型
f = fittype(' a1*exp(b1*x + e1) + a2*x*exp(b2*x+e2) '); %
% 拟合数据
%triple sieve: 2^0.367n-37.15; database size: 3.2 * 2^0.2075n
start_point = [0.05, 0.1,0.2544, 0.1438, 0, 0]; % 初始拟合参数 , 0.292*log(2), 0.208*log(2)
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower',[0,0,0.2544, 0.1438,-inf,-inf], 'Upper', [1,1,0.2544, 0.1438,inf,inf], 'StartPoint', start_point); % 拟合参数选项
[SimTpump, gof] = fit(x, y, f, fit_options); % 进行拟合
% 输出拟合结果
disp(SimTpump);
% 绘制拟合曲线和原始数据
plot(SimTpump, x, y);
xlabel('n'), ylabel('Walltime/sec'), % title('multi-exponential fit');
%a1 =     0.05
%a2 =     0.1
%b1 =     0.2544  (固定在边界处) b1/ln2 = 0.367
%b2 =     0.1438  (固定在边界处) b2/ln2 = 0.2075
%e1 =     -18.79  e1/ln(2) = -27.108
%e2 =     -11.22  e2/ln(2) = -16.187

% 设置文字大小，同时影响坐标轴标注、图例、标题等。
set(gca,'FontSize',14); 

% 调整坐标轴刻度字号大小
set(gca,'fontsize',14);%'linewidth',4,,'fontname','Times'
legend('original data', "SimTpump(n) = 0.05 * 2^{0.367n-27.108}+0.1*n*2^{0.2075n-16.187}",  'Location','NorthWest');

exportgraphics(gcf,'/Users/summer/Documents/GitHub/MyStudy/Artical_in_Public/[2024][Crypto]ProPnjBKZ/5.[20240213]crypto2024/experiment/practical_cost_model_in_dd/SimTpump.png','Resolution',600)