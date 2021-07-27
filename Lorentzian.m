%% MvsH_Derivative.m
close all;
clear

%get sample name
prompt = 'Sample: ';
samplename = input(prompt, 's');

%generate file name
s0 = '.dat';
filename = strcat(samplename, s0);

fid = fopen(filename);
linenum = 6;
C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
fclose(fid);

D = [C{:}];
D = [D{:}];
E = D(7:11);

if E == 'MPMS3'
    file_format_date = 1; %(1 = new file format after 2019)
    
elseif E == 'SQUID'
    file_format_date = 2; %(2 = old file format before 2019)
end

clear C D E

if file_format_date == 1 %new file format from after 2019
    table_MvsH = readtable(filename, 'HeaderLines', 40, 'CommentStyle', {'(', ')'});
    
    % get sample mass
    fid = fopen(filename);
    for i = 1:1:24
        tline = fgetl(fid);
    end
    fclose(fid);
    expression = ',';
    splitstring = regexp(tline, expression, 'split');
    massstr = cell2mat(splitstring(2));
    mass = str2num(massstr);
    mass = mass./1000; %convert mg to g

    % get sample MW
    fid = fopen(filename);
    for i = 1:1:26
        tline = fgetl(fid);
    end
    fclose(fid);
    expression = ',';
    splitstring = regexp(tline, expression, 'split');
    MWstr = cell2mat(splitstring(2));
    MW = str2num(MWstr);


elseif file_format_date == 2 %old file format from before 2019
table_MvsH = readtable(filename, 'HeaderLines', 27, 'CommentStyle', {'(', ')'});

% get sample mass
fid = fopen(filename);
for i = 1:1:11
    tline = fgetl(fid);
end
fclose(fid);
expression = ',';
splitstring = regexp(tline, expression, 'split');
massstr = cell2mat(splitstring(2));
mass = str2num(massstr);
mass = mass./1000; %convert mg to g

% get sample MW
fid = fopen(filename);
for i = 1:1:13
    tline = fgetl(fid);
end
fclose(fid);
expression = ',';
splitstring = regexp(tline, expression, 'split');
MWstr = cell2mat(splitstring(2));
MW = str2num(MWstr);

end

% calculate moles
mol = mass./MW;

% determine measurement mode (1 = VSM, 2 = DC)
if isnan(table_MvsH.DCMomentFreeCtr_emu_(1)) 
     M2_emu_g_ = table_MvsH.Moment_emu_./mol./5585; %convert to bohr magneton
     m_mode = 1;
else 
     M2_emu_g_ = table_MvsH.DCMomentFreeCtr_emu_./mass;
     m_mode = 2;
end
   
% collect variables in one table
table_MvsH_short = table(round(table_MvsH.Temperature_K_, 2, 'significant'), table_MvsH.MagneticField_Oe_./10000, M2_emu_g_, 'VariableNames', {'Temperature_K_', 'MagneticField_Oe_', 'M_emu_g_'});

% identifies unique temps
temps = unique(table_MvsH_short.Temperature_K_);

% prompt user for temperature
fprintf('Available Temps: \n');
disp(temps)

prompt = 'Which temperature do you want to fit? ';
User_temp = input(prompt);

% prompt user for forward or reverse scan
prompt = ('Forward (1) or reverse (2) scan? (input 1 or 2): ');
scan_direction = input(prompt);

% save scan direction for later as sd_text
if scan_direction == 1
    sd_text = '_Forward';
else
    sd_text = '_Reverse';
end

% select data for given temperature
data_MvsH_short = table2array(table_MvsH_short);
firstColumn = data_MvsH_short(:,1);
data_MvsH_short = data_MvsH_short(firstColumn==User_temp, :);
table_MvsH_short = array2table(data_MvsH_short, 'VariableNames', {'Temperature_K_', 'MagneticField_Oe_', 'M_emu_g_'});

% change bounds depending on measurement mode
if m_mode == 1
  loop_start = 1;
    
    % deletes repeated points so the interpolation works
    [C,ia] = unique(table_MvsH_short.MagneticField_Oe_, 'stable');
    B = table_MvsH_short(ia,:);
    table_MvsH_short = B;
    loop_end = height(table_MvsH_short);
    
    clear B C ia;
    % delete last point if odd number of points so the forward = reverse
    if mod(loop_end,2) == 1 
    data_MvsH_short = table2array(table_MvsH_short);
    data_MvsH_short(end,:) = [];
    table_MvsH_short = array2table(data_MvsH_short, 'VariableNames', {'Temperature_K_', 'MagneticField_Oe_', 'M_emu_g_'});
    loop_end = height(table_MvsH_short);
    end
    
    % forward
    if scan_direction == 1
        derivative_start = (round(height(table_MvsH_short) ./ 2)) + 1 ;
        derivative_end = round(height(table_MvsH_short));
    % reverse
    elseif scan_direction == 2
        derivative_start = 1;
        derivative_end = (round(height(table_MvsH_short) ./ 2));
    end
  
  % for VSM, smooth MvsH with moving median and window of 100
  table_MvsH_short.M_emu_g_ = smoothdata(table_MvsH_short.M_emu_g_, 'movmedian', 100);
  
% DC mode  
elseif m_mode == 2
  loop_start = 1;
  loop_end = height(table_MvsH_short);
  
    if scan_direction == 1
      % forward
        derivative_start = 1;
        derivative_end = round(height(table_MvsH_short) ./ 2);
    elseif scan_direction == 2
      % reverse
        derivative_start = round(height(table_MvsH_short) ./ 2);
        derivative_end = loop_end;
    end
    
  % if DC mode, repeat the -7T point to make the forward/reverse scans equal in length (this is based on Kyle's and Ben's sequence)
  data_MvsH_short = table2array(table_MvsH_short);
  midpoint = data_MvsH_short(round(loop_end./2),:);
  data_MvsH_short = [data_MvsH_short(1:round(loop_end./2), :); midpoint; data_MvsH_short(round(loop_end./2)+1:end, :)];
  table_MvsH_short = array2table(data_MvsH_short, 'VariableNames', {'Temperature_K_', 'MagneticField_Oe_', 'M_emu_g_'});
  
  loop_end = height(table_MvsH_short);
  
  if scan_direction == 1
    derivative_end = (round(height(table_MvsH_short) ./ 2));
  elseif scan_direction == 2
    derivative_end = loop_end;
  end
end
    
% interpolation for MvsH loop forward
H_equalspacing_forward = (linspace(table_MvsH_short.MagneticField_Oe_(1),table_MvsH_short.MagneticField_Oe_(loop_end/2), loop_end/2))' ;
M_interpolated_forward = interp1(table_MvsH_short.MagneticField_Oe_(1:loop_end/2), table_MvsH_short.M_emu_g_(1:loop_end/2), H_equalspacing_forward);

% interpolation for MvsH loop reverse
H_equalspacing_reverse = (linspace(table_MvsH_short.MagneticField_Oe_((loop_end/2)+1),table_MvsH_short.MagneticField_Oe_(loop_end), loop_end/2))' ;
M_interpolated_reverse = interp1(table_MvsH_short.MagneticField_Oe_((loop_end/2)+1:loop_end), table_MvsH_short.M_emu_g_((loop_end/2)+1:loop_end), H_equalspacing_reverse);

% add the forward and reverse H and M values into the table
table_MvsH_short.H_equalspacing = vertcat(H_equalspacing_forward,H_equalspacing_reverse);
table_MvsH_short.M_interpolated = vertcat(M_interpolated_forward, M_interpolated_reverse);

% differentiate
H = table_MvsH_short.H_equalspacing(derivative_start:derivative_end);
M = table_MvsH_short.M_interpolated(derivative_start:derivative_end);

dM = gradient(M)./gradient(H);

%% Fitting
% Fit with 1, 2, or 3 peaks

% prepare data
[xData, yData] = prepareCurveData(H, dM);

% get number of peaks
prompt = 'How many peaks do you want to fit? ';
Num_peaks = input(prompt);

% fit parameters (bounds and starting point - must adjust start points)
if Num_peaks == 1
    ft = fittype( 'y_01 + ((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2))', 'independent', 'x', 'dependent', 'y' );
    
    % fit options
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.MaxFunEvals = 6000;
    opts.MaxIter = 4000;
    
    % A (0, inf), w (0, inf), x (-7 T, 7 T), y (0 inf)
    opts.Lower = [0 0 -7 -Inf];
    opts.Upper = [Inf Inf 7 Inf];
    opts.StartPoint = [50 0.25 0 0];
    
elseif Num_peaks == 2

    ft = fittype( 'y_01 + ((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2))+ y_02 + ((2*A_2)/pi)*(w_2/(4*(x-x_c2)^2+(w_2)^2))', 'independent', 'x', 'dependent', 'y' );
    
    % fit options
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.MaxFunEvals = 6000;
    opts.MaxIter = 4000;
    
    % A1, A2, W1, W2, Xc1, Xc2, y1, y2
    opts.Lower = [0 0 0 0 -7 -7 -Inf -Inf];
    opts.Upper = [Inf Inf Inf Inf 7 7 Inf Inf];
    opts.StartPoint = [50 50 0.25 0.25 0 0.5 0 0];

elseif Num_peaks == 3
    ft = fittype('y_01 + ((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2))+ y_02 + ((2*A_2)/pi)*(w_2/(4*(x-x_c2)^2+(w_2)^2))+ y_03 + ((2*A_3)/pi)*(w_3/(4*(x-x_c3)^2+(w_3)^2))', 'independent', 'x', 'dependent', 'y' );
    
    % fit options
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.MaxFunEvals = 6000;
    opts.MaxIter = 4000;
    
    % A1, A2, A3, W1, W2, W3, Xc1, Xc2, Xc3, y1, y2, y3
    opts.Lower = [0 0 0 0 0 0 -7 -7 -7 -Inf -Inf -Inf];
    opts.Upper = [Inf Inf Inf Inf Inf Inf 7 7 7 Inf Inf Inf];
    opts.StartPoint = [50 50 50 0.25 0.25 0.25 0 0.5 2 0 0 0];
end

% fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% grab coefficients for later
coeffs = coeffvalues(fitresult);

% manipulating the temperature to save within file name and in legend
temp_text_file = strcat(num2str(User_temp) ,'K');
temp_text_legend = strcat(num2str(User_temp) ,' K');

% save fitresult as a .mat file
fit_data = strcat(samplename, '_', temp_text_file, sd_text,  '.mat');
save(fit_data,'fitresult','ft', 'gof', 'coeffs')

% saving coeffs for plotting individual peaks
if Num_peaks == 1
    A_1 = coeffs(1);
    w_1 = coeffs(2); 
    x_c1 = coeffs(3);
    y_01 = coeffs(4);
    
elseif Num_peaks == 2
    A_1 = coeffs(1);
    A_2 = coeffs(2);

    w_1 = coeffs(3);
    w_2 = coeffs(4); 

    x_c1 = coeffs(5);
    x_c2 = coeffs(6);

    y_01 = coeffs(7);
    y_02 = coeffs(8);

elseif Num_peaks == 3
    A_1 = coeffs(1);
    A_2 = coeffs(2);
    A_3 = coeffs(3);

    w_1 = coeffs(4);
    w_2 = coeffs(5); 
    w_3 = coeffs(6);

    x_c1 = coeffs(7);
    x_c2 = coeffs(8);
    x_c3 = coeffs(9);

    y_01 = coeffs(10);
    y_02 = coeffs(11);
    y_03 = coeffs(12);
end

% supress warning for header names
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

%% Figure 1
% Only MvsH
yyaxis left
plot(table_MvsH_short.MagneticField_Oe_ , table_MvsH_short.M_emu_g_, 'Linewidth', 2.5)

% automatically adjust y bounds
ylim ([min(table_MvsH_short.M_interpolated)-max(table_MvsH_short.M_interpolated)*0.1, max(table_MvsH_short.M_interpolated+max(table_MvsH_short.M_interpolated)*0.1)])

% create ylabel
if m_mode == 1
    ylabel({'{\itM} (\mu_B mol^{-1})'});
elseif m_mode == 2
    ylabel({'{\itM} (emu g^{-1})'});
end
    
% create xlabel
xlabel({'{\itH} (T)'});

% specify figure properties
ax = gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.Box = 'on';
ax.LineWidth = 1.5;

yyaxis right
ax.YTickLabel = []; 

legend (temp_text_legend, 'Location',[0.2 0.6 0.12 0.075],'EdgeColor','none');
print(gcf, strcat(samplename,'_', temp_text_file, sd_text, '_MvsH'), '-dpng')
savefig([samplename,'_', temp_text_file, sd_text '_MvsH'  '.fig'])

%% Figure 2
% MvsH curve w/ derivative and fit
figure(2);
yyaxis left
plot(table_MvsH_short.H_equalspacing, table_MvsH_short.M_interpolated,'LineWidth', 2.5)

% automatically adjust y bounds
ylim ([min(table_MvsH_short.M_interpolated)-max(table_MvsH_short.M_interpolated)*0.1, max(table_MvsH_short.M_interpolated+max(table_MvsH_short.M_interpolated)*0.1)])

% create ylabel
if m_mode == 1
    ylabel({'{\itM} (\mu_B mol^{-1})'});
elseif m_mode == 2
    ylabel({'{\itM} (emu g^{-1})'});
end

% set y bounds for dM plot
yyaxis right

% plot derivative points with lorentzian fit
p = plot(fitresult, xData, yData);
set (p, 'MarkerSize', 10, 'Color', [0.8500, 0.3250, 0.0980]);
hold on

ylim ([min(dM)-max(dM)*0.05, max(dM+max(dM)*0.2)]);

% create ylabel
if m_mode == 1
    ylabel({'{\itM} (\mu_B mol^{-1})'});
elseif m_mode == 2
    ylabel({'{\itM} (emu g^{-1})'});
end

% create xlabel
xlabel({'{\itH} (T)'});

% specify figure properties
ax = gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.Box = 'on';
ax.LineWidth = 1.5;

legend (temp_text_legend,'Derivative', 'Lorentzian Fit', 'Location',[0.2 0.6 0.2 0.1],'EdgeColor','none');
print(gcf, strcat(samplename,'_', temp_text_file, sd_text, '_MvsH_Fit'), '-dpng')
savefig([samplename,'_', temp_text_file, sd_text '_MvsH_Fit'  '.fig'])

%% Figure 3
% MvsH curve w/ individual peaks contributing to derivative 
figure(3);
yyaxis left
plot(table_MvsH_short.H_equalspacing, table_MvsH_short.M_interpolated,'LineWidth', 2.5)
hold on

% automatically adjust y bounds
ylim ([min(table_MvsH_short.M_interpolated)-max(table_MvsH_short.M_interpolated)*0.1, max(table_MvsH_short.M_interpolated+max(table_MvsH_short.M_interpolated)*0.1)])

% create ylabel
if m_mode == 1
    ylabel({'{\itM} (\mu_B mol^{-1})'});
elseif m_mode == 2
    ylabel({'{\itM} (emu g^{-1})'});
end

% set y bounds for dM plot
yyaxis right

xlim([-7, 7]);
p_2 = plot(fitresult);
set (p_2, 'Color', [0.8500, 0.3250, 0.0980]);
hold on

xlim([-8, 8]);
ylim ([min(dM)-max(dM)*0.05, max(dM+max(dM)*0.2)]);

% create ylabel
if m_mode == 1
    ylabel({'{\it\chi} (\mu_B mol^{-1} T^{-1})'});
elseif m_mode == 2
    ylabel({'{\it\chi} (emu g^{-1} T^{-1})'});
end

% create xlabel
xlabel({'{\itH} (T)'});

% specify figure properties
ax = gca;
ax.FontName = 'Arial';
ax.FontWeight = 'Bold';
ax.FontSize = 13;
ax.Box = 'on';
ax.LineWidth = 1.5;

% need this to plot the lorentzian peak as a function w/ coeffs in it
syms x

%plotting individual peaks
if Num_peaks == 1
    fplot(y_01 + ((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2)),[-7 7],'--m', 'LineWidth', 1);
    
elseif Num_peaks == 2
    fplot(((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2)),[-7 7],'--m', 'LineWidth', 1);
    hold on
    fplot(((2*A_2)/pi)*(w_2/(4*(x-x_c2)^2+(w_2)^2)),[-7 7], '--','Color', [0 0.6 0.4], 'LineWidth', 1);
    
elseif Num_peaks == 3
    fplot(((2*A_1)/pi)*(w_1/(4*(x-x_c1)^2+(w_1)^2)),[-7 7],'--m', 'LineWidth', 1);
    hold on
    fplot(((2*A_2)/pi)*(w_2/(4*(x-x_c2)^2+(w_2)^2)),[-7 7], '--','Color', [0 0.6 0.4], 'LineWidth', 1);
    hold on
    fplot(((2*A_3)/pi)*(w_3/(4*(x-x_c3)^2+(w_3)^2)),[-7 7], '--c', 'LineWidth', 1);
end

legend (temp_text_legend,'Lorentzian Fit','Peak 1','Peak 2', 'Peak 3', 'Location',[0.2 0.6 0.2 0.1],'EdgeColor','none');
print(gcf, strcat(samplename,'_', temp_text_file, sd_text, '_MvsH_Peaks'), '-dpng')
savefig([samplename,'_', temp_text_file, sd_text '_MvsH_Peaks'  '.fig'])

clearvars ans ax expression fid firstColumn H_equalspacing_forward H_equalspacing_reverse i M_interpolated_forward...
          M_interpolated_reverse massstr midpoint MWstr p p_2 prompt s0 splitstring tline temp_text_file temp_text_legend User_temp...
          derivative_end derivative_start