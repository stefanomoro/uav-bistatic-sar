%% Windowed dopppler
dop_win_size = 2^nextpow2(round(1 / PRI));     % window time of 1s at least
win_idx = 1:dop_win_size;
N_cycle = floor(size(RC_Df_fixed,2)/dop_win_size)*2-1;

slow_freq = linspace(-PRF/2,PRF/2,dop_win_size);
speed_ax = slow_freq * lambda/2 * 3.6;
x_idx = and(speed_ax<20,speed_ax>-20);
y_idx = and(R_ax>=-10,R_ax<100);

x = speed_ax(x_idx);
y = R_ax(y_idx);


fig = figure;
images ={}; 
RC_shifted = circshift(RC_Df_fixed,samp_shift,1);
for i = 1:N_cycle
    RC_doppler = fftshift(fft(RC_shifted(:,win_idx),[],2),2);
%     RC_doppler = 10*log10(RC_doppler);
    for ii = 1:size(RC_doppler,1)
        RC_doppler(ii,:) = RC_doppler(ii,:)./ (sum(abs(RC_Df_fixed(ii,:)))) ;
    end
    AA = abs(RC_doppler(y_idx,x_idx));

    
    win_idx = win_idx + dop_win_size/2;
    imagesc(x,y,AA);
    title("Range-speed"),xlabel("speed [km/h]"),ylabel("Range [m]")
    drawnow
    frame = getframe(fig);
    images{i} = frame2im(frame);
end
close
%% CREATE GIF
filename = strcat(experiment_name,'.gif'); % Specify the output file name
for idx = 1:N_cycle
    [A,map] = rgb2ind(images{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
    end
end