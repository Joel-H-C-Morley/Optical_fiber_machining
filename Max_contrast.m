% function ContrastTest

%% Startup Processes
t0 = tic;
uEye_camera(0); % initializing camera
uEye_camera(8, 2); %exposure time
pause(0.5)

% For attocube
clear AtCube

% C885.turnOnXY()
if ( ~exist ( 'AtCube', 'var' ) || ~isa ( AtCube, 'AtCube_control' ) )
    AtCube = AtCube_control();
else
    error('There exis the var. "AtCube"')
end
AtCube.turnOn()

% start = -1.370; %start position in mm
% % stepsize = 0.0002; %stepsize in mm
% n = 20; %number of frames
% range = 0.03; %range in mm
% z_list = linspace(start,start+range,n);

sample_c = [548,968]; % sample centre, pixels
sample_a = 10; % sample width/height, pixels

% for i=1:n
% [err, M, w, h] = uEye_camera(1); % fetch image
% M(M<0)=M(M<0)+256;
% M = reshape(M, 3, w*h);
% I = reshape(M, 3, w, h);
% I = permute(I, [3 2 1])/255;
% % IG(:,:,i)=I(:,:,2); % take only the green channel
% Im_BG = 0.7*I(:,:,2)+0.3*I(:,:,3);
% % ImBG_autoCont = rescale(Im_BG,0,256);
% Sample = Im_BG(sample_c(1)-0.5*sample_a:sample_c(1)+0.5*sample_a,sample_c(2)-0.5*sample_a:sample_c(2)+0.5*sample_a,:);
% contrast(i) = permute(mean(mean(Sample,1),2),[1 3 2]);
% AtCube.move_z(z_list(i));
% end
% figure
% plot(z_list,contrast)

lambd = 570*10^(-6);
% [peakval,peak] = max((contrast-contrast(1)).^2)
steps = lambd/5; %range in mm
start = AtCube.getPosition_z(); %start position in mm
n = 40; %number of frames
z_list = linspace(start,start+(steps*n),n);

for i=1:n
[err, M, w, h] = uEye_camera(1); % fetch image
M(M<0)=M(M<0)+256;
M = reshape(M, 3, w*h);
I = reshape(M, 3, w, h);
I = permute(I, [3 2 1])/255;
% IG(:,:,i)=I(:,:,2); % take only the green channel
Im_BG = 0.7*I(:,:,2)+0.3*I(:,:,3);
% ImBG_autoCont = rescale(Im_BG,0,256);
Sample = Im_BG(sample_c(1)-0.5*sample_a:sample_c(1)+0.5*sample_a,sample_c(2)-0.5*sample_a:sample_c(2)+0.5*sample_a,:);
contrast(i) = permute(mean(mean(Sample,1),2),[1 3 2]);
AtCube.move_z(z_list(i));
end

figure
plot(z_list,contrast,z_list,movmean(abs(contrast-contrast(1)),7))
% %%%------------------ Fit to moving average  ------------------%%%%

% mov_av = movmean(abs(contrast-contrast(1)),7);
% Gauss = @(ps,x) ps(1)*exp(-((x-ps(2))/ps(3)).^2/2);
% Int = [max(mov_av),start+(steps*n/2),0.0005];
% [params,resnorm,residual,exitflag,output]=lsqcurvefit(Gauss,...
%                                                     Int,...
%                                                     z_list,...
%                                                     mov_av,...
%                                                     [0.9*Int(1),1.1*Int(2),0.9*Int(3)],...
%                                                     [1.1*Int(1),0.9*Int(2),1.1*Int(3)]);
% figure
% plot(z_list,mov_av,'+',z_list,Gauss(params,z_list))
% 
% AtCube.move_z(params(2)-4*lambd);
% uEye_camera(4);
% AtCube.delete()

%%%------------------ Fit to cos func ------------------%%%%
figure
[M,In] = max(contrast);
z_focus = z_list(In-10:In+10);
contrast_focus = contrast(In-10:In+10);
plot(z_focus, contrast_focus,'+')
int = @(c,x) cos(2*pi.*(x-c)./(570e-6/2));
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    [cf,resnorm,residual,exitflag,output]=lsqcurvefit(int,...
                                                    [z_list(In)],...
                                                    z_focus,...
                                                    normalize(contrast_focus),...
                                                    [1.1*z_list(In)],...
                                                    [0.9*z_list(In)],...
                                                    options);
    plot(z_focus,int(cf,z_focus),z_focus,normalize(contrast_focus));

    
%%%------------------ Interference Model ------------------%%%%

div=linspace(-0.5,0.5,n)';
LArray=@(l) arrayfun(@(div0) l(1)+l(2)*div0,div); % v(1)= average l, v(2) = velocity spread dl
int = @(l,x) trapz(LArray(l),cos(l(4)+2*pi.*(x-l(3))./(LArray(l)/2))); % Intensity signal
ampM = @(l,x) (mean(maxk(int(l,x),3))-mean(mink(int(l,x),3))) % calculate amplitude of model
finalint = @(l,x) int(l,x)/max(int(l,x)); % Normalise model
finalint = @(l,x) offset+(ampD/ampM(l,x)*int(l,x)); % Scale model to match data

%%%------------------------ Fitting -----------------------------%%%% 

LInitB = [550e-6,85e-6,0.10024,1*pi]; % Initial values
LInitG = [510e-6,75e-6,0.0995,1*pi];
LInitR = [585e-6,70e-6,0.1,0*pi];
LInitM = [550e-6,100e-6,0.1,1*pi]

TFR = endsWith(Filename,'R.tif')
TFG = endsWith(Filename,'G.tif')
TFB = endsWith(Filename,'B.tif')

if TFR == true
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    [lR,resnorm,residual,exitflag,output]=lsqcurvefit(finalint,...
                                                    LInitR,...
                                                    z_list,...
                                                    contrast,...
                                                    [0.98*LInitR(1),0.9*LInitR(2),0.10015,0*pi],...
                                                    [1.02*LInitR(1),1.1*LInitR(2),0.1004,0.09*pi],...
                                                    options);
    plot(z_list,finalint(lR,z_list),z_list,contrast);
else
    if TFG == true
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        [lG,resnorm,residual,exitflag,output]=lsqcurvefit(finalint,...
                                                    LInitG,...
                                                    z_list,...
                                                    contrast,...
                                                    [0.9*LInitG(1),0.9*LInitG(2),0.0990,0],...
                                                    [1.1*LInitG(1),1.1*LInitG(2),0.1000,2*pi],...
                                                    options);
        plot(z_list,finalint(lG,z_list),z_list,contrast);
    else
        if TFB == true
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        [lB,resnorm,residual,exitflag,output]=lsqcurvefit(finalint,...
                                                    LInitB,...
                                                    z_list,...
                                                    contrast,...
                                                    [0.9*LInitB(1),0.9*LInitB(2),0.0995,0],...
                                                    [1.1*LInitB(1),1.2*LInitB(2),0.1005,2*pi],...
                                                    options);
        plot(z_list,finalint(lB,z_list),z_list,contrast);
        else
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        [lM,resnorm,residual,exitflag,output]=lsqcurvefit(finalint,...
                                                    LInitM,...
                                                    z_list,...
                                                    contrast,...
                                                    [0.9*LInitM(1),0.9*LInitM(2),0.0995,0],...
                                                    [1.1*LInitM(1),1.2*LInitM(2),0.1005,2*pi],...
                                                    options);
        plot(z_list,finalint(lB,z_list),z_list,contrast);
        end
end
end
                                                
hold on
plot(z_list,finalint(l,z_list))
scatter(z_list,finalcontrast)
hold off

dt = toc(t0)
%end