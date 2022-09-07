function [r, flat] = curv_rad(Sample)
G = @(ps,x) ps(1)*exp(-((x-ps(2))/ps(3)).^2/2)+ps(4); % Gaussian function to fit the data to

[pks,locs,w,p] = findpeaks(-(Sample(:,2)-max(Sample(:,2))),... % find position and width of milling shape
                            'MinPeakDistance',250,...
                            'MinPeakProminence',0.1*max(Sample(:,2)));

if size(pks,1) == 0 % if reconstructing a flat surface, don't proceed with fitting a curve
    flat = 1;
    fig_slice = newfig('Cross Section');
        set(gcf,'Position',[750 610 1150 400])
        plt = findobj(fig_slice, 'type', 'axes');
        if isempty(plt)
            plt = axes(fig_slice);
        end
    disp('No peaks found');
    hold on
    plot(plt,Sample(:,1),Sample(:,2))
    hold off

xlabel('x (um)')
ylabel('Height (um)')
else % if peaks are found crop data to focus on the central part of depression (1.5*sigma
flat = 0;
sigma_idx = round(w/2);
data_fit= Sample(round(locs-1.5*sigma_idx):locs+round(1.5*sigma_idx),2);
x_fit = Sample(locs-1.5*sigma_idx:locs+1.5*sigma_idx,1);
    
%% fit a Gaussian
% Estimate parameters
x0 = Sample(locs);
A = -p;
% sigma = Sample(locs(2),1)-x0;
sigma = (0.329*w)/2;
bg = max(Sample(:,2));
 
Init_params = [A x0 sigma bg];%  initial parameters for A, x0, sigma,bg
[fparams,resnorm,residual]=lsqcurvefit(G,...
                                       Init_params,...
                                       x_fit,...
                                       data_fit,...
                                       [2.7*Init_params(1),-8,0.7*Init_params(3),0.8*Init_params(4)],...
                                       [0.8*Init_params(1),8,2*Init_params(3),1.7*Init_params(4)]);
fig_slice = newfig('Cross Section');
        set(gcf,'Position',[750 620 1150 400])
        plt_cross_section = findobj(fig_slice, 'type', 'axes');
        if isempty(plt_cross_section)
            plt_cross_section = axes(fig_slice);
        end
        hold on
plot(plt_cross_section,Sample(:,1),Sample(:,2),'+',x_fit,G(fparams,x_fit))
hold off

%% draw the equivalent circle or radius
r=-(fparams(3)^2)/fparams(1);%um
circle_cent_x=fparams(2);
circle_cent_y=min(G(fparams,x_fit))+r;
circle_limit_low = (3*pi/2)-0.7*(atan(w/4/r));
circle_limit_high = (3*pi/2)+0.7*(atan(w/4/r));
theta = linspace(circle_limit_low ,circle_limit_high,100);
hold on
plot(plt_cross_section,circle_cent_x + r*cos(theta),circle_cent_y + r*sin(theta),'-')
hold off
text(0.8*circle_cent_x,Sample(200,2),sprintf('r= %.3f um',r))%print radius
xlabel('x (um)')
ylabel('Height (um)')
end
end