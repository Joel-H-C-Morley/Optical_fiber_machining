function [slice_x, slice_y, I1, surf_offset] = Fiber_reconstruction(lambd,lam_step, z_start_pos, AtCube, RCheckBox, GCheckBox, BCheckBox)
%% Loop for data acquisition
% AtCube.move_z(z_start_pos+(lam_step*lambd));
%z_start = z_start_pos+(lam_step*lambd); %start position in mm
z_start = z_start_pos+(lam_step*lambd/2); %start position in mm
stepsize = lambd/8; %stepsize in mm
n=6; %number of frames
range = (n-1)*stepsize; %range in mm
z_list = linspace(z_start,z_start+range,n);
    
for k = 1:6 %loop to move sample in z
    c = z_list(k);
    AtCube.move_z(c);
    z=AtCube.getPosition_z();
    for i=1:4 %number of frames to acquire and avearge
        ImBG_autoCont = grabImage(RCheckBox, GCheckBox, BCheckBox);
    end
    MeanI(:,:,k) = mean(ImBG_autoCont,3);   
end
Im_opt_cont = im2double(rescale(MeanI,0,256));
AtCube.move_z(z_start);

%% Crop image around fibre
image_size = 500; %Choose cropping area
T = 548-(image_size/2);
B = 548+(image_size/2);
L = 968-(image_size/2);
R = 968+(image_size/2);
IR_crop = Im_opt_cont(T:B,L:R,:); %crop and isolate fibre face
sig = 5;
FilterSize = 21;
IR_filt=imgaussfilt(uint8(mean(IR_crop,3)),sig,'FilterSize',FilterSize); %% Filter
Im_bw = imfill(imbinarize(IR_filt,147/255),'holes'); %Create binary image to use centroiding
%% Create phase map
I1=IR_crop(:,:,1);
I2=IR_crop(:,:,2);
I3=IR_crop(:,:,3);
I4=IR_crop(:,:,4);
I5=IR_crop(:,:,5);
I6=IR_crop(:,:,6);
phi2 = -atan(((3*I2-4*I4+I6)./(I1-4*I3+3*I5))); %calculate phase map

%% Unwrap and create image of surface with no background
Frame = ones(size(phi2,1),size(phi2,2));
[PU,PC,N_unwrap,t_unwrap]=CPULSI(2*phi2,Frame,100,0.1,250,250,false);
%[PU,PC,N_unwrap,t_unwrap]=CPULSI(1*phi2,Frame,100,0.1,250,250,false);
IR_NaN = double(Im_bw);%convert binary image to double
IR_NaN(IR_NaN<1)=nan;%set 0s to NaN
sig = 3;
FilterSize = 11;
p_unwr_filt=flip(imgaussfilt(0.5*PU.*IR_NaN,sig,'FilterSize',FilterSize)); %% Filter
%p_unwr_filt=imgaussfilt(PU.*IR_NaN,sig,'FilterSize',FilterSize); %% Filter
%% 3D plot
pixel=0.329; % Pixel FOV in um
[X,Y]=meshgrid(-image_size/2:image_size/2,-image_size/2:image_size/2);   % Generate 2D meshgrid full
x=X*pixel; %in um                   
y=Y*pixel; %in um
%surface = (p_unwr_filt./(2.*pi)).*lambd; % convert phase to height in mm
surface = (p_unwr_filt./(2.*pi)).*lambd/2; % convert phase to height in mm
surf_offset = surface-min(min(surface));
fig_surf = newfig('Fiber Surface');
        set(gcf,'Position',[1250 250 600 400])
        plt_fiber_surf = findobj(fig_surf, 'type', 'axes');
        if isempty(plt_fiber_surf)
            plt_fiber_surf = axes(fig_surf);
        end
surf(plt_fiber_surf,x,y,1000*surf_offset) % Show distance calibrated results
daspect([1 1 0.05])
shading interp
% axis equal
colorbar;
xlabel('um')
ylabel('um')
zlabel('um')
data_line_x = 1000*surf_offset(:,250);
data_line_y = 1000*surf_offset(250,:);
data_line_x(isnan(data_line_x)) = 0;%set NaNs to 0 in cross section
data_line_y(isnan(data_line_y)) = 0;%set NaNs to 0 in cross sectionslice = [x(1,:)',data_line'];
slice_x = [x(1,:)' data_line_x];
slice_y = [x(1,:)' data_line_y'];
end
