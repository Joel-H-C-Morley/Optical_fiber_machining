function Surface_reconstruction_fcn(lambd,lam_step, z_start_pos, AtCube, RCheckBox, GCheckBox, BCheckBox) % Function to reconstruct th surface, not the same as used for the fiber
    %hdata = guidata(fig);
    %% 6 frame Loop for data acquisition
    %AtCube.move_z(z_start_pos+(lam_step*lambd));
    AtCube.move_z(z_start_pos+(lam_step*lambd/2));
    zstart = AtCube.getPosition_z(); %start position in mm
    stepsize = lambd/8; %stepsize in mm
    n=6; %number of frames
    range = (n-1)*stepsize; %range in mm
    z_list = linspace(zstart,zstart+range,n);
    
    for k = 1:6 %loop to move fiber in z
        c = z_list(k);
        AtCube.move_z(c);
%         z=AtCube.getPosition_z();
            for i=1:2 %number of frames to acquire and avearge
                Im_mono = grabImage(RCheckBox, GCheckBox, BCheckBox);
            end
        MeanI(:,:,k) = mean(Im_mono,3);   
    end
    AtCube.move_z(zstart);
    Im_opt_cont = im2double(rescale(MeanI,0,256));
    %% Crop image around fibre
    image_size = 500; %Choose cropping area
    T = 548-(image_size/2);
    B = 548+(image_size/2);
    L = 968-(image_size/2);
    R = 968+(image_size/2);
    IR_crop = Im_opt_cont(T:B,L:R,:); %crop and isolate fibre face by making background 0

    %% Create phase map
    I1=IR_crop(:,:,1);
    I2=IR_crop(:,:,2);
    I3=IR_crop(:,:,3);
    I4=IR_crop(:,:,4);
    I5=IR_crop(:,:,5);
    I6=IR_crop(:,:,6);
    phi2 = -atan(((3*I2-4*I4+I6)./(I1-4*I3+3*I5))); %calculate phase map
    
    %% Unwrap and create image of surface
    Frame = ones(size(phi2,1),size(phi2,2));
    [PU,~,~,~]=CPULSI(2*phi2,Frame,100,0.0001,250,250,false);
    phase_unwrapped=0.5*PU;%Restore correct phase & remove background
    sig = 3;
    FilterSize = 15;
    Im_filt=flip(im2double(imgaussfilt(phase_unwrapped,sig,'FilterSize',FilterSize))); %% Filter

    %% 3D plot
    pixel=0.329; % Pixel FOV in um
    [X,Y]=meshgrid(-image_size/2:image_size/2,-image_size/2:image_size/2);   % Generate 2D meshgrid full
    x=X*pixel; %in um                   
    y=Y*pixel; %in um
    surface = (Im_filt./(2.*pi)).*lambd/2; % convert phase to height in mm
    surf_offset = surface-min(min(surface));
    fig_surf = newfig('Surface');
        set(gcf,'Position',[1100 320 600 400])
        plt_surf = findobj(fig_surf, 'type', 'axes');
        if isempty(plt_surf)
            plt_surf = axes(fig_surf);
        end
    surf(plt_surf,x,y,1000*surf_offset) % Show distance calibrated results
    daspect([1 1 0.05])
    shading interp
    % axis equal
    colorbar;
    xlabel('um')
    ylabel('um')
    zlabel('um')
    slice = [x(1,:)',1000*surf_offset(:,250)];
    try 
        curv_rad(slice);
    catch ME
        disp('Could not fit curvature')
    end

end