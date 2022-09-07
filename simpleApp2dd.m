function simpleApp2dd(varargin)
%% Startup Processes

uEye_camera(0); % initializing camera
uEye_camera(8, 5); %exposure time
pause(0.5)

% For PI
clear C885
PI_Controller = PI_GCS_Controller();
PI_Controller.Destroy();
% handles to X, Y, Z motor axes
if ( ~exist ( 'C885', 'var' ) || ~isa ( C885, 'PI_control' ) )
    C885 = PI_control();
else
    error('There exis the var. "C885"')
end
% C885.turnOnXY()

% For attocube
clear AtCube
if ( ~exist ( 'AtCube', 'var' ) || ~isa ( AtCube, 'AtCube_control' ) )
    AtCube = AtCube_control();
else
    error('There exis the var. "AtCube"')
end
AtCube.turnOn()

%FPGA setup
bit_file = 'CO2_pulse_control.bit';% load bit file
mex_ok_interface('open')
mex_ok_interface('configure', bit_file)
mex_ok_interface('swi', 6, 3);% set default levels
mex_ok_interface('uwi');

%initial GUI settings
fig = uifigure;
fig.Color = [0.94 0.94 0.94];
% fig.Position = [100 100 1571 686];
fig.Position = [385 370 1530 678];
fig.Name = 'Fiber Machining';

%% Begin assembling data structure for configuration
hdata = guidata(fig);
setting_file = 'shoot_settings.mat';
load(setting_file, 'hdata');
hdata.misc.z_start_pos = AtCube.getPosition_z();
hdata.misc.grid = 0;
hdata.misc.shoot = 0;

hdata.pulse.laser_width = 0;
hdata.pulse.laser_shutter_delay1 = 2;
hdata.pulse.laser_shutter_delay2 = 2;

guidata(fig,hdata)

% Create RCheckBox
RCheckBox = uicheckbox(fig);
RCheckBox.Text = 'R';
RCheckBox.Position = [961 635 31 37];
RCheckBox.Value = 1;
% Create GCheckBox
GCheckBox = uicheckbox(fig);
GCheckBox.Text = 'G';
GCheckBox.Position = [961 611 31 37];
GCheckBox.Value = 1;
% Create BCheckBox
BCheckBox = uicheckbox(fig);
BCheckBox.Text = 'B';
BCheckBox.Position = [961 585 31 37];
BCheckBox.Value = 1;
Im = grabImage(RCheckBox, GCheckBox, BCheckBox);
hdata.misc.Image = uint8(Im);

% Initialise FPGA
mex_ok_interface('swi', 0, hdata.pulse.shoot_trig_delay*10);
mex_ok_interface('swi', 1, hdata.pulse.trig_pulse_delay*10);
mex_ok_interface('swi', 2, hdata.pulse.laser_shutter_delay1*10);
mex_ok_interface('swi', 3, hdata.pulse.laser_low_time1*10);
mex_ok_interface('swi', 4, hdata.pulse.laser_width*10);
mex_ok_interface('swi', 5, hdata.pulse.end_delay*10);
mex_ok_interface('swi', 7, hdata.pulse.laser_low_time2*10);
mex_ok_interface('swi', 8, hdata.pulse.laser_shutter_delay2*10);
mex_ok_interface('uwi');

% % set initial PWM Duty cyle
% power_v = 0.103*hdata.pow.pwr;
% daqmx_write_scalar_voltage(hdata.pow.daqmx_dev, power_v, hdata.pow.daqmx_Vmin, hdata.pow.daqmx_Vmax, hdata.pow.daqmx_timeout);

%% Callbacks

        % Button pushed function: SaveButton
        function SaveImagecallback(src,event)
            dir = 'C:\Users\co2la\Documents\Data';
            dir2 = 'X:\Data\WG3';
            time = clock;
            year = sprintf('%.0f',time(1,1));
            month = sprintfc('%02d', time(1,2));
            day = sprintfc('%02d', time(1,3));
            hour = sprintf('%.0f',time(1,4));
            min = sprintfc('%02d', time(1,5));
            folder = append(year, month, day);
            file = append(hour, min);

            if ~exist(char(fullfile(dir, folder)), 'dir')
            mkdir(char(fullfile(dir, folder)))
            end
            if ~exist(char(fullfile(dir2, folder)), 'dir')
            mkdir(char(fullfile(dir2, folder)))
            end

            pltsurf = findobj('type','figure','name','Fiber Surface');
            pltcross = findobj('type','figure','name','Cross Section');

            if ~isempty(findobj('type','figure','name','Fiber Surface'))
            imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir, folder, file),'.tif')))
            saveas(pltsurf, char(strcat(fullfile(dir, folder, file),'surf.fig')))
            Surface = hdata.surf;
            save(char(strcat(fullfile(dir, folder, file),'.mat')),'Surface')
            saveas(pltcross, char(strcat(fullfile(dir, folder, file),'cross.fig')))
            else
            imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir, folder, file),'pre.tif')))
            end
            if ~isempty(findobj('type','figure','name','Fiber Surface'))
            imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir2, folder, file),'.tif')))
            saveas(pltsurf, char(strcat(fullfile(dir2, folder, file),'surf.fig')))
            Surface = hdata.surf;
            save(char(strcat(fullfile(dir2, folder, file),'.mat')),'Surface')
            saveas(pltcross, char(strcat(fullfile(dir2, folder, file),'cross.fig')))
            else
            imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir2, folder, file),'pre.tif')))
            end
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(fig, event)
            stop(tm);
            save('shoot_settings.mat', 'hdata');
            mex_ok_interface('close')
            uEye_camera(4);
            C885.delete()
            AtCube.delete()
            close all
            closereq();
        end
    
        % Button pushed function: AutoCentreButton
        function AutoCentreButtonPushed(fig, event)
            %drawcircle(gridax,'Center',[968,548],'Radius',200,'Color','k','FaceAlpha',0.01,'LineWidth',0.1,'MarkerSize',0.1);
            Im_0 = hdata.misc.Image;
            IR_filt=imgaussfilt(Im_0 ,5,'FilterSize',21); %% Filter
            [centers, radii] = imfindcircles(IR_filt, [180 250]);

            if isempty(centers)
                disp('No fiber found')
            else
                if size(centers(1,:))>1
                    disp('Centre may not be correct, try autofocus first')
                end
            x_start = C885.getPosition_x();
            y_start = C885.getPosition_y();

            pixel=0.000329; % Pixel FOV in mm
            diff_x = pixel*(968-centers(1,1));
            diff_y = pixel*(548-centers(1,2));
            C885.move_x(x_start-diff_x)
            C885.move_y(y_start-diff_y)
            end
        end

        % Button pushed function: AutoFocusButton
        function AutoFocusButtonPushed(fig, event)
            stop (tm)
            step_w = 0.0015;
            range_w = 0.045 ;
            focus = AtCube.getPosition_z();
            startz = focus-(range_w/2);
            zdata = linspace(startz,startz+range_w,range_w/step_w);
            AtCube.move_z(startz);

            sample_c = [548,968]; % sample centre, pixels
            sample_a = 10; % sample width/height, pixels

                for i = 1:size(zdata,2)
                    I = grabImage(RCheckBox, GCheckBox, BCheckBox);
                    Sample = I(sample_c(2)-0.5*sample_a:sample_c(2)+0.5*sample_a,sample_c(1)-0.5*sample_a:sample_c(1)+0.5*sample_a);
                    Imcont3(i) = mean(mean(Sample,1),2);
                    AtCube.move_z(zdata(i));
                end
            figure        
            plot(zdata, Imcont3)
            hold on
            plot(zdata,movmean(Imcont3,5))
                          
            [val, Ind] = max(Imcont3-movmean(Imcont3,5));
             disp(zdata(Ind))
            AtCube.move_z(zdata(Ind));
            start(tm)
        end
     
        % Button pushed function:  Max Contrast
        function Find_max_contrast(src,event)
            
            hdata = guidata(fig);
            sample_c = [548,968]; % sample centre, pixels
            sample_a = 10; % sample width/height, pixels
            hdata.misc.lambd = 570*10^(-6);
            steps = hdata.misc.lambd/5; %range in mm
            zstart = AtCube.getPosition_z(); %start position in mm
            n = 40; %number of frames
            z_list = linspace(zstart,zstart+(steps*n),n);

            for i=1:n
                Imag = grabImage(RCheckBox, GCheckBox, BCheckBox);
                Sample = Imag(sample_c(1)-0.5*sample_a:sample_c(1)+0.5*sample_a,sample_c(2)-0.5*sample_a:sample_c(2)+0.5*sample_a,:);
                contrast(i) = permute(mean(mean(Sample,1),2),[1 3 2]);
                AtCube.move_z(z_list(i));
            end
            
            mov_av = movmean(abs(contrast-contrast(1)),7);
            Gauss = @(ps,x) ps(1)*exp(-((x-ps(2))/ps(3)).^2/2);
            Int = [max(mov_av),zstart+(steps*n/2),0.0005];
            [params,resnorm,residual,exitflag,output]=lsqcurvefit(Gauss,...
                                                    Int,...
                                                    z_list,...
                                                    mov_av,...
                                                    [0.9*Int(1),1.1*Int(2),0.9*Int(3)],...
                                                    [1.1*Int(1),0.9*Int(2),1.1*Int(3)]);
            fig_cont = newfig('Contrast');
            set(gcf,'Position',[320 320 400 400])
            plt_cont = findobj(fig_cont, 'type', 'axes');
            if isempty(plt_cont)
                plt_cont = axes(fig_cont);
            end
            plot(plt_cont,z_list,normalize(contrast),z_list,normalize(Gauss(params,z_list)))
            hdata.misc.z_start_pos = params(2);
            disp(hdata.misc.z_start_pos)
%             AtCube.move_z(hdata.misc.z_start_pos);
            hdata.misc.max_contrast_status = 1;
            guidata(fig,hdata)
        end
  
        % Value changed function: StartingpositionfrommaxcontrastSpinner
        function StartingpositionfrommaxcontrastSpinnerValueChanged(fig, event)
            value = StartingpositionfrommaxcontrastSpinner.Value;
            hdata = guidata(fig);
            hdata.misc.lam_step = value;
            guidata(fig,hdata)
        end

        % Button pushed function: Reconstruct Fiber
        function Reconstruct_FiberButtonPushed(src,event)
            stop(tm)
            hdata = guidata(fig);
            if hdata.misc.max_contrast_status == 0
                hdata.misc.z_start_pos = AtCube.getPosition_z();
            end
            
                if TabGroup2.SelectedTab == SurfaceTab
                    Surface_reconstruction_fcn();
                else
                    [slice_x, slice_y, I1, hdata.surf] = Fiber_reconstruction(hdata.misc.lambd, hdata.misc.lam_step, hdata.misc.z_start_pos, AtCube);
                    close(findobj('type','figure','name','Cross Section'))
            try
                    curv_rad(slice_x);
                    curv_rad(slice_y);
            catch ME
            	display('Could not fit or measure curvature')
            end
                    fig_Im = newfig('Image');
                    set(gcf,'Position',[840 250 450 400])
                    Im1 = findobj(fig_Im, 'type', 'axes');
                    if isempty(Im1)
                        Im1 = axes(fig_Im);                      
                    end
                    imagesc(Im1,I1) % Show distance calibrated results
                    axis image %Make sure axes are equal sizes
                end
            guidata(fig,hdata)
            start(tm)
        end
    
        % Value changed function: GridSwitch
        function GridSwitchValueChanged(src, event)
             switch src.Value
                case 'On'
                    gridax.XGrid = 'on';
                    gridax.YGrid = 'on';
                    set(allchild(gridax), 'Visible', 'on')
                    delete(hdata.mill_plot);
                case 'Off'
                    gridax.XGrid = 'off';
                    gridax.YGrid = 'off';
                    set(allchild(gridax), 'Visible', 'off')
            end
        end
    
        % Value changed function: DirectionSwitch
        function DirectionSwitchValueChanged(fig, event)
            hdata = guidata(fig);
            switch DirectionSwitch.Value
                case 'Winter'
                     hdata.mill.direction = 0;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
                case 'Summer'
                     hdata.mill.direction = 1;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
            end 
        end
    
        % Value changed function: StyleSwitch
        function StyleSwitchValueChanged(fig, event)
            hdata = guidata(fig);
            switch StyleSwitch.Value
                case 'Radial'
                     hdata.mill.style = 0;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
                     hdata.mill.style
                case 'C.centric'
                     hdata.mill.style = 1;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
                     hdata.mill.style
            end 
        end

        % Value changed function: RandomSwitch
        function RandomSwitchValueChanged(fig, event)
            hdata = guidata(fig);
            switch RandomSwitch.Value
                case 'Off'
                     hdata.mill.random = 0;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
                case 'On'
                     hdata.mill.random = 1;
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
            end 
        end
               
        % Value changed function: StepSizeEditField
        function StepSizeEditFieldValueChanged(fig, ~)
            value = StepSizeEditField.Value;
            hdata = guidata(fig);
            hdata.misc.step_size = 0.001*value;
            guidata(fig,hdata)
        end
    
        % Button pushed function: m1Button
        function Buttonm1Pushed(fig, event)
            value = StepSizeEditField.Value;
            hdata = guidata(fig);
            f = fix(abs(value) .* 10 .^ (-floor(log10(abs(value)))));
            hdata.misc.step_size = 0.001*(value-(value/f));
            StepSizeEditField.Value = 1000*hdata.misc.step_size;
            guidata(fig,hdata)
        end
    
        % Button pushed function: p1Button
        function Buttonp1Pushed(fig, event)
            value = StepSizeEditField.Value;
            hdata = guidata(fig);
            f = fix(abs(value) .* 10 .^ (-floor(log10(abs(value)))));
            hdata.misc.step_size = 0.001*(value+(value/f));
            StepSizeEditField.Value = 1000*hdata.misc.step_size;
            guidata(fig,hdata)
        end
    
        % Button pushed function: d10Button
        function Buttond10Pushed(fig, ~)
            value = StepSizeEditField.Value;
            hdata = guidata(fig);
            hdata.misc.step_size = 0.001*(0.1*value);
            StepSizeEditField.Value = 1000*hdata.misc.step_size;
            guidata(fig,hdata)
        end
    
        % Button pushed function: x10Button
        function Buttonx10Pushed(fig, event)
            value = StepSizeEditField.Value;
            hdata = guidata(fig);
            hdata.misc.step_size = 0.001*(10*value);
            StepSizeEditField.Value = 1000*hdata.misc.step_size;
            guidata(fig,hdata)
        end
    
        % Button pushed function: UpButton
        function UpButtonPushed(fig, event)
             y_pos = C885.getPosition_y();
             C885.move_y(y_pos-hdata.misc.step_size);
             yEditField.Value = C885.getPosition_y();
        end

        % Button pushed function: DownButton
        function DownButtonPushed(fig, event)
             y_pos = C885.getPosition_y();
             C885.move_y(y_pos+hdata.misc.step_size);
             yEditField.Value = C885.getPosition_y();
        end

        % Button pushed function: RightButton
        function RightButtonPushed(fig, event)
            x_pos = C885.getPosition_x();
            C885.move_x(x_pos+hdata.misc.step_size);
            xEditField.Value = C885.getPosition_x();
        end
    
        % Button pushed function: LeftButton
        function LeftButtonPushed(fig, event)
            x_pos = C885.getPosition_x();
            C885.move_x(x_pos-hdata.misc.step_size);
            xEditField.Value = C885.getPosition_x();
        end

        % Button pushed function: InButton
        function InButtonPushed(fig, event)
            z_pos = AtCube.getPosition_z();
            AtCube.move_z(z_pos+hdata.misc.step_size);
            zEditField.Value = AtCube.getPosition_z();
        end

        % Button pushed function: OutButton
        function OutButtonPushed(fig, event)
            z_pos = AtCube.getPosition_z();
            AtCube.move_z(z_pos-hdata.misc.step_size);
            zEditField.Value = AtCube.getPosition_z();
        end
    
        % Button pushed function: MovetoTargetButton
        function MovetoTargetButtonPushed(fig, event)
            %if ShootLamp.Color = 'green'
             %  disp('Wait for shooting to finish')
            %else
                BSLamp.Color = 'red';
                ShootPosLamp.Color = 'red';
                C885.move_x(hdata.target.x_pos);
                C885.move_y(hdata.target.y_pos);
                AtCube.move_z(hdata.target.z_pos);
                C885.move_z(23.9);
                MovetoBeamButton.Value = 0;
                disp([hdata.target.x_pos, hdata.target.y_pos, hdata.target.z_pos])
            %end
        end
    
        % Button pushed function: SaveTargetposButton
        function SaveTargetposButtonPushed(fig, event)
            hdata = guidata(fig);
            hdata.target.x_pos = C885.getPosition_x();
            hdata.target.y_pos = C885.getPosition_y();
            hdata.target.z_pos = AtCube.getPosition_z();
            disp([hdata.target.x_pos, hdata.target.y_pos, hdata.target.z_pos])
            guidata(fig,hdata);
        end

        % Button pushed function: MovetoFiberButton
        function MovetoFiberButtonPushed(fig, event)
            BSLamp.Color = 'red';
            ShootPosLamp.Color = 'red';
            C885.move_x(hdata.fiber.x_pos);
            C885.move_y(hdata.fiber.y_pos);
            AtCube.move_z(hdata.fiber.z_pos);
            C885.move_z(23.9);
            disp([hdata.fiber.x_pos, hdata.fiber.y_pos, hdata.fiber.z_pos])
            MovetoBeamButton.Value = 0
        end
    
        % Button pushed function: SaveFiberposButton
        function SaveFiberposButtonPushed(fig, event)
            hdata = guidata(fig);
            hdata.fiber.x_pos = C885.getPosition_x();
            hdata.fiber.y_pos = C885.getPosition_y();
            hdata.fiber.z_pos = AtCube.getPosition_z();
            disp([hdata.fiber.x_pos, hdata.fiber.y_pos, hdata.fiber.z_pos])
            guidata(fig,hdata);
        end
    
        % Button pushed function: Move to ShootPositionButton
        function ShootPositioncallback(fig,event)
            C885.move_z(0)
            BSLamp.Color = 'green';
            diff_x = hdata.beam.x_pos-hdata.knife.x_pos;
            diff_y = hdata.beam.y_pos-hdata.knife.y_pos;
            diff_z = hdata.beam.z_pos-hdata.knife.z_pos;
            x_start = C885.getPosition_x();
            y_start = C885.getPosition_y();
            z_start = AtCube.getPosition_z();
            disp([x_start, y_start, z_start])
            C885.move_x(x_start+diff_x);
            C885.move_y(y_start+diff_y);
            AtCube.move_z(z_start+diff_z+hdata.pow.offset);
            %AtCube.move_z(z_start+diff_z);
            ShootPosLamp.Color = 'green';
%             stop(tm);
%             start(pow_tm);
        end
    
        % Value changed function: Zoffset
        function Spin_offsetValueChanged(fig, event, hdata)
            value = Spin_offset.Value;
            hdata = guidata(fig);
            hdata.pow.offset = value;
            disp(hdata.pow.offset)
            guidata(fig,hdata);
        end
    
        % Value changed function: PulseSpinner
        function Spin_pulseValueChanged(fig, event)
            value = Spin_pulse.Value;
            hdata = guidata(fig);
            hdata.pulse.lowtime1 = value/2;
            hdata.pulse.lowtime2 = value/2;
            % set pulse width
            mex_ok_interface('swi', 4, hdata.pulse.laser_width*10);
            mex_ok_interface('uwi');
            guidata(fig,hdata);
        end
       
        % Value changed function: NoshotsSpinner
        function NoshotsSpinnerValueChanged(fig, event)
            value = NoshotsSpinner.Value;
            hdata = guidata(fig);
            hdata.misc.num_shots = value;
            hdata.misc.num_shots
            guidata(fig,hdata);
        end
    
        % Value changed function: NoringsSpinner
        function NoringsSpinnerValueChanged(fig, event, hdata)
            value = NoringsSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.n_rings = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end
    
        % Value changed function: FreqSpinner
        function FreqSpinnerValueChanged(fig, event, hdata)
            value = FreqSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.n_shots = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end
    
        % Value changed function: RotSpinner
        function RotSpinnerValueChanged(fig, event, hdata)
            value = RotationSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.rot = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end
    
        % Value changed function: RotSpinner
        function RadSpinnerValueChanged(fig, event, hdata)
            value = RadiusSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.rad = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end
    
        % Value changed function: AngleSpinner
        function AngleSpinnerValueChanged(fig, event)
            value = AngleSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.ecc_angle = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end

        % Value changed function: RadiiratioSpinner
        function RadiiratioSpinnerValueChanged(fig, event, hdata)
            value = RadiiratioSpinner.Value;
            hdata = guidata(fig);
            hdata.mill.ecc_degree = value;
            guidata(fig,hdata);
            dotmillingschematic();
            plot_milling();
        end

        % Value changed function: beam_x
        function beam_xValueChanged(fig, event)
            value = beam_xEditField.Value
            hdata = guidata(fig);
            hdata.beam.x_pos = value;
            guidata(fig,hdata)
        end
    
        % Value changed function: beam_y
        function beam_yValueChanged(fig, event)
            value = beam_yEditField.Value
            hdata = guidata(fig);
            hdata.beam.y_pos = value;
            guidata(fig,hdata)
        end
    
        % Value changed function: beam_z
        function beam_zValueChanged(fig, event)
            value = beam_zEditField.Value
            hdata = guidata(fig);
            hdata.beam.z_pos = value;
            guidata(fig,hdata)
        end
    
        % Value changed function: Grid Shooting
        function GridCheckBoxValueChanged(fig, event)
            hdata = guidata(fig);
            value = GridCheckBox.Value;
            if value == 0
                hdata.misc.grid = 0;
            else
                hdata.misc.grid = 1;
            end
            guidata(fig,hdata)
        end
           
        % Value changed function: ShootingSwitch
        function ShootingSwitchValueChanged(fig, event)
            hdata = guidata(fig);
            switch ShootingSwitch.Value
                case 'Single'
                     hdata.misc.shoot = 0;
                     guidata(fig,hdata)
                     delete(hdata.mill_plot);
                case 'Milling'
                     hdata.misc.shoot = 1;
                     guidata(fig,hdata)
                     plot_milling();
            end
        end
        
        % Button pushed function: ShootButton
        function Shootcallback(fig,event)
            hdata = guidata(fig);
            hdata.misc.shoot
            if AtCube.getPosition_z() < 3+hdata.beam.z_pos
                if C885.getPosition_z() == 0
                    ShootLamp.Color = 'green';
                    if hdata.misc.grid == 0
                        if hdata.misc.shoot == 0
                            for i=1:hdata.misc.num_shots
                                mex_ok_interface('ati', 64, 1);% activate trigger
                                mex_ok_interface('ati', 64, 2);% activate trigger
                                pause((0.001*hdata.pulse.laser_width)+0.5)
                            end
                            else
                            z_0 = AtCube.getPosition_z();
                            x_shoot = C885.getPosition_x();
                            y_shoot = C885.getPosition_y();
                            %AtCube.move_z(z_0+hdata.pow.offset);
                            for i=1:hdata.misc.num_shots
                                for s = 1:hdata.mill.totalshots
                                    C885.move_x(x_shoot+(hdata.mill.x_dot(s)*0.000329));
                                    C885.move_y(y_shoot+(hdata.mill.y_dot(s)*0.000329));
                                    mex_ok_interface('ati', 64, 1);% activate trigger
                                    mex_ok_interface('ati', 64, 2);% activate trigger
                                    pause((0.001*hdata.pulse.laser_width)+1)
                                end
                            end    
                        end
                    else
                        XPOS = C885.getPosition_x();
                        YPOS = C885.getPosition_y();
                        for i=1:hdata.misc.num_shots
                            i
                            for iy=1:4
                                for ix=1:4
                                    %C885.move_x(XPOS+ix*0.1);
                                    mex_ok_interface('ati', 64, 1);% activate trigger
                                    mex_ok_interface('ati', 64, 2);% activate trigger
                                    pause(0.5)
                                end
                                C885.move_y(YPOS+iy*0.1);
                            end
                        end
                    end
                end
                        
                else
                    disp('System not in shooting configuration')
            end
                ShootLamp.Color = 'red';
                if TabGroup2.SelectedTab == SurfaceTab
                    MovetoTargetButtonPushed()
                else
                    MovetoFiberButtonPushed()
                end
                
        end
                    
%% Helper functions
function ImBG_autoCont = grabImage(RCheckBox, GCheckBox, BCheckBox)
    [err, M, w, h] = uEye_camera(1); % fetch image
        if err>0
            error('Capturing an image failed!')
        end
    M(M<0)=M(M<0)+256; %Restructure incoming pixel values to fall between 0 and 256
    M = reshape(M, 3, w*h); %Create 3D matrix for 2D image in R, G and B
    I = reshape(M, 3, w, h);
    I = permute(I, [3 2 1]);
    RVal = RCheckBox.Value; % Get checkbox values for which colour channels to use
    GVal = GCheckBox.Value;
    BVal = BCheckBox.Value;
    Im_BG = 0.55*RVal*I(:,:,2)+0.35*GVal*I(:,:,3)+0.1*BVal*I(:,:,1); % Create monochrome image from the 3 colour channels
    ImBG_autoCont = rescale(Im_BG,0,256); % Rescale pixel values to maximise contrast, keep floating point precision
end

function updateImage(src, evt)
    hdata = guidata(fig);
    Im_opt_cont = grabImage(RCheckBox, GCheckBox, BCheckBox); % Grab processed image from camera
    hdata.misc.Image = uint8(Im_opt_cont); % Convert pixel values to integers for use as an image. Save image to GUI data to be used for surface reconstruction etc.
    imshow(hdata.misc.Image, 'Parent', Imageax); % Update the GUI axes with new image
    %set(Imageax, 'CData', Im_opt_cont);
    guidata(fig,hdata);
end

% function [FM] = focus_scan(frames,zdata)
%     ImBG_autoCont = ones(1096,1936,frames);
%     for i = 1:size(zdata,2)
%       for j = 1:frames
%         ImBG_autoCont(:,:,j) = grabImage();
%       end
%         Imcont3 = mean(ImBG_autoCont,3);
%         IR_filt=imgaussfilt(Imcont3,4,'FilterSize',11); %% Filter
%         AtCube.move_z(zdata(i));
%         FM(i) = fmeasure(Imcont3,'CONT',[545 968 500 500]);
%     end
% end

function Surface_reconstruction_fcn % Function to reconstruct th surface, not the same as used for the fiber
    hdata = guidata(fig);
    %% 6 frame Loop for data acquisition
    AtCube.move_z(hdata.misc.z_start_pos+(hdata.misc.lam_step*hdata.misc.lambd));
    zstart = AtCube.getPosition_z(); %start position in mm
    stepsize = hdata.misc.lambd/8; %stepsize in mm
    n=6; %number of frames
    range = (n-1)*stepsize; %range in mm
    z_list = linspace(zstart,zstart+range,n);
    
    for k = 1:6 %loop to move fiber in z
        c = z_list(k);
        AtCube.move_z(c);
%         z=AtCube.getPosition_z();
            for i=1:2 %number of frames to acquire and avearge
                Im_BG = grabImage(RCheckBox, GCheckBox, BCheckBox);
            end
        MeanI(:,:,k) = mean(Im_BG,3);   
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
    Im_filt=im2double(imgaussfilt(phase_unwrapped,sig,'FilterSize',FilterSize)); %% Filter

    %% 3D plot
    pixel=0.329; % Pixel FOV in um
    [X,Y]=meshgrid(-image_size/2:image_size/2,-image_size/2:image_size/2);   % Generate 2D meshgrid full
    x=X*pixel; %in um                   
    y=Y*pixel; %in um
    surface = (Im_filt./(2.*pi)).*hdata.misc.lambd; % convert phase to height in mm
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
    guidata(fig,hdata)
end

function [x_dot, y_dot] = dotmillingschematic()
    hdata = guidata(fig);
    fiber_r = hdata.mill.rad/0.329; % Fiber radius in mm (approx)
    r_list = round(linspace(0,fiber_r,2+hdata.mill.n_rings));
    r_list(1) = [];
    r_list(size(r_list,2)) = [];
    deg_of_ecc = hdata.mill.ecc_degree;
    phi = hdata.mill.ecc_angle; % angle of ellipse
    rx = r_list.*deg_of_ecc;
    ry = r_list;
    for r_idx = 1:hdata.mill.n_rings % Generate values for x and y, each loop is for each ring of milling pattern
    r = @(theta) (rx*ry*(r_list(r_idx)^2))./(sqrt(((rx*r_list(r_idx))^2)*sin((phi*pi/180)+theta).^2+((ry*r_list(r_idx))^2)*cos((phi*pi/180)+theta).^2));
    theta(r_idx,:)=linspace(0+(hdata.mill.rot*r_idx),...
                            (2*pi+(hdata.mill.rot*r_idx))-((2*pi+(hdata.mill.rot*r_idx))/hdata.mill.n_shots),...
                            hdata.mill.n_shots);
    dyn_r = sqrt((rx(r_idx).*cos(theta(r_idx,:))).^2+(ry(r_idx).*sin(theta(r_idx,:))).^2);
    xp(r_idx,:)= (rx(r_idx).*cos(theta(r_idx,:)).*cos(phi*pi/180))-(ry(r_idx).*sin(theta(r_idx,:)).*sin(phi*pi/180));
    yp(r_idx,:)= (rx(r_idx).*cos(theta(r_idx,:)).*sin(phi*pi/180))+(ry(r_idx).*sin(theta(r_idx,:)).*cos(phi*pi/180));
    end

    if hdata.mill.style == 0
        x_dot = permute(xp, [1 2]); % Shoot radially
        y_dot = permute(yp, [1 2]);
        else
        x_dot = permute(xp, [2 1]); % Shoot concentrically
        y_dot = permute(yp, [2 1]);
    end
    hdata.mill.x_dot = -x_dot(:); 
    hdata.mill.y_dot = -y_dot(:);
    if hdata.mill.direction == 0 % Determine milling direction
        else
        hdata.mill.x_dot = flip(hdata.mill.x_dot);
        hdata.mill.y_dot = flip(hdata.mill.y_dot);
    end 
    if hdata.mill.random == 0 % Make milling pattern random or maximally separated
        else
        %rand_list = randperm(length(hdata.mill.x_dot));
%         hdata.mill.x_dot = hdata.mill.x_dot(rand_list);
%         hdata.mill.y_dot = hdata.mill.y_dot(rand_list); 

        ring = length(hdata.mill.x_dot)/hdata.mill.n_rings;
        section =ring/3;
        idx = ones(1,length(hdata.mill.x_dot));
        for i = 2:length(hdata.mill.x_dot)
            if idx(i-1)+section > ring
                if i > ring 
                	if i > 2*ring
                        idx(i) = idx(i-ring)+ring;
                        else
                        idx(i) = idx(i-ring)+ring;
                    end
                    else
                    idx(i) = i-(2*((i-1)/3));
                end
                else
                idx(i) = idx(i-1)+section;
            end
        end
%         section = length(hdata.mill.x_dot)/3;
%         for i = 2:length(hdata.mill.x_dot)
%             if idx(i-1)+section > length(hdata.mill.x_dot)
%                 idx(i) = i-(2*((i-1)/3));
%             else
%                 idx(i) = idx(i-1)+section;
%             end
%         end
        hdata.mill.x_dot = hdata.mill.x_dot(idx);
        hdata.mill.y_dot = hdata.mill.y_dot(idx);
    end 
    guidata(fig,hdata);    
end

function plot_milling()
hdata = guidata(fig);
hdata.mill.totalshots = hdata.mill.n_rings*hdata.mill.n_shots;
f = linspace(1,10,hdata.mill.totalshots);
xc=968;
yc=548;
hdata.mill_plot = scatter(gridax, (xc+hdata.mill.x_dot(:))',(yc+hdata.mill.y_dot(:))',[],'+','LineWidth',2,'CData',f);
drawcircle(gridax,'Center',[968,548],'Radius',200,'Color','k','FaceAlpha',0.01,'LineWidth',0.1,'MarkerSize',0.1);
guidata(fig,hdata)
end

%% GUI layout
            % Create Image, must use plot axes for the image in order to also draw the
            % alignment lines
            Imageax = uiaxes(fig);
            gridax = uiaxes(fig);
            Imageax.Position = [314 89 893 511];
            Imageax.XLim = [0 1936]; % Set limits of axes
            Imageax.YLim = [0 1096];
            gridax.Position = [314 89 893 511];
            gridax.XLim = [0 1936]; % Set limits of axes
            gridax.YLim = [0 1096];
            gridax.Color = 'None';
            Imageax.TickLength = [0 0];
            Imageax.XTick = 968;
            Imageax.XTickLabel = '';
            Imageax.YTick = 548;
            Imageax.YTickLabel = '';
            gridax.TickLength = [0 0];
            gridax.XTick = 968;
            gridax.XTickLabel = '';
            gridax.YTick = 548;
            gridax.YTickLabel = '';
            gridax.XGrid = 'on';
            gridax.YGrid = 'on';
            gridax.GridAlpha = 1;
            gridax.GridLineStyle = ':';
            imshow(hdata.misc.Image, 'Parent', Imageax);
            drawcircle(gridax,'Center',[968,548],'Radius',195,'Color','k','FaceAlpha',0.01,'LineWidth',0.1,'MarkerSize',0.1);
            tm = timer('TimerFcn', @updateImage,'Period', 0.2,'ExecutionMode','fixedSpacing','BusyMode','drop');
            start(tm);
               
            % Create SaveButton
            SaveButton = uibutton(fig,'push','ButtonPushedFcn',@SaveImagecallback);
            SaveButton.Text = 'Save';
            SaveButton.Position = [351 603 180 60];
            
            % Create AutoCentreButton
            AutoCentreButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoCentreButtonPushed);
            AutoCentreButton.FontWeight = 'bold';
            AutoCentreButton.FontColor = [0.0745 0.6235 1];
            AutoCentreButton.Position = [573 603 180 60];
            AutoCentreButton.Text = {'Auto-Centre'; ''};
            
%             % Create AutoFocusButton
%             AutoFocusButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoFocusButtonPushed);
%             AutoFocusButton.FontWeight = 'bold';
%             AutoFocusButton.FontColor = [0.0745 0.6235 1];
%             AutoFocusButton.Position = [794 603 180 60];
%             AutoFocusButton.Text = {'Auto-Focus'; ''};
   
            % Create CloseButton
            CloseButton = uibutton(fig, 'push','ButtonPushedFcn', @CloseButtonPushed);
            CloseButton.Position = [1015 603 180 60];
            CloseButton.Text = 'Close';
            
            % Create MaxContrastButton
            MaxContrastButton = uibutton(fig, 'push', 'ButtonPushedFcn', @Find_max_contrast);
            MaxContrastButton.Position = [351 25 180 60];
            MaxContrastButton.Text = 'Max Contrast';
           
            % Create ReconstructFiberButton
            ReconstructFiberButton = uibutton(fig,'ButtonPushedFcn', @Reconstruct_FiberButtonPushed);
            ReconstructFiberButton.FontSize = 20;
            ReconstructFiberButton.Position = [794 25 180 60];
            ReconstructFiberButton.Text = 'Reconstruct Fiber';
            
            % Create StartingpositionfrommaxcontrastSpinnerLabel
            StartingpositionfrommaxcontrastSpinnerLabel = uilabel(fig);
            StartingpositionfrommaxcontrastSpinnerLabel.HorizontalAlignment = 'right';
            StartingpositionfrommaxcontrastSpinnerLabel.Position = [561 33 103 47];
            StartingpositionfrommaxcontrastSpinnerLabel.Text = 'Starting position from max contrast';
            % Create StartingpositionfrommaxcontrastSpinner
            StartingpositionfrommaxcontrastSpinner = uispinner(fig,'ValueChangedFcn',@StartingpositionfrommaxcontrastSpinnerValueChanged);
            StartingpositionfrommaxcontrastSpinner.Limits = [-10 10];
            StartingpositionfrommaxcontrastSpinner.Position = [672 28 67 52];
            StartingpositionfrommaxcontrastSpinner.Value = hdata.misc.lam_step;
            StartingpositionfrommaxcontrastSpinner.Step = 0.5;
            % Create Label
            Label = uilabel(fig);
            Label.FontSize = 16;
            Label.Position = [748 25 47 58];
            Label.Text = 'λ';
            
            % Create GridSwitchLabel
            GridSwitchLabel = uilabel(fig);
            GridSwitchLabel.HorizontalAlignment = 'center';
            GridSwitchLabel.Position = [1067,28,45,20];
            GridSwitchLabel.Text = 'Grid';
            % Create GridSwitch
            GridSwitch = uiswitch(fig, 'slider','ValueChangedFcn',@GridSwitchValueChanged);
            GridSwitch.Position = [1071 53 45 20];
            GridSwitch.Value = 'On';
            
            % Create XYZstagesPanel
            XYZstagesPanel = uipanel(fig);
            XYZstagesPanel.TitlePosition = 'centertop';
            XYZstagesPanel.Title = 'XYZ stages';
            XYZstagesPanel.Position = [11 28 283 643];
            
            % Create UpButton
            UpButton = uibutton(XYZstagesPanel, 'push','ButtonPushedFcn', @UpButtonPushed);
            UpButton.Position = [81 408 128 49];
            UpButton.Text = 'Up';
            % Create DownButton
            DownButton = uibutton(XYZstagesPanel, 'push','ButtonPushedFcn', @DownButtonPushed);
            DownButton.Position = [81 220 128 49];
            DownButton.Text = 'Down';           
            % Create LeftButton
            LeftButton = uibutton(XYZstagesPanel, 'push','ButtonPushedFcn', @LeftButtonPushed);
            LeftButton.Position = [14 271 48 134];
            LeftButton.Text = 'Left';            
            % Create RightButton
            RightButton = uibutton(XYZstagesPanel, 'push','ButtonPushedFcn', @RightButtonPushed);
            RightButton.Position = [226 271 48 134];
            RightButton.Text = 'Right';
            % Create OutButton
            OutButton = uibutton(XYZstagesPanel, 'push','ButtonPushedFcn', @OutButtonPushed);
            OutButton.Position = [70 309 71 57];
            OutButton.Text = 'Out';
            % Create InButton
            InButton = uibutton(XYZstagesPanel, 'push', 'ButtonPushedFcn', @InButtonPushed);
            InButton.Position = [148 309 71 57];
            InButton.Text = 'In';
            
            % Create StepSizePanel
            StepSizePanel = uipanel(XYZstagesPanel);
            StepSizePanel.TitlePosition = 'centertop';
            StepSizePanel.Title = 'Step Size';
            StepSizePanel.Position = [23 491 243 125];         
            % Create StepSizeEditField
            StepSizeEditField = uieditfield(StepSizePanel, 'numeric', 'ValueChangedFcn', @StepSizeEditFieldValueChanged);
            StepSizeEditField.Position = [69 57 64 44];
            StepSizeEditField.Value = hdata.misc.step_size*1000;           
            % Create umLabel
            umLabel = uilabel(StepSizePanel);
            umLabel.Position = [146 49 57 60];
            umLabel.Text = 'um';    
            % Create m1Button
            m1Button = uibutton(StepSizePanel, 'push','ButtonPushedFcn',@Buttonm1Pushed);
            m1Button.Position = [13 9 40 40];
            m1Button.Text = '-1';            
            % Create p1Button
            p1Button = uibutton(StepSizePanel, 'push','ButtonPushedFcn',@Buttonp1Pushed);
            p1Button.Position = [74 9 40 40];
            p1Button.Text = '+1';
            % Create d10Button
            d10Button = uibutton(StepSizePanel, 'push','ButtonPushedFcn',@Buttond10Pushed);
            d10Button.Position = [135 9 40 40];
            d10Button.Text = '/10';                        
            % Create x10Button
            x10Button = uibutton(StepSizePanel, 'push','ButtonPushedFcn',@Buttonx10Pushed);
            x10Button.Position = [196 9 40 40];
            x10Button.Text = 'x10';
                        
            % Create xEditFieldLabel
            xEditFieldLabel = uilabel(XYZstagesPanel);
            xEditFieldLabel.HorizontalAlignment = 'right';
            xEditFieldLabel.Position = [23 174 13 22];
            xEditFieldLabel.Text = 'x';
            % Create xEditField
            xEditField = uieditfield(XYZstagesPanel, 'numeric');
            xEditField.Position = [44 169 53 35];
            xEditField.Value = C885.getPosition_x();
            xEditField.ValueDisplayFormat = '%.4f';
            % Create yEditFieldLabel
            yEditFieldLabel = uilabel(XYZstagesPanel);
            yEditFieldLabel.HorizontalAlignment = 'right';
            yEditFieldLabel.Position = [108 174 15 22];
            yEditFieldLabel.Text = 'y';
            % Create yEditField
            yEditField = uieditfield(XYZstagesPanel, 'numeric');
            yEditField.Position = [128 169 53 35];
            yEditField.Value = C885.getPosition_y();
            yEditField.ValueDisplayFormat = '%.3f';
            % Create zEditFieldLabel
            zEditFieldLabel = uilabel(XYZstagesPanel);
            zEditFieldLabel.HorizontalAlignment = 'right';
            zEditFieldLabel.Position = [194 172 15 22];
            zEditFieldLabel.Text = 'z';
            % Create zEditField
            zEditField = uieditfield(XYZstagesPanel, 'numeric');
            zEditField.Position = [216 168 53 35];
            zEditField.Value = AtCube.getPosition_z();
            zEditField.ValueDisplayFormat = '%.3f';
            
            % Create TabGroup2
            TabGroup2 = uitabgroup(XYZstagesPanel);
            TabGroup2.Position = [28 9 228 144];
            % Create SurfaceTab
            SurfaceTab = uitab(TabGroup2);
            SurfaceTab.Title = 'Surface';
            % Create MovetotargetButton
            MovetoTargetButton = uibutton(SurfaceTab, 'push','ButtonPushedFcn',@MovetoTargetButtonPushed);
            MovetoTargetButton.Position = [69 70 97 37];
            MovetoTargetButton.Text = 'Move';
            % Create SaveTargetpositionButton
            SaveTargetposButton = uibutton(SurfaceTab, 'push','ButtonPushedFcn',@SaveTargetposButtonPushed);
            SaveTargetposButton.Position = [69 16 97 37];
            SaveTargetposButton.Text = 'Save position';
            % Create FiberTab
            FiberTab = uitab(TabGroup2);
            FiberTab.Title = 'Fiber';
            % Create MovetofiberButton
            MovetofiberButton = uibutton(FiberTab, 'push','ButtonPushedFcn',@MovetoFiberButtonPushed);
            MovetofiberButton.Position = [69 70 97 37];
            MovetofiberButton.Text = 'Move';
            % Create SaveFiberposButton
            SaveFiberposButton = uibutton(FiberTab, 'push','ButtonPushedFcn',@SaveFiberposButtonPushed);
            SaveFiberposButton.Position = [69 16 97 37];
            SaveFiberposButton.Text = 'Save position';        
                       
            % Create ShootingPanel
            ShootingPanel = uipanel(fig);
            ShootingPanel.TitlePosition = 'centertop';
            ShootingPanel.Title = 'Shooting';
            ShootingPanel.Position = [1215 16 306 647];
                       
            % Create MovetoBeamButton
            MovetoBeamButton = uibutton(ShootingPanel, 'state','ValueChangedFcn',@ShootPositioncallback);
            MovetoBeamButton.Position = [29 557 143 54];
            MovetoBeamButton.Text = 'Move to Beam';
            MovetoBeamButton.Value = 0;
       
            % Create z_axis Offset SpinnerLabel
            Spin_offset_Label = uilabel(ShootingPanel);
            Spin_offset_Label.HorizontalAlignment = 'right';
            Spin_offset_Label.Position = [59 514 35 27];
            Spin_offset_Label.Text = {'Z-axis'; 'offset (mm)'};;
            % Create z_axis Offset Spinner
            Spin_offset = uispinner(ShootingPanel,'ValueChangedFcn', @Spin_offsetValueChanged);
            Spin_offset.Step = 0.1;
            Spin_offset.Limits = [-4 4];
            Spin_offset.Position = [100 510 58 35];
            Spin_offset.Value = hdata.pow.offset;
           
            % Create Pulse SpinnerLabel
            Spin_pulse_Label = uilabel(ShootingPanel);
            Spin_pulse_Label.HorizontalAlignment = 'right';
            Spin_pulse_Label.Position = [34 463 60 42];
            Spin_pulse_Label.Text = {'Pulse '; 'Width (ms)'};
            % Create Pulse Spinner
            Spin_pulse = uispinner(ShootingPanel,'ValueChangedFcn', @Spin_pulseValueChanged);
            Spin_pulse.Step = 1;
            Spin_pulse.Limits = [0 4000];
            Spin_pulse.Position = [100 467 58 35];
            Spin_pulse.Value = hdata.pulse.laser_width;

            % Create NoshotsSpinnerLabel
            NoshotsSpinnerLabel = uilabel(ShootingPanel);
            NoshotsSpinnerLabel.HorizontalAlignment = 'right';
            NoshotsSpinnerLabel.Position = [40 434 56 22];
            NoshotsSpinnerLabel.Text = 'No. shots';
            % Create NoshotsSpinner
            NoshotsSpinner = uispinner(ShootingPanel, 'ValueChangedFcn', @NoshotsSpinnerValueChanged);
            NoshotsSpinner.Limits = [1 Inf];
            NoshotsSpinner.Position = [100 426 58 35];
            NoshotsSpinner.Value = hdata.misc.num_shots;

            % Create MillingPanel
            MillingPanel = uipanel(ShootingPanel);
            MillingPanel.TitlePosition = 'centertop';
            MillingPanel.Title = 'Milling';
            MillingPanel.Position = [12 107 285 306];      
            % Create NoringsSpinnerLabel
            NoringsSpinnerLabel = uilabel(MillingPanel);
            NoringsSpinnerLabel.HorizontalAlignment = 'right';
            NoringsSpinnerLabel.Position = [14 218 53 22];
            NoringsSpinnerLabel.Text = 'No. rings';
            % Create NoringsSpinner
            NoringsSpinner = uispinner(MillingPanel, 'ValueChangedFcn', @NoringsSpinnerValueChanged);
            NoringsSpinner.Position = [72 212 58 35];
            NoringsSpinner.Step = 1;
            NoringsSpinner.Value = hdata.mill.n_rings;
            % Create FreqSpinner_2Label
            FreqSpinnerLabel = uilabel(MillingPanel);
            FreqSpinnerLabel.HorizontalAlignment = 'right';
            FreqSpinnerLabel.Position = [37 159 31 22];
            FreqSpinnerLabel.Text = 'Freq.';
            % Create FreqSpinner_2
            FreqSpinner = uispinner(MillingPanel, 'ValueChangedFcn', @FreqSpinnerValueChanged);
            FreqSpinner.Position = [72 153 58 35];
            FreqSpinner.Step = 3;
            FreqSpinner.Value = hdata.mill.n_shots;
            % Create RotationSpinnerLabel
            RotationSpinnerLabel = uilabel(MillingPanel);
            RotationSpinnerLabel.HorizontalAlignment = 'right';
            RotationSpinnerLabel.Position = [17 101 50 22];
            RotationSpinnerLabel.Text = 'Rotation';
            % Create RotationSpinner
            RotationSpinner = uispinner(MillingPanel, 'ValueChangedFcn', @RotSpinnerValueChanged);
            RotationSpinner.Position = [72 95 58 35];
            RotationSpinner.Step = 0.05;
            RotationSpinner.Value = hdata.mill.rot;
            % Create RadiusSpinnerLabel
            RadiusSpinnerLabel = uilabel(MillingPanel);
            RadiusSpinnerLabel.HorizontalAlignment = 'right';
            RadiusSpinnerLabel.Position = [22 40 43 22];
            RadiusSpinnerLabel.Text = 'Radius';
            % Create RadiusSpinner
            RadiusSpinner = uispinner(MillingPanel, 'ValueChangedFcn', @RadSpinnerValueChanged);
            RadiusSpinner.Position = [69 33 61 37];
            RadiusSpinner.Step = 1;
            RadiusSpinner.Value = hdata.mill.rad;
            % Create DirectionSwitchLabel
            DirectionSwitchLabel = uilabel(MillingPanel);
            DirectionSwitchLabel.HorizontalAlignment = 'center';
            DirectionSwitchLabel.Position = [182 232 53 22];
            DirectionSwitchLabel.Text = 'Direction';
            % Create DirectionSwitch
            DirectionSwitch = uiswitch(MillingPanel, 'slider','ValueChangedFcn',@DirectionSwitchValueChanged);
            DirectionSwitch.Items = {'Winter', 'Summer'};
            DirectionSwitch.Position = [185 254 45 20];
            DirectionSwitch.Value = 'Winter';
            % Create StyleSwitchLabel
            StyleSwitchLabel = uilabel(MillingPanel);
            StyleSwitchLabel.HorizontalAlignment = 'center';
            StyleSwitchLabel.Position = [181 167 53 22];
            StyleSwitchLabel.Text = 'Style';
            % Create StyleSwitch
            StyleSwitch = uiswitch(MillingPanel, 'slider','ValueChangedFcn',@StyleSwitchValueChanged);
            StyleSwitch.Items = {'Radial', 'C.centric'};
            StyleSwitch.Position = [183 194 45 20];
            StyleSwitch.Value = 'Radial';
            % Create RandomSwitchLabel
            RandomSwitchLabel = uilabel(MillingPanel);
            RandomSwitchLabel.HorizontalAlignment = 'center';
            RandomSwitchLabel.Position = [181 111 53 22];
            RandomSwitchLabel.Text = 'Non-consecutive';
            % Create RandomSwitch
            RandomSwitch = uiswitch(MillingPanel, 'slider','ValueChangedFcn',@RandomSwitchValueChanged);
            RandomSwitch.Items = {'Off', 'On'};
            RandomSwitch.Position = [186 134 45 20];
            RandomSwitch.Value = 'Off';

            % Create EccentricityPanel
            EccentricityPanel = uipanel(MillingPanel);
            EccentricityPanel.TitlePosition = 'centertop';
            EccentricityPanel.Title = 'Eccentricity';
            EccentricityPanel.Position = [139 8 140 104];
            % Create AngleSpinnerLabel
            AngleSpinnerLabel = uilabel(EccentricityPanel);
            AngleSpinnerLabel.HorizontalAlignment = 'right';
            AngleSpinnerLabel.Position = [28 11 36 22];
            AngleSpinnerLabel.Text = 'Angle';
            % Create AngleSpinner
            AngleSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @AngleSpinnerValueChanged);
            AngleSpinner.Limits = [0 270];
            AngleSpinner.Position = [75 5 58 35];
            AngleSpinner.Value = hdata.mill.ecc_angle;
            % Create RadiiratioSpinnerLabel
            RadiiratioSpinnerLabel = uilabel(EccentricityPanel);
            RadiiratioSpinnerLabel.HorizontalAlignment = 'right';
            RadiiratioSpinnerLabel.Position = [2 53 60 22];
            RadiiratioSpinnerLabel.Text = 'Radii ratio';
            % Create RadiiratioSpinner
            RadiiratioSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @RadiiratioSpinnerValueChanged);
            RadiiratioSpinner.Step = 0.05;
            RadiiratioSpinner.Limits = [0 2];
            RadiiratioSpinner.Position = [75 47 58 35];
            RadiiratioSpinner.Value = hdata.mill.ecc_degree;
            
            % Create BeamPositionPanel
            BeamPositionPanel = uipanel(ShootingPanel);
            BeamPositionPanel.TitlePosition = 'centertop';
            BeamPositionPanel.Title = 'Beam Position';
            BeamPositionPanel.Position = [197 426 100 191];
            % Create beam_xEditFieldLabel
            beam_xEditFieldLabel = uilabel(BeamPositionPanel);
            beam_xEditFieldLabel.HorizontalAlignment = 'right';
            beam_xEditFieldLabel.Position = [16 122 10 22];
            beam_xEditFieldLabel.Text = 'x';
            % Create beam_xEditField
            beam_xEditField = uieditfield(BeamPositionPanel, 'numeric','ValueChangedFcn', @beam_xValueChanged);
            beam_xEditField.Position = [31 116 60 35];
            beam_xEditField.Value = hdata.beam.x_pos;
            beam_xEditField.ValueDisplayFormat = '%.4f';
            % Create beam_yEditFieldLabel
            beam_yEditFieldLabel = uilabel(BeamPositionPanel);
            beam_yEditFieldLabel.HorizontalAlignment = 'right';
            beam_yEditFieldLabel.Position = [16 74 10 22];
            beam_yEditFieldLabel.Text = 'y';
            % Create beam_yEditField
            beam_yEditField = uieditfield(BeamPositionPanel, 'numeric','ValueChangedFcn', @beam_yValueChanged);
            beam_yEditField.Position = [31 66 60 35];
            beam_yEditField.Value = hdata.beam.y_pos;
            beam_yEditField.ValueDisplayFormat = '%.4f';
            % Create beam_zEditFieldLabel
            beam_zEditFieldLabel = uilabel(BeamPositionPanel);
            beam_zEditFieldLabel.HorizontalAlignment = 'right';
            beam_zEditFieldLabel.Position = [16 23 10 22];
            beam_zEditFieldLabel.Text = 'z';
            % Create beam_zEditField
            beam_zEditField = uieditfield(BeamPositionPanel, 'numeric', 'ValueChangedFcn', @beam_zValueChanged);
            beam_zEditField.Position = [31 16 60 35];
            beam_zEditField.Value = hdata.beam.z_pos;
            beam_zEditField.ValueDisplayFormat = '%.4f';
            
            % Create BSLampLabel
            BSLampLabel = uilabel(ShootingPanel);
            BSLampLabel.HorizontalAlignment = 'right';
            BSLampLabel.Position = [12 21 19 22];
            BSLampLabel.Text = 'BS';
            % Create BSLamp
            BSLamp = uilamp(ShootingPanel);
            BSLamp.Position = [38 19 26 26];
            if C885.getPosition_z() == 0
                    BSLamp.Color = 'green';
                else
                    BSLamp.Color = 'red';
            end

            % Create ShootPosLampLabel
            ShootPosLabel = uilabel(ShootingPanel);
            ShootPosLabel.HorizontalAlignment = 'right';
            ShootPosLabel.Position = [9 55 29 26];
            ShootPosLabel.Text = 'Fiber';
            % Create FiberLamp
            ShootPosLamp = uilamp(ShootingPanel);
            ShootPosLamp.Position = [40 57 24 24];
            if AtCube.getPosition_z() < 1+hdata.beam.z_pos
                    ShootPosLamp.Color = 'green';
                else
                    ShootPosLamp.Color = 'red';
            end
            
            % Create ShootButton
            ShootButton = uibutton(ShootingPanel, 'push', 'ButtonPushedFcn', @Shootcallback);
            ShootButton.Text = 'Shoot';
            ShootButton.Position = [93 14 155 45];
            
            % Create GridCheckBox
            GridCheckBox = uicheckbox(ShootingPanel, 'ValueChangedFcn', @GridCheckBoxValueChanged);
            GridCheckBox.Text = 'Grid';
            GridCheckBox.Position = [89 70 45 22];
            GridCheckBox.Value = hdata.misc.grid;
            
%             % Create ShootingSwitchLabel
%             ShootingSwitchLabel = uilabel(ShootingPanel);
%             ShootingSwitchLabel.HorizontalAlignment = 'center';
%             ShootingSwitchLabel.Position = [193 71 45 20];
%             ShootingSwitchLabel.Text = 'Shooting';
            % Create ShootingSwitch
            ShootingSwitch = uiswitch(ShootingPanel, 'slider', 'ValueChangedFcn', @ShootingSwitchValueChanged);
            ShootingSwitch.Items = {'Single', 'Milling'};
            ShootingSwitch.Position = [193 71 45 20];
            ShootingSwitch.Value = 'Single';
            
            % Create Shoot Lamp
            ShootLamp = uilamp(ShootingPanel);
            ShootLamp.Position = [263 29 21 21];
            ShootLamp.Color = 'red';
end

