function simpleApp2dd(varargin)
%% Startup Processes

uEye_camera(0); % initializing camera
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
fig.Tag = 'Fiber Machining';
fig.HandleVisibility = 'on'

%% Begin assembling data structure for configuration
hdata = guidata(fig);
setting_file = 'shoot_settings.mat';
load(setting_file, 'hdata');
hdata.misc.z_start_pos = AtCube.getPosition_z();
uEye_camera(8, hdata.misc.exp); %exposure time

hdata.pulse.trig_pulse_delay = 10;
hdata.pulse.shoot_trig_delay = 0;
hdata.pulse.end_delay = 1;
hdata.pulse.laser_shutter_delay1 = 20;
hdata.pulse.laser_shutter_delay2 = 600;

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
            cleave_tension = CleaveTensionSpinner.Value;

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
                imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir2, folder, file),'.tif')))
                saveas(pltsurf, char(strcat(fullfile(dir, folder, file),'surf.fig')))
                saveas(pltsurf, char(strcat(fullfile(dir2, folder, file),'surf.fig')))
                Surface = hdata.surf;
                save(char(strcat(fullfile(dir, folder, file),'.mat')),'Surface')
                save(char(strcat(fullfile(dir2, folder, file),'.mat')),'Surface')
                saveas(pltcross, char(strcat(fullfile(dir, folder, file),'cross.fig')))
                saveas(pltcross, char(strcat(fullfile(dir2, folder, file),'cross.fig')))
                if hdata.misc.flat == 1
                    imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir2, 'Cleaving Face',num2str(cleave_tension), file),'.tif')))
                    saveas(pltsurf, char(strcat(fullfile(dir2, 'Cleaving Face',num2str(cleave_tension), file),'surf.fig')))
                    save(char(strcat(fullfile(dir2, 'Cleaving Face',num2str(cleave_tension), file),'.mat')),'Surface')
                    saveas(pltcross, char(strcat(fullfile(dir2, 'Cleaving Face',num2str(cleave_tension), file),'cross.fig')))
                end
            else
                imwrite(uint8(hdata.misc.Image(298:798,718:1218)), char(strcat(fullfile(dir, folder, file),'pre.tif')))
            end
        end

        % Value changing function: CleaveTensionSpinner
        function CleaveTensionSpinnerValueChanged(app, event)
            value = CleaveTensionSpinner.Value;
            hdata = guidata(fig);
            hdata.misc.cleave_tension = value;
            guidata(fig,hdata)
            
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
    
        % Value changing function: ExposuremsSpinner
        function ExposuremsSpinnerValueChanged(app, event)
            value = ExposuremsSpinner.Value;
            hdata = guidata(fig);
            hdata.misc.exp = value;
            uEye_camera(8, hdata.misc.exp); %exposure time
            guidata(fig,hdata)
            
        end
    
        % Button pushed function:  Max Contrast
        function Find_max_contrast(src,event)
            
            hdata = guidata(fig);
            sample_c = [548,968]; % sample centre, pixels
            sample_a = 10; % sample width/height, pixels
            hdata.misc.lambd = 570*10^(-6); %wavelegth*2?
            %hdata.misc.lambd = 283*10^(-6); %
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
            
            
            
%             N = normalize(contrast);
%             Nt = N.';
%             z_listt = z_list.';
%             interference = cat(2, z_listt, Nt);
%             dir = 'C:\Users\co2la\Documents\Data';
%             time = clock;
%             year = sprintf('%.0f',time(1,1));
%             month = sprintfc('%02d', time(1,2));
%             day = sprintfc('%02d', time(1,3));
%             hour = sprintf('%.0f',time(1,4));
%             min = sprintfc('%02d', time(1,5));
%             folder = append(year, month, day);
%             file = append(hour, min);
                 
            
            %writematrix(interference,char(strcat(fullfile(dir, folder, file),'normalize.csv')));
            %export(normalize(contrast),'XLSFile','normalize.xlsx')
            
            
            
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

        % Button pushed function: Reconstruct
        function Reconstruct_FiberButtonPushed(src,event)
            stop(tm)
            hdata = guidata(fig);
            if hdata.misc.max_contrast_status == 0
                hdata.misc.z_start_pos = AtCube.getPosition_z();
            end
            
                if TabGroup2.SelectedTab == SurfaceTab
                    Surface_reconstruction(hdata.misc.lambd, hdata.misc.lam_step, hdata.misc.z_start_pos, AtCube, RCheckBox, GCheckBox, BCheckBox);
                else
                    [slice_x, slice_y, I1, hdata.surf] = Fiber_reconstruction(hdata.misc.lambd, hdata.misc.lam_step, hdata.misc.z_start_pos, AtCube, RCheckBox, GCheckBox, BCheckBox);
                    close(findobj('type','figure','name','Cross Section'))
            try
                    [r_curv, hdata.misc.flat] = curv_rad(slice_x);
                    [r_curv, hdata.misc.flat] = curv_rad(slice_y);
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
        function GridSwitchValueChanged(src, ~)
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
                     hdata.mill.random = 'Off';
                     guidata(fig,hdata)
                     dotmillingschematic();
                     plot_milling();
                case 'On'
                     hdata.mill.random = 'On';
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
                C885.move_z(24.13);
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
            C885.move_z(24.13);
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
            hdata = guidata(fig);
            value = Spin_pulse.Value+hdata.pulse.laser_shutter_delay1+hdata.pulse.laser_shutter_delay2;
            hdata.pulse.laser_low_time1 = value/2;
            hdata.pulse.laser_low_time2 = value/2;
            % set pulse width
            mex_ok_interface('swi', 3, hdata.pulse.laser_low_time1*10);
            mex_ok_interface('swi', 7, hdata.pulse.laser_low_time2*10);
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
    
        % Value changed function: ShotDelaySpinner
        function ShotDelaySpinnerValueChanged(app, event)
            value = ShotDelaySpinner.Value;
            hdata = guidata(fig);
            hdata.mill.shot_delay = value;
            guidata(fig,hdata);            
        end

        % Button pushed function: MillingPatternButton
        function MillingPatternButtonPushed(app, event)
            hdata = guidata(fig);
            [hdata.mill.x_dot, hdata.mill.y_dot, hdata.mill.totalshots] = Mill_pattern_GUI(fig)
            guidata(fig,hdata);
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
%                      hdata.misc.shoot = 0;
%                      guidata(fig,hdata)
                     delete(hdata.mill_plot);
                case 'Milling'
%                      hdata.misc.shoot = 1;
%                      guidata(fig,hdata)
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
                        switch ShootingSwitch.Value
                            case 'Single'
                                for i=1:hdata.misc.num_shots
                                    tic;
%                                     mex_ok_interface('ati', 64, 1);% activate trigger
                                    mex_ok_interface('ati', 64, 2);% activate trigger
                                    pause((0.002*hdata.pulse.laser_low_time2)+hdata.mill.shot_delay)
                                    toc
                                    
                                end
                            case 'Milling'
                                z_0 = AtCube.getPosition_z();
                                x_shoot = C885.getPosition_x();
                                y_shoot = C885.getPosition_y();
                                %AtCube.move_z(z_0+hdata.pow.offset);
                                for i=1:hdata.misc.num_shots
                                    for s = 1:hdata.mill.totalshots
                                        C885.move_x(x_shoot+(hdata.mill.x_dot(s)*0.000329));
                                        C885.move_y(y_shoot+(hdata.mill.y_dot(s)*0.000329));
%                                         mex_ok_interface('ati', 64, 1);% activate trigger
                                        mex_ok_interface('ati', 64, 2);% activate trigger
                                        pause((0.002*hdata.pulse.laser_low_time2)+hdata.mill.shot_delay)
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
                                    C885.move_x(XPOS+ix*0.1);
                                    mex_ok_interface('ati', 64, 1);% activate trigger
                                    mex_ok_interface('ati', 64, 2);% activate trigger
                                    pause((0.002*hdata.pulse.laser_low_time2)+hdata.mill.shot_delay)
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
    Im_mono = 0.55*RVal*I(:,:,2)+0.35*GVal*I(:,:,3)+0.1*BVal*I(:,:,1); % Create monochrome image from the 3 colour channels
    ImBG_autoCont = rescale(Im_mono,0,256); % Rescale pixel values to maximise contrast, keep floating point precision
end

function updateImage(src, evt)
    hdata = guidata(fig);
    Im_opt_cont = grabImage(RCheckBox, GCheckBox, BCheckBox); % Grab processed image from camera
    hdata.misc.Image = uint8(Im_opt_cont); % Convert pixel values to integers for use as an image. Save image to GUI data to be used for surface reconstruction etc.
    imshow(hdata.misc.Image, 'Parent', Imageax); % Update the GUI axes with new image
    %set(Imageax, 'CData', Im_opt_cont);
    guidata(fig,hdata);
end

function plot_milling()
milldata = guidata(fig);
f = linspace(1,10,hdata.mill.totalshots);
xc=968;
yc=548;
hdata.mill_plot = scatter(gridax, (xc+hdata.mill.x_dot(:))',(yc+hdata.mill.y_dot(:))',[],'+','LineWidth',2,'CData',f);
drawcircle(gridax,'Center',[968,548],'Radius',200,'Color','k','FaceAlpha',0.01,'LineWidth',0.1,'MarkerSize',0.1);
guidata(fig,hdata)
end
hdata.myFun1=@plot_milling;
setappdata(fig,'fun_handles',hdata);
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
            
            % Create CleaveTensionSpinnerLabel
            CleaveTensionSpinnerLabel = uilabel(fig);
            CleaveTensionSpinnerLabel.HorizontalAlignment = 'right';
            CleaveTensionSpinnerLabel.Position = [325 608 38 28];
            CleaveTensionSpinnerLabel.Text = {'Cleave'; 'Tension'};
            % Create CleaveTensionSpinner
            CleaveTensionSpinner = uispinner(fig, 'ValueChangedFcn',@CleaveTensionSpinnerValueChanged);
            CleaveTensionSpinner.Step = 5;
            CleaveTensionSpinner.Limits = [100 300];
            CleaveTensionSpinner.Position = [370 608 55 28];
            CleaveTensionSpinner.Value = hdata.misc.cleave_tension;
            
            % Create SaveButton
            SaveButton = uibutton(fig,'push','ButtonPushedFcn',@SaveImagecallback);
            SaveButton.Text = 'Save';
            SaveButton.Position = [437 603 125 60];
            
            % Create AutoCentreButton
            AutoCentreButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoCentreButtonPushed);
            AutoCentreButton.FontWeight = 'bold';
            AutoCentreButton.FontColor = [0.0745 0.6235 1];
            AutoCentreButton.Position = [573 603 140 60];
            AutoCentreButton.Text = {'Auto-Centre'; ''};
            
%             % Create AutoFocusButton
%             AutoFocusButton = uibutton(fig, 'push','ButtonPushedFcn', @AutoFocusButtonPushed);
%             AutoFocusButton.FontWeight = 'bold';
%             AutoFocusButton.FontColor = [0.0745 0.6235 1];
%             AutoFocusButton.Position = [794 603 180 60];
%             AutoFocusButton.Text = {'Auto-Focus'; ''};
            
            % Create ExposuremsSpinnerLabel
            ExposuremsSpinnerLabel = uilabel(fig);
            ExposuremsSpinnerLabel.HorizontalAlignment = 'right';
            ExposuremsSpinnerLabel.Position = [767 608 56 27];
            ExposuremsSpinnerLabel.Text = {'Exposure'; '(ms)'};
            % Create ExposuremsSpinner
            ExposuremsSpinner = uispinner(fig, 'ValueChangedFcn',@ExposuremsSpinnerValueChanged);
            ExposuremsSpinner.Position = [830 607 62 30];
            ExposuremsSpinner.Step = 0.2;
            ExposuremsSpinner.Value = hdata.misc.exp;
            
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
            ReconstructFiberButton.Text = 'Reconstruct';
            
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
            Spin_pulse.Value = (2*hdata.pulse.laser_low_time1)-(hdata.pulse.laser_shutter_delay1+hdata.pulse.laser_shutter_delay2);
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
            % Create ShotDelaySpinnerLabel
            ShotDelaySpinnerLabel = uilabel(ShootingPanel);
            ShotDelaySpinnerLabel.HorizontalAlignment = 'right';
            ShotDelaySpinnerLabel.Position = [32 388 64 22];
            ShotDelaySpinnerLabel.Text = 'Shot Delay';
            % Create ShotDelaySpinner
            ShotDelaySpinner = uispinner(ShootingPanel, 'ValueChangedFcn', @ShotDelaySpinnerValueChanged);
            ShotDelaySpinner.Position = [100 380 58 35];
            ShotDelaySpinner.Value = hdata.mill.shot_delay;
            ShotDelaySpinner.Step = 0.1;
            
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
            
            % Create KEPositionPanel
            KEPositionPanel = uipanel(ShootingPanel);
            KEPositionPanel.TitlePosition = 'centertop';
            KEPositionPanel.Title = 'KE Position';
            KEPositionPanel.Position = [197 230 100 183];
            % Create beam_xEditFieldLabel
            knife_xEditFieldLabel = uilabel(KEPositionPanel);
            knife_xEditFieldLabel.HorizontalAlignment = 'right';
            knife_xEditFieldLabel.Position = [16 122 10 22];
            knife_xEditFieldLabel.Text = 'x';
            % Create knife_xEditField
            knife_xEditField = uieditfield(KEPositionPanel, 'numeric','ValueChangedFcn', @knife_xValueChanged);
            knife_xEditField.Position = [31 116 60 35];
            knife_xEditField.Value = hdata.knife.x_pos;
            knife_xEditField.ValueDisplayFormat = '%.4f';
            % Create knife_yEditFieldLabel
            knife_yEditFieldLabel = uilabel(KEPositionPanel);
            knife_yEditFieldLabel.HorizontalAlignment = 'right';
            knife_yEditFieldLabel.Position = [16 74 10 22];
            knife_yEditFieldLabel.Text = 'y';
            % Create knife_yEditField
            knife_yEditField = uieditfield(KEPositionPanel, 'numeric','ValueChangedFcn', @knife_yValueChanged);
            knife_yEditField.Position = [31 66 60 35];
            knife_yEditField.Value = hdata.knife.y_pos;
            knife_yEditField.ValueDisplayFormat = '%.4f';
            % Create knife_zEditFieldLabel
            knife_zEditFieldLabel = uilabel(KEPositionPanel);
            knife_zEditFieldLabel.HorizontalAlignment = 'right';
            knife_zEditFieldLabel.Position = [16 23 10 22];
            knife_zEditFieldLabel.Text = 'z';
            % Create knife_zEditField
            knife_zEditField = uieditfield(KEPositionPanel, 'numeric', 'ValueChangedFcn', @knife_zValueChanged);
            knife_zEditField.Position = [31 16 60 35];
            knife_zEditField.Value = hdata.knife.z_pos;
            knife_zEditField.ValueDisplayFormat = '%.4f';
            
            % Create MillingPatternButton
            MillingPatternButton = uibutton(ShootingPanel, 'push','ButtonPushedFcn', @MillingPatternButtonPushed);;
            MillingPatternButton.Position = [26 131 254 69];
            MillingPatternButton.Text = 'Milling Pattern';
            
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

