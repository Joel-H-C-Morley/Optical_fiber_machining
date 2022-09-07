function [shot_pos_x, shot_pos_y, n_pos] = Mill_pattern_GUI(varargin)
    
            % Create UIFigure and hide until all components are created
            UIFigure = uifigure('Visible', 'off');
            UIFigure.Position = [100 100 680 350];
            UIFigure.Name = 'MATLAB UIFigure';
            
    milldata = guidata(UIFigure);
    
    setting_file = 'mill_settings.mat';
    load(setting_file, 'milldata');

      milldata.line = 'Points';
%     milldata.freq = [4 5 6 4 6 6];
%     milldata.rad = [10 20 30 40 60 80];
%     milldata.rot = [0.05 0.1 0.15 0.2 0.25 0.3];
%     milldata.n_rings = 4;
%     milldata.ratio = 1;
%     milldata.ecc_angle = 0;
%     milldata.x2 = ones(1,6);
%     milldata.x3 = ones(1,6);
%     milldata.x4 = ones(1,6);
%     milldata.x5 = ones(1,6);
%     milldata.x6 = ones(1,6);
%     milldata.x1 = ones(1,6);
%     milldata.y2 = ones(1,6);
%     milldata.y3 = ones(1,6);
%     milldata.y4 = ones(1,6);
%     milldata.y5 = ones(1,6);
%     milldata.y6 = ones(1,6);
%     milldata.y1 = ones(1,6);
    shot_pos_x = milldata.x_dot;
    shot_pos_y = milldata.y_dot;
    n_pos = length(milldata.x_dot);
    guidata(UIFigure,milldata);

        % Button pushed function: SaveButton
        function SaveButtonPushed(UIFigure, event)
            save(setting_file, 'milldata');
            n_pos = sum(milldata.freq);
            h = findobj('Tag','Fiber Machining');
            hdata = guidata(h);
            hdata.mill.x_dot = 3.29*milldata.x_dot;
            hdata.mill.y_dot = 3.29*milldata.y_dot;
            hdata.mill.totalshots = length(milldata.x_dot);
            func=getappdata(h,'fun_handles');
            guidata(h,hdata);
            fprintf('Frequency is: [');
            fprintf('%g, ', milldata.freq(1:milldata.n_rings-1));
            fprintf('%g]\n', milldata.freq(milldata.n_rings));
            fprintf('Radius is: [');
            fprintf('%g, ', milldata.rad(1:milldata.n_rings-1));
            fprintf('%g]\n', milldata.rad(milldata.n_rings));
            fprintf('Frequency is: [');
            fprintf('%g, ', milldata.rot(1:milldata.n_rings-1));
            fprintf('%g]\n', milldata.rot(milldata.n_rings));
        end
    
        % Value changed function: RingDropDown
        function RingDropDownValueChanged(UIFigure, event)
            value = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch value
                case {1, 2, 3, 4, 5, 6}
                   FrequencySpinner.Value = milldata.freq(value);
                   RadiusSpinner.Value = milldata.rad(value);
                   RotationSpinner.Value = milldata.rot(value);
                   milldata.ring = value;
                   guidata(UIFigure,milldata);   
            end
        end

        % Value changed function: NoringsSpinner
        function NoringsSpinnerValueChanged(UIFigure, event)
            value = NoringsSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.n_rings = value;
            guidata(UIFigure,milldata)
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: FrequencySpinner
        function FrequencySpinnerValueChanged(UIFigure, event)
            value = FrequencySpinner.Value;
            ring_num = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch ring_num
                case {1, 2, 3, 4, 5, 6}
                    milldata.freq(ring_num) = value;
                    guidata(UIFigure,milldata);
            end
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling(); 
        end

        % Value changed function: RadiusSpinner
        function RadiusSpinnerValueChanged(UIFigure, event)
            value = RadiusSpinner.Value;
            ring_num = RingDropDown.Value;
            milldata = guidata(UIFigure);
            switch ring_num
                case {1, 2, 3, 4, 5, 6}
                    milldata.rad(ring_num) = value;
                    guidata(UIFigure,milldata);
            end
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: RotationSpinner
        function RotationSpinnerValueChanged(UIFigure, event)
            value = RotationSpinner.Value;
            milldata = guidata(UIFigure);
            ring_num = RingDropDown.Value;
            switch ring_num
                case {1, 2, 3, 4, 5, 6}
                    milldata.rot(ring_num) = value;
                    guidata(UIFigure,milldata);
            end
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: Consecutive
        function ConsecutiveCheckBoxValueChanged(UIFigure, event)
            value = ConsecutiveCheckBox.Value;
            milldata = guidata(UIFigure);
            milldata.consecutive = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();            
        end

        % Value changed function: clockwise
        function ClockwiseCheckBoxValueChanged(UIFigure, event)
            value = ClockwiseCheckBox.Value;
            milldata = guidata(UIFigure);
            milldata.clockwise = value;
            guidata(UIFigure,milldata);
            [millx, milly] = dotmillingschematic(RingDropDown.Value);
            plot_milling();             
        end
 
        % Value changed function: LineSwitch
        function LineSwitchValueChanged(UIFigure, event)
            value = LineSwitch.Value;
        plot_milling();            
        end

        % Value changed function: RatioSpinner
        function RatioSpinnerValueChanged(UIFigure, event)
            value = RatioSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.ratio = value;
            guidata(UIFigure,milldata);
            for idx = 1:RingDropDown.Value
                [millx, milly] = dotmillingschematic(idx);
                plot_milling();
            end             
        end

        % Value changed function: AngleSpinner
        function AngleSpinnerValueChanged(UIFigure, event)
            value = AngleSpinner.Value;
            milldata = guidata(UIFigure);
            milldata.ecc_angle = value;
            guidata(UIFigure,milldata);
            for idx = 1:RingDropDown.Value
                [millx, milly] = dotmillingschematic(idx);
                plot_milling();
            end
        end

    %% Helper functions
    
    function [xp, yp] = dotmillingschematic(ring_num)
    milldata = guidata(UIFigure);  
    for j = 1:milldata.n_rings
        
        r = milldata.rad(j); % Fiber radius in mm
        deg_of_ecc = milldata.ratio; %eccentricity ratio
        phi = milldata.ecc_angle; % angle of ellipse
        rx = r.*deg_of_ecc; %x radius given eccentricity
        ry = r;
        r = @(theta) (rx*ry*(r^2))./(sqrt(((rx*r)^2)*sin((phi*pi/180)+theta).^2+((ry*r)^2)*cos((phi*pi/180)+theta).^2)); % Distance from origin of point (radius) based on the angle of rotation about the origin (constant if no eccentricity)
        theta=linspace(0+(((pi/180)*milldata.rot(j))),...% Angle of rotation about origin required between each point on ring
                            (2*pi+((pi/180)*milldata.rot(j)))-((2*pi+((pi/180)*milldata.rot(j)))/milldata.freq(j)),...
                            milldata.freq(j));
        xp = (rx.*cos(theta).*cos(phi*pi/180))-(ry.*sin(theta).*sin(phi*pi/180)); % x co-ordinates for positions on the ring
        yp = (rx.*cos(theta).*sin(phi*pi/180))+(ry.*sin(theta).*cos(phi*pi/180)); % y co-ordinates for positions on the ring
        milldata.x_pos{j} = xp;
        milldata.y_pos{j} = yp;

        if milldata.consecutive == 1 % Make milling pattern random or maximally separated
            else
            ring_freq = length(milldata.x_pos{j});
            section =ring_freq/3;
            idx = ones(1,length(milldata.x_pos{j}));
            for i = 2:length(milldata.x_pos{j})
                if idx(i-1)+section > ring_freq
                    if i > ring_freq 
                            idx(i) = idx(i-ring_freq)+ring_freq;
                        else
                        idx(i) = i-(2*((i-1)/3));
                    end
                    else
                    idx(i) = idx(i-1)+section;
                end

            end
            milldata.x_pos{j} = milldata.x_pos{j}(idx);
            milldata.y_pos{j} = milldata.y_pos{j}(idx);
        end      
    end
    
    finalarray_x=[];
    finalarray_y=[];
    for idx = 1:milldata.n_rings % put x and y values from a table of values separating rings into an x and y list to be exported to the the GUI data
        finalarray_x = [finalarray_x, milldata.x_pos{idx}];
        finalarray_y = [finalarray_y, milldata.y_pos{idx}];
        milldata.x_dot = finalarray_x;
        milldata.y_dot = finalarray_y;
        if milldata.clockwise == 1 % Determine milling direction   
        else
            milldata.x_dot = flip(milldata.x_dot);
            milldata.y_dot = flip(milldata.y_dot);
        end
    end    
    guidata(UIFigure,milldata);   
    end

    function plot_milling()
        milldata = guidata(UIFigure);
        totalshots = length(milldata.x_dot);
        f = linspace(0,10,length( milldata.x_dot));
        xc=0;
        yc=0;
        
        %hold( UIAxes, 'on');
        switch LineSwitch.Value
            case 'Line'
            milldata.mill_plot = plot(UIAxes, milldata.x_dot, milldata.y_dot','LineWidth', 0.5, 'Marker', '+', 'MarkerSize',5, 'MarkerEdgeColor', 'red');
            case 'Points'
            milldata.mill_plot = scatter(UIAxes, milldata.x_dot, milldata.y_dot',[],'+','LineWidth',2,'CData',f);
        end
        drawcircle(UIAxes,'Center',[0,0],'Radius',60,'Color','w','FaceAlpha',0.1,'LineWidth',2,'MarkerSize',2);
        %hold( UIAxes, 'off');
        guidata(UIFigure,milldata);
    end
    
    %% GUI Layout

        % Create FrequencySpinnerLabel
        FrequencySpinnerLabel = uilabel(UIFigure);
        FrequencySpinnerLabel.HorizontalAlignment = 'right';
        FrequencySpinnerLabel.Position = [34 163 62 22];
        FrequencySpinnerLabel.Text = 'Frequency';
        % Create FrequencySpinner
        FrequencySpinner = uispinner(UIFigure, 'ValueChangedFcn',@FrequencySpinnerValueChanged);
        FrequencySpinner.Position = [108 154 76 42];
        FrequencySpinner.Limits = [1 100];
        FrequencySpinner.Step = 3;
        FrequencySpinner.Value = milldata.freq(1);

        % Create RadiusSpinnerLabel
        RadiusSpinnerLabel = uilabel(UIFigure);
        RadiusSpinnerLabel.HorizontalAlignment = 'right';
        RadiusSpinnerLabel.Position =  [25 105 71 22];
        RadiusSpinnerLabel.Text = 'Radius (um)';
        % Create RadiusSpinner
        RadiusSpinner = uispinner(UIFigure, 'ValueChangedFcn', @RadiusSpinnerValueChanged);
        RadiusSpinner.Position = [108 95 76 42];
        RadiusSpinner.Limits = [1 100];
        RadiusSpinner.Value = milldata.rad(1);

        % Create RotationSpinnerLabel
        RotationSpinnerLabel = uilabel(UIFigure);
        RotationSpinnerLabel.HorizontalAlignment = 'right';
        RotationSpinnerLabel.Position = [14 43 82 22];
        RotationSpinnerLabel.Text = 'Rotation (deg)';
        % Create RotationSpinner
        RotationSpinner = uispinner(UIFigure, 'ValueChangedFcn', @RotationSpinnerValueChanged);
        RotationSpinner.Position = [108 33 76 42];
        RotationSpinner.Limits = [0 360];
        RotationSpinner.Step = [1];
        RotationSpinner.Value = milldata.rot(1);

        % Create NoringsSpinnerLabel
        NoringsSpinnerLabel = uilabel(UIFigure);
        NoringsSpinnerLabel.HorizontalAlignment = 'right';
        NoringsSpinnerLabel.Position = [42 295 54 22];
        NoringsSpinnerLabel.Text = 'No. rings';
        % Create NoringsSpinner
        NoringsSpinner = uispinner(UIFigure, 'ValueChangedFcn', @NoringsSpinnerValueChanged);
        NoringsSpinner.Position = [108 285 76 42];
        RadiusSpinner.Limits = [1 100];
        NoringsSpinner.Value = milldata.n_rings;

        % Create ShootingOrderPanel
        ShootingOrderPanel = uipanel(UIFigure);
        ShootingOrderPanel.TitlePosition = 'centertop';
        ShootingOrderPanel.Title = 'Shooting Order';
        ShootingOrderPanel.Position = [536 216 129 107];
        % Create ConsecutiveCheckBox
        ConsecutiveCheckBox = uicheckbox(ShootingOrderPanel, 'ValueChangedFcn', @ConsecutiveCheckBoxValueChanged);
        ConsecutiveCheckBox.Text = 'Consecutive';
        ConsecutiveCheckBox.Position = [16 59 98 22];
        % Create ClockwiseCheckBox
        ClockwiseCheckBox = uicheckbox(ShootingOrderPanel, 'ValueChangedFcn', @ClockwiseCheckBoxValueChanged);
        ClockwiseCheckBox.Text = 'Clockwise';
        ClockwiseCheckBox.Position = [16 35 98 22];
        % Create LineSwitch
        LineSwitch = uiswitch(ShootingOrderPanel, 'slider', 'ValueChangedFcn', @LineSwitchValueChanged);
        LineSwitch.Items = {'Points', 'Line'};
        LineSwitch.Position = [48 9 45 20];
        LineSwitch.Value = milldata.line;

        % Create EccentricityPanel
        EccentricityPanel = uipanel(UIFigure);
        EccentricityPanel.TitlePosition = 'centertop';
        EccentricityPanel.Title = 'Eccentricity';
        EccentricityPanel.Position = [536 96 129 102];
        % Create RatioSpinnerLabel
        RatioSpinnerLabel = uilabel(EccentricityPanel);
        RatioSpinnerLabel.HorizontalAlignment = 'right';
        RatioSpinnerLabel.Position = [19 43 34 22];
        RatioSpinnerLabel.Text = 'Ratio';
        % Create RatioSpinner
        RatioSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @RatioSpinnerValueChanged);
        RatioSpinner.Position = [57 43 57 22];
        RatioSpinner.Value = milldata.ratio;
        % Create AngleSpinnerLabel
        AngleSpinnerLabel = uilabel(EccentricityPanel);
        AngleSpinnerLabel.HorizontalAlignment = 'right';
        AngleSpinnerLabel.Position = [15 9 36 22];
        AngleSpinnerLabel.Text = 'Angle';
        % Create AngleSpinner
        AngleSpinner = uispinner(EccentricityPanel, 'ValueChangedFcn', @AngleSpinnerValueChanged);
        AngleSpinner.Position = [56 9 58 22];
        AngleSpinner.Value = milldata.ecc_angle;
        AngleSpinner.Step = 0.05;

        % Create RingDropDownLabel
        RingDropDownLabel = uilabel(UIFigure);
        RingDropDownLabel.HorizontalAlignment = 'right';
        RingDropDownLabel.Position = [44 230 30 22];
        RingDropDownLabel.Text = 'Ring';
        % Create RingDropDown
        RingDropDown = uidropdown(UIFigure, 'ValueChangedFcn', @RingDropDownValueChanged);
        RingDropDown.Items = {'1', '2', '3', '4', '5', '6'};
        RingDropDown.ItemsData = [1 2 3 4 5 6];
        RingDropDown.Position = [108 230 76 22];
        RingDropDown.Value = 1;

        % Create SaveButton
        SaveButton = uibutton(UIFigure, 'push', 'ButtonPushedFcn', @SaveButtonPushed);
        SaveButton.Position = [537 29 129 44];
        SaveButton.Text = 'Save';

        % Create UIAxes
        UIAxes = uiaxes(UIFigure);
        UIAxes.XLim = [-80 80];
        UIAxes.YLim = [-80 80];
        UIAxes.XAxisLocation = 'origin';
        UIAxes.XTick = [-60 -40 -20 0 20 40 60];
        UIAxes.XTickLabel = '';
        UIAxes.YAxisLocation = 'origin';
        UIAxes.YTick = [-60 -40 -20 0 20 40 60];
        UIAxes.YTickLabel = '';
        UIAxes.Position = [197 20 324 310];
        UIAxes.YTickLabel = '';
        UIAxes.Color = 'black';
        plot_milling()

        % Show the UIFigure after all components are created
        UIFigure.Visible = 'on';

        % Construct UIFigure
        function UIFigure = Mill_pattern

            % Create UIFigure and components
            createComponents(UIFigure)

            % Register the UIFigure with UIFigure Designer
            registerUIFigure(UIFigure, UIFigure)

            if nargout == 0
                clear UIFigure
            end
        end

        % Code that executes before UIFigure deletion
        function delete(UIFigure)

            % Delete UIFigure when UIFigure is deleted
            delete(UIFigure)
        end
end