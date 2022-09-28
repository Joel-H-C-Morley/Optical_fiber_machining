function [] = D2_fit_GUI(recon_surface)
    

% Create D2fit_UIFigure and hide until all components are created
    D2fit_UIFigure = uifigure('Visible', 'off');
    D2fit_UIFigure.AutoResizeChildren = 'off';
    D2fit_UIFigure.Position = [200 100 901 805];
    D2fit_UIFigure.Name = '2D Fitting';

    fitdata = guidata(D2fit_UIFigure);
    setting_file = 'fit_settings.mat';
    load(setting_file, 'fitdata');
    if ~isempty(recon_surface)
    else
        fitdata.surface = recon_surface
    end
    guidata(D2fit_UIFigure,fitdata);

    %% Functions
    
    % Button pushed function: OpenButton
    function OpenButtonPushed(D2fit_UIFigure, event)
        fitdata = guidata(D2fit_UIFigure);
        [file,path] = uigetfile('C:\Users\Joel Morley\Documents\20220906','*.mat')
        surf_data = load(fullfile(path,file));
        FileEditField.Value = fullfile(path,file);
        fitdata.surface = cell2mat(struct2cell(surf_data));
        guidata(D2fit_UIFigure,fitdata);
    end

    function SaveButtonPushed(D2fit_UIFigure, event)
    % to be completed depending on requirements
    end

    function RangeumEditFieldValueChanged(D2fit_UIFigure, event)
        fitdata = guidata(D2fit_UIFigure);
        fitdata.range = RangeumEditField.Value;
        guidata(D2fit_UIFigure,fitdata);
    end

    function FitButtonPushed(D2fit_UIFigure, event)
        %% Initialise variables and matrices
        fitdata = guidata(D2fit_UIFigure);
        Z = fitdata.surface*1000; %height : mm to um
        [X,Y] = meshgrid(-(size(Z,1)-1)/2:(size(Z,2)-1)/2);
        xdata = zeros(size(X,1),size(Y,2),2);
        pixel=0.329; % Pixel FOV in um
        xdata(:,:,1) = X*pixel;%in um
        xdata(:,:,2) = Y*pixel;%in um

        widthpx = round(fitdata.range/pixel); %range for fitting (pixel)
        [Xfit,Yfit] = meshgrid(-widthpx:widthpx); %mesh for fitting
        xdata_fit = zeros(widthpx*2+1,widthpx*2+1,2);
        d2 = (Xfit).^2 + (Yfit).^2 <= widthpx.^2; %remove values outside of fitting area
        D2 = double(d2);
        D2(D2<1)=nan;
        xdata_fit(:,:,1) = Xfit*pixel.*D2;%fitting range in um, exclude outside
        xdata_fit(:,:,2) = Yfit*pixel.*D2;%fitting range in um, exclude outside

        %% Fitting
        Z_fit = Z((size(Z,1)-1)/2-widthpx:(size(Z,1)-1)/2+widthpx,(size(Z,1)-1)/2-widthpx:(size(Z,1)-1)/2+widthpx).*D2; %crop z to fitting range
        Zdata = Z_fit(~isnan(Z_fit));

        % initial fit parameters and limits for Gaussian
        xg0 = [1            ,0   ,500         ,0   ,500         , 0, 10]; %xg0 = [Amp,x0,wx,y0,wy,angle,fi]; Inital guess parameters for gaussian fit
        lbg = [0.01*xg0(1) ,-50 ,0.01*xg0(3) ,-50 ,0.01*xg0(5) , -pi/4, 0.01*xg0(7)]; %lower boundry for gaussian fit
        ubg = [100*xg0(1)   ,50  ,5*xg0(3)    ,50  , 5*xg0(5)   ,pi/4 ,1000*xg0(7)]; %upper boundry for gaussian fit
        [xg,resnormg,residualg,exitflagg] = lsqcurvefit(@D2GaussFunctionRot, xg0, xdata_fit, Zdata, lbg, ubg);
        RfGx = (xg(3)^2)/xg(1);%radius calculated from gaussian fit
        RfGy = (xg(5)^2)/xg(1);%radius calculated from gaussian fit
        varg = var(residualg(:),1); %calculate distribution
        difg = Z_fit - D2GaussFunctionRot2(xg,xdata_fit); %residual
        figlength = length(difg);
        c_axis = round(figlength/2); %get index for central axis
        fitdata.rot = xg(6)*180/pi();
        fitdata.ecc = xg(3)/xg(5);

        %% Plotting
        D1_fit = D2GaussFunctionRot2(xg,xdata_fit);

        X_G = @(ps,x) -ps(1)*exp(-((x-ps(2))/ps(3)).^2/2)+ps(7); % Gaussian function to fit the data to
        Y_G = @(ps,x) -ps(1)*exp(-((x-ps(4))/ps(5)).^2/2)+ps(7); % Gaussian function to fit the data to
        Y_r=-(xg(3)^2)/-xg(1);%um
        X_r=-(xg(5)^2)/-xg(1);%um
        Y_circle_cent_x=xg(2);
        X_circle_cent_x=xg(4);
        X_circle_cent_z=min(D1_fit(c_axis,:))+X_r;
        Y_circle_cent_z=min(D1_fit(:,c_axis))+Y_r;
        X_circle_limit_low = (3*pi/2)-(2*atan((Z_fit(c_axis,1)-min(Z_fit(c_axis,:)))/fitdata.range));
        X_circle_limit_high = (3*pi/2)+(2*atan((Z_fit(c_axis,end)-min(Z_fit(c_axis,:)))/fitdata.range));
        Y_circle_limit_low = (3*pi/2)-(2*atan((Z_fit(1,c_axis)-min(Z_fit(:,c_axis)))/fitdata.range));
        Y_circle_limit_high = (3*pi/2)+(2*atan((Z_fit(end,c_axis)-min(Z_fit(:,c_axis)))/fitdata.range));
        Y_theta = linspace(Y_circle_limit_low ,Y_circle_limit_high,figlength);
        X_theta = linspace(X_circle_limit_low ,X_circle_limit_high,figlength);
        cla(CS_x_axes)
        cla(CS_y_axes)
        hold(CS_x_axes,'on')
        plot(CS_x_axes, X_circle_cent_x + X_r*cos(X_theta),X_circle_cent_z + X_r*sin(X_theta),'.-')
        plot(CS_x_axes, xdata_fit(c_axis,:,1), D1_fit(c_axis,:), xdata_fit(c_axis,:,1), Z_fit(c_axis,:),'+')
        hold(CS_x_axes,'off')
        hold(CS_y_axes,'on')
        plot(CS_y_axes, xdata_fit(c_axis,:,1), D1_fit(:,c_axis), xdata_fit(c_axis,:,1), Z_fit(:,c_axis), '+')
        plot(CS_y_axes, Y_circle_cent_x + Y_r*cos(Y_theta),Y_circle_cent_z + Y_r*sin(Y_theta),'.-')
        hold(CS_y_axes,'off')
        fitdata.xr = X_r;
        fitdata.yr = Y_r;
        fitdata.xc = xg(2);
        fitdata.yc = xg(4);
        RotationradEditField.Value = fitdata.rot;
        EccentricityratioEditField.Value = fitdata.ecc;
        XCenterumEditField.Value = fitdata.xc;
        YCenterumEditField.Value = fitdata.yc;
        RoCxumEditField.Value = fitdata.xr;
        RoCyumEditField.Value = fitdata.yr;

        xa = linspace(-fitdata.range,fitdata.range,figlength);
        ya = linspace(-fitdata.range,fitdata.range,figlength);
        imagesc(Residual_axes, xa, ya, difg);
        strg = ['Residual for Gaussian fit : ' num2str(varg)];
        title(Residual_axes, strg)
        
        contour3(Surface_axes, xdata(:,:,1),xdata(:,:,2),Z, 30, 'k', 'LineWidth', 1.5) %plot data
%         hold(Surface_axes,'off')
%         title('Reconstructed Contours with fitted surface')
        surface(Surface_axes, xdata_fit(:,:,1),xdata_fit(:,:,2),D2GaussFunctionRot2(xg,xdata_fit),'EdgeColor','none') %plot fit
        axis(Surface_axes, [-(fitdata.range+20) (fitdata.range+20) -(fitdata.range+20) (fitdata.range+20) -0.5 2 ])
        view(Surface_axes,90,0)
        daspect(Surface_axes, [1 1 0.05])
%         hold(Surface_axes,'off')
        guidata(D2fit_UIFigure,fitdata);
        
    end
    
    function Fg = D2GaussFunctionRot(x,xdata)
        xaxisrot= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
        yaxisrot= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
        x0rot = x(4)*cos(x(6)) - x(2)*sin(x(6));
        y0rot = x(4)*sin(x(6)) + x(2)*cos(x(6));
        xaxisrot = xaxisrot (~isnan (xaxisrot));
        yaxisrot = yaxisrot (~isnan (yaxisrot));
        L = length(xaxisrot);
        xyaxis = zeros(L,2);
        xyaxis(:,1) = xaxisrot;
        xyaxis(:,2) = yaxisrot;
        Fg = -x(1)*exp(   -((xyaxis(:,1)-x0rot).^2/(2*x(3)^2) + (xyaxis(:,2)-y0rot).^2/(2*x(5)^2) )    )+x(7);
    end

    function Fg2 = D2GaussFunctionRot2(x,xdata)
        xdatarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
        xdatarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
        x0rot = x(4)*cos(x(6)) - x(2)*sin(x(6));
        y0rot = x(4)*sin(x(6)) + x(2)*cos(x(6));
        Fg2 = -x(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    )+x(7);
    end
    %% GUI Layout

            % Create InputPanel
            InputPanel = uipanel(D2fit_UIFigure);
            InputPanel.TitlePosition = 'centertop';
            InputPanel.Title = 'Input';
            InputPanel.FontWeight = 'bold';
            InputPanel.FontSize = 14;
            InputPanel.Position = [51 712 423 75];

            % Create OpenButton
            OpenButton = uibutton(InputPanel, 'push', 'ButtonPushedFcn', @OpenButtonPushed);
            OpenButton.Position = [30 12 100 22];
            OpenButton.Text = 'Open';

            % Create FileEditFieldLabel
            FileEditFieldLabel = uilabel(InputPanel);
            FileEditFieldLabel.HorizontalAlignment = 'right';
            FileEditFieldLabel.Position = [156 12 25 22];
            FileEditFieldLabel.Text = 'File';

            % Create FileEditField
            FileEditField = uieditfield(InputPanel, 'text');
            FileEditField.Position = [196 12 208 22];
            FileEditField.Value = 'Current Reconstruction';

            % Create ResultsPanel
            ResultsPanel = uipanel(D2fit_UIFigure);
            ResultsPanel.TitlePosition = 'centertop';
            ResultsPanel.Title = 'Results';
            ResultsPanel.FontWeight = 'bold';
            ResultsPanel.FontSize = 14;
            ResultsPanel.Position = [51 27 423 129];

            % Create SaveFitButton
            SaveFitButton = uibutton(ResultsPanel, 'push', 'ButtonPushedFcn', @SaveButtonPushed);
            SaveFitButton.Position = [50 9 100 22];
            SaveFitButton.Text = 'Save Fit';

            % Create RoCxumEditFieldLabel
            RoCxumEditFieldLabel = uilabel(ResultsPanel);
            RoCxumEditFieldLabel.HorizontalAlignment = 'right';
            RoCxumEditFieldLabel.Position = [27 68 42 28];
            RoCxumEditFieldLabel.Text = {'RoC, x'; ' (um)'};

            % Create RoCxumEditField
            RoCxumEditField = uieditfield(ResultsPanel, 'numeric');
            RoCxumEditField.ValueDisplayFormat = '%.0f';
            RoCxumEditField.Position = [84 74 37 22];
            RoCxumEditField.Value = fitdata.xr;

            % Create RoCyumEditFieldLabel
            RoCyumEditFieldLabel = uilabel(ResultsPanel);
            RoCyumEditFieldLabel.HorizontalAlignment = 'right';
            RoCyumEditFieldLabel.Position = [27 37 42 28];
            RoCyumEditFieldLabel.Text = {'RoC, y'; ' (um)'};

            % Create RoCyumEditField
            RoCyumEditField = uieditfield(ResultsPanel, 'numeric');
            RoCyumEditField.ValueDisplayFormat = '%.0f';
            RoCyumEditField.Position = [84 43 37 22];
            RoCyumEditField.Value = fitdata.yr;

            % Create RotationradEditFieldLabel
            RotationradEditFieldLabel = uilabel(ResultsPanel);
            RotationradEditFieldLabel.HorizontalAlignment = 'right';
            RotationradEditFieldLabel.Position = [292 68 50 28];
            RotationradEditFieldLabel.Text = {'Rotation'; '(deg)'};

            % Create RotationradEditField
            RotationradEditField = uieditfield(ResultsPanel, 'numeric');
            RotationradEditField.ValueDisplayFormat = '%.1f';
            RotationradEditField.Position = [357 74 37 22];
            RotationradEditField.Value = fitdata.rot;

            % Create EccentricityratioEditFieldLabel
            EccentricityratioEditFieldLabel = uilabel(ResultsPanel);
            EccentricityratioEditFieldLabel.HorizontalAlignment = 'right';
            EccentricityratioEditFieldLabel.Position = [275 37 67 28];
            EccentricityratioEditFieldLabel.Text = {'Eccentricity'; 'ratio'};

            % Create EccentricityratioEditField
            EccentricityratioEditField = uieditfield(ResultsPanel, 'numeric');
            EccentricityratioEditField.ValueDisplayFormat = '%.2f';
            EccentricityratioEditField.Position = [357 43 37 22];
            EccentricityratioEditField.Value = fitdata.ecc;

            % Create XCenterumEditFieldLabel
            XCenterumEditFieldLabel = uilabel(ResultsPanel);
            XCenterumEditFieldLabel.HorizontalAlignment = 'right';
            XCenterumEditFieldLabel.Position = [156 68 53 28];
            XCenterumEditFieldLabel.Text = {'X Center'; '(um)'};

            % Create XCenterumEditField
            XCenterumEditField = uieditfield(ResultsPanel, 'numeric');
            XCenterumEditField.ValueDisplayFormat = '%.0f';
            XCenterumEditField.Position = [223 74 37 22];
            XCenterumEditField.Value = 0;

            % Create YCenterumEditFieldLabel
            YCenterumEditFieldLabel = uilabel(ResultsPanel);
            YCenterumEditFieldLabel.HorizontalAlignment = 'right';
            YCenterumEditFieldLabel.Position = [156 37 53 28];
            YCenterumEditFieldLabel.Text = {'Y Center'; '(um)'};

            % Create YCenterumEditField
            YCenterumEditField = uieditfield(ResultsPanel, 'numeric');
            YCenterumEditField.ValueDisplayFormat = '%.0f';
            YCenterumEditField.Position = [223 43 37 22];
            YCenterumEditField.Value = 0;
            
            % Create FitButton
            FitButton = uibutton(D2fit_UIFigure, 'push', 'ButtonPushedFcn', @FitButtonPushed);
            FitButton.Position = [295 663 112 31];
            FitButton.Text = 'Fit';

            % Create RangeumEditFieldLabel
            RangeumEditFieldLabel = uilabel(D2fit_UIFigure);
            RangeumEditFieldLabel.HorizontalAlignment = 'right';
            RangeumEditFieldLabel.Position = [171 663 44 28];
            RangeumEditFieldLabel.Text = {' Range'; ' (um)'};

            % Create RangeumEditField
            RangeumEditField = uieditfield(D2fit_UIFigure, 'numeric', 'ValueChangedFcn', @RangeumEditFieldValueChanged);
            RangeumEditField.Limits = [0.001 100];
            RangeumEditField.Position = [227 669 53 22];
            RangeumEditField.Value = fitdata.range;

            % Create CS_x_axes
            CS_x_axes = uiaxes(D2fit_UIFigure);
            title(CS_x_axes, 'Cross Section, x')
            xlabel(CS_x_axes, 'x (um)')
            ylabel(CS_x_axes, 'Height (um)')
            CS_x_axes.PlotBoxAspectRatio = [2.30666666666667 1 1];
            CS_x_axes.BoxStyle = 'full';
            CS_x_axes.ClippingStyle = 'rectangle';
            CS_x_axes.Clipping = 'off';
            CS_x_axes.Position = [20 410 454 265];

            % Create CS_y_axes
            CS_y_axes = uiaxes(D2fit_UIFigure);
            title(CS_y_axes, 'Cross Section, y')
            xlabel(CS_y_axes, 'y (um)')
            ylabel(CS_y_axes, 'Height (um)')
            zlabel(CS_y_axes, 'Z')
            CS_y_axes.PlotBoxAspectRatio = [2.29139072847682 1 1];
            CS_y_axes.Clipping = 'off';
            CS_y_axes.Position = [20 164 454 265];

            % Create Residual_axes
            Residual_axes = uiaxes(D2fit_UIFigure);
            title(Residual_axes, 'Residual for Gaussian fit')
            xlabel(Residual_axes, 'x (um)')
            ylabel(Residual_axes, 'y (um)')
            zlabel(Residual_axes, 'Z')
            Residual_axes.DataAspectRatio = [1 1 1];
            Residual_axes.PlotBoxAspectRatio = [1 1 1];
            Residual_axes.Clipping = 'off';
            Residual_axes.Position = [491 412 400 380];

            % Create Surface_axes
            Surface_axes = uiaxes(D2fit_UIFigure);
            title(Surface_axes, 'Reconstructed Contours with fitted surface')
            xlabel(Surface_axes, 'x (um)')
            ylabel(Surface_axes, 'y (um)')
            zlabel(Surface_axes, 'Height (um)')
            Surface_axes.DataAspectRatio = [1 1 1];
            Surface_axes.PlotBoxAspectRatio = [1 1 1];
            Surface_axes.Clipping = 'on';
            Surface_axes.Position = [491 8 400 380];

        % Show the UIFigure after all components are created
        D2fit_UIFigure.Visible = 'on';

        % Construct UIFigure
        function D2fit_UIFigure = Mill_pattern

            % Create UIFigure and components
            createComponents(D2fit_UIFigure)

            % Register the UIFigure with UIFigure Designer
            registerUIFigure(D2fit_UIFigure, D2fit_UIFigure)

            if nargout == 0
                clear UIFigure
            end
        end

        % Code that executes before UIFigure deletion
        function delete(D2fit_UIFigure)

            % Delete UIFigure when UIFigure is deleted
            delete(D2fit_UIFigure)
        end
end