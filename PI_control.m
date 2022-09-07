classdef PI_control 
    %PI_control class is for control the PI C-885 MotionMaster easily.
    %
    properties (SetAccess = protected)
        PI_Controller
        PIdevice
        axis_x = '1';
        axis_y = '5';
        axis_z = '3';
        xyzMinMax = [80 155 0 13 0 25]
    end
    
    methods
        function obj = PI_control()
            % PI_control Construct an instance of this class
            % Here they cheack the information then try to connect via
            % TCPIP and initialize them
            
            obj.PI_Controller = PI_GCS_Controller ();
            %class(obj.PI_Controller)
            
            ip = '192.168.0.3';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
            port = 50000;           % Is 50000 for almost all PI controllers
            Controller = obj.PI_Controller;
            obj.PIdevice = Controller.ConnectTCPIP ( ip, port ) ;
        
            % Query controller identification string
            %obj.PIdevice.qIDN()
            %class(obj.PIdevice)
            % initialize PIdevice object for use in MATLAB
            %PIdevice = obj.PIdevice.InitializeController ();
        end        
        
        function turnOnXYZ(obj)
            % switch servo on for axis
            switchOn    = 1;
            switchOff   = 0;
            obj.PIdevice.SVO ( obj.axis_x, switchOn );
            obj.PIdevice.SVO ( obj.axis_y, switchOn );
            obj.PIdevice.SVO ( obj.axis_z, switchOn );

            % reference axis
            obj.PIdevice.FRF ( obj.axis_x );  % find reference
            obj.PIdevice.FRF ( obj.axis_y );  % find reference
            obj.PIdevice.FRF ( obj.axis_z );  % find reference
            bReferencing = 1;            
            fprintf(1,'Stage is referencing\n')
            % wait for referencing to finish
            while(0 ~= obj.PIdevice.qFRF ( obj.axis_x ) == 0 )                        
                pause(0.1);           
                %[char ( 8 ) * ones( 1, 7 ), '.']
            end
            fprintf(1, 'x axis was turned on\n');
            while(0 ~= obj.PIdevice.qFRF ( obj.axis_y ) == 0 )                        
                pause(0.1);           
                %[char ( 8 ) * ones( 1, 7 ), '.']
            end
            fprintf(1,'y axis was turned on\n')
            while(0 ~= obj.PIdevice.qFRF ( obj.axis_z ) == 0 )                        
                pause(0.1);           
                %[char ( 8 ) * ones( 1, 7 ), '.']
            end
            fprintf(1,'z axis was turned on\n')
        end
%%%------------------Move x axis------------------%%%        
        function move_x(obj, x)
            if ( (x < obj.xyzMinMax(1)) || (x > obj.xyzMinMax(2)))
                disp(' Commanded position out of range!')
            else
                obj.PIdevice.MOV ( obj.axis_x, x );
                % disp ( 'Stage is moving')
                % wait for motion to stop
                while(0 ~= obj.PIdevice.IsMoving ( obj.axis_x ) )
                    pause ( 0.1 );
                    %[char ( 8 ) * ones( 1, 7 ), '.']
                end
            end
        positonReached = obj.PIdevice.qPOS(obj.axis_x);
        end
%%%------------------Move y axis------------------%%%          
        function move_y(obj, y)
            if ( (y < obj.xyzMinMax(3)) || (y > obj.xyzMinMax(4)))
                disp(' Commanded position out of range!')
            else
                obj.PIdevice.MOV ( obj.axis_y, y );
                disp ( 'Stage is moving')
                % wait for motion to stop
                while(0 ~= obj.PIdevice.IsMoving ( obj.axis_y ) )
                    pause ( 0.1 );
                    %[char ( 8 ) * ones( 1, 7 ), '.']
                end
            end
        positonReached = obj.PIdevice.qPOS(obj.axis_y);
        end
%%%------------------Move z axis------------------%%%                  
        function move_z(obj, z)
            if ( (z < obj.xyzMinMax(5)) || (z > obj.xyzMinMax(6)))
                disp(' Commanded position out of range!')
            else
                obj.PIdevice.MOV ( obj.axis_z, z );
                % disp ( 'Stage is moving')
                % wait for motion to stop
                while(0 ~= obj.PIdevice.IsMoving ( obj.axis_z ) )
                    pause ( 0.1 );
                    %[char ( 8 ) * ones( 1, 7 ), '.']
                end
            end
        positonReached = obj.PIdevice.qPOS(obj.axis_z);
        end
%%%------------------Get x, y, z position------------------%%%
        function x = getPosition_x(obj)
            x = obj.PIdevice.qPOS(obj.axis_x);
        end
        
        function y = getPosition_y(obj)
            y = obj.PIdevice.qPOS(obj.axis_y);
        end
        
        function y = getPosition_z(obj)
            y = obj.PIdevice.qPOS(obj.axis_z);
        end
        
        function delete(obj)
            obj.PIdevice.CloseConnection();
            obj.PI_Controller.Destroy();
            clear PI_Controller;
            clear PIdevice; 
            fprintf(1, 'Object was deleted\n')
        end
    end
end
