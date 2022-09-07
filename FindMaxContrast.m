clear all
tic

%%%-------------- Stage preamble --------------%%%
if ( ~exist ( 'AtCube', 'var' ) || ~isa ( AtCube, 'AtCube_control' ) )
    AtCube = AtCube_control();
else
    error('There exis the var. "AtCube"')
end

AtCube.turnOn()

%%%-------------- Camera preamble --------------%%%
uEye_camera(0); % initializing camera

%%-------------- Scan z axis for contrast--------------%%%
% DataFolder = ['C:\Users\co2la\Documents\Data\' datestr(now,formatOut) '\ContrastTest'];
centre = 0.031;
start = 0.025; %start position in mm
n=50; %number of frames
stepsize = ((centre-start)*2)/n; %stepsize in mm


range = (n-1)*stepsize; %range in mm
z_list = linspace(start,start+range,n);

for i = 1:size(z_list,2)
    c = z_list(i);
    AtCube.move_z(c);
    z=AtCube.getPosition_z();
    Filename = [num2str(c,'%.8f') 'B.tif'];

for ex = 1:1
    exp_list = [15];
%% set frame rate
FR = uEye_camera(6, 5);
%% set exposure time
uEye_camera(8, exp_list(ex));
pause(0.2); 
cent_pixel=ones(n,4,3,10);
%% Take an image and shape matrix 
% for k = 1:(150/exp_list(ex))
    for k = 1:2
[err, M, w, h] = uEye_camera(1); % fetch image
M = reshape(M+128, 3, w*h);
I = reshape(M, 3, w, h);
I = permute(I, [3 2 1]);
% cent_pixel(i,ex,1,k)=I(size(I,1)/2,size(I,2)/2,1);
% cent_pixel(i,ex,2,k)=I(size(I,1)/2,size(I,2)/2,2);
% cent_pixelbox(i,ex,3,k)=I(size(I,1)/2,size(I,2)/2,3);

sample_c= [570,950]; % sample centre, pixels
sample_a=20; % sample width/height, pixels
Sample = I(sample_c(1)-0.5*sample_a:sample_c(1)+0.5*sample_a,sample_c(2)-0.5*sample_a:sample_c(2)+0.5*sample_a,:);
contrast = permute(mean(mean(Sample,1),2),[1 3 2]);

cent_pixelbox(i,ex,3,k)=contrast(3);
end
% cent_pixelA(i,ex,1)=mean(cent_pixel(i,ex,1,:),4);
% cent_pixelA(i,ex,2)=mean(cent_pixel(i,ex,2,:),4);
cent_pixelAbox(i,ex,3)=mean(cent_pixelbox(i,ex,3,:),4);
end

end
% plot(z_list,normalize(cent_pixel));

% for i=1:1
% figure
% plot(z_list,normalize(cent_pixel(:,1,i)),z_list,normalize(cent_pixel(:,2,i))-2,z_list,normalize(cent_pixel(:,3,i))-3,z_list,normalize(cent_pixel(:,4,i))-4)
% end
for i=3:3
figure
plot(z_list,cent_pixelAbox(:,1,i),z_list,cent_pixelAbox(:,2,i)-0,z_list,cent_pixelAbox(:,3,i)-0)
end
for i=3:3
figure
plot(z_list,normalize(cent_pixelAbox(:,1,i)),z_list,normalize(cent_pixelAbox(:,2,i))-2,z_list,normalize(cent_pixelAbox(:,3,i))-3)
end
uEye_camera(4); % close camera
AtCube.delete()
toc