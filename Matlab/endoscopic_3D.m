clc;
clear all;
close all;
%% step1: Camera calibration

% Setting calibration board
num_x = 11; % number of circles in the x direction of the calibration board
num_y = 9; % number of circles in the y direction of the calibration board
dist_circ = 25; % The distance between two adjacent circles on the calibration board
disp('Start camera calibration...');

% Load the pixel coordinates of the corner points on the calibration board detected in the camera image
% These coordinates are stored in the camera_imagePoints.mat file, downloaded from the online open source file
load('camera_imagePoints.mat'); % load the camera image points of the centers

% The world coordinate system of corner points
worldPoints = generateCheckerboardPoints([num_y+1,num_x+1], dist_circ);

% The pixel coordinates in the image and the generated world coordinates are used to calibrate the camera
% The estimateCameraParameters function is used for camera calibration
% The internal and external parameters and calibration errors are obtained.
[cameraParams,imagesUsed,estimationErrors] =estimateCameraParameters(imagePoints,worldPoints, 'EstimateTangentialDistortion', true);

% Displays a reprojection error map of the camera calibration to evaluate the accuracy of the calibration.
figure(1);showReprojectionErrors(cameraParams);title('Reprojection error in camera calibration');
% Generates external parameters for the camera calibration, including the rotation and shift matrices of the camera.
figure(2);showExtrinsics(cameraParams);title('External parameters for camera calibration'); 

% save Rc_1, Tc_1, KK, inv_KK,
% Pose of the camera with respect to the origin of the world coordinate system (rotation + translation = external parameter)
Rc_1 = cameraParams.RotationMatrices(:,:,1); % Rotation matrix of the camera，it's a three-dimensional array
Rc_1 = Rc_1'; % transpose
Tc_1 = cameraParams.TranslationVectors(1, :); % The translation vector of the camera，it's a two-dimensional array
Tc_1 = Tc_1'; % transpose

% camera internal parameter
% CameraParams. IntrinsicMatrix is get from the camera calibration process
% It contains such as focal length, main point position, and distortion parameters.
KK = cameraParams.IntrinsicMatrix';
save('CamCalibResult.mat', 'Rc_1', 'Tc_1', 'KK');

%% step2: Projector calibration
% Similar to the process of camera calibration

disp('Start projector calibration...');

% Load the pixel coordinates of the corner points on the calibration board detected in the projector image
% These coordinates are stored in the projector_imagePoints.mat file, downloaded from the online open source file
load('projector_imagePoints.mat'); % load the projector image points of the centers

% The world coordinate system of corner points
worldPoints = generateCheckerboardPoints([num_y+1,num_x+1], dist_circ);
[prjParams,imagesUsed,estimationErrors] = estimateCameraParameters(prjPoints,worldPoints,'EstimateTangentialDistortion', true);
figure(3); showReprojectionErrors(prjParams); title('Reprojection error of the projector');
figure(4); showExtrinsics(prjParams);title('External parameters of the projector');

% save parameter(rotation + translation = external parameter)
Rc_1 = prjParams.RotationMatrices(:,:,1);
Rc_1 = Rc_1';
Tc_1 = prjParams.TranslationVectors(1, :);
Tc_1 = Tc_1';
KK = prjParams.IntrinsicMatrix';
save('PrjCalibResult.mat', 'Rc_1', 'Tc_1', 'KK')

%% step3: input parameters

% Set parameters according to the size of the test picture and the specifications of the hardware device.
% Make sure the coordinates match between the camera, the projector and the test picture
width = 640; % camera width
height = 480; % camera height
prj_width = 912; % projector width

%camera: Projection matrix Pc
load('CamCalibResult.mat');
Kc = KK; % camera internal parameter
Ac = Kc * [Rc_1, Tc_1]; % Calculates the camera projection matrix

%projector: Projection matrix Ap
load('PrjCalibResult.mat');
Kp = KK; % projector internal parameter
Ap = Kp * [Rc_1, Tc_1]; % Calculates the projection matrix of the projector


% The fringe frequency 64, which is also the spacing (one cycle consists of 64 pixels), is used to calculate the absolute phase
% The other frequencies 1 and 8 are used to wrap the phase unfolding
%Corresponding to the phase information of the test picture
f = 64; % Fringe frequency
load('up_test_obj.mat');
up_test_obj = up_test_obj / f; % Normalize the phase to between [0, 2pi]
figure; imshow(up_test_obj / (2 * pi)); % The normalized phase image is displayed
colorbar; % Displays the color bar
title("Phase diagram, freq=" + num2str(f));
figure; mesh(up_test_obj); % Displays phase images in the form of a three-dimensional grid
colorbar;
title("Three-dimensional grid phase diagram, freq=" + num2str(f));

% Calculates the projector coordinates
x_p = up_test_obj / (2 * pi) * prj_width;

%% step4: 3D reconstruction 

% Initialize the matrices Xws, Yws, and Zws for storing three-dimensional world coordinates.
Xws = nan(height, width);
Yws = nan(height, width);
Zws = nan(height, width);

% Start a nested loop, iterating through each pixel in a two-dimensional image
for y = 1:height
for x = 1:width
% Check whether the pixels at coordinates (y, x) in the up_test_obj matrix contain valid phase values
if ~isnan(up_test_obj(y, x))
% Calculate the adjusted pixel coordinates uc, vc, and up by subtracting 1
uc = x - 1;
vc = y - 1;
up = (x_p(y, x) - 1);


% Eq. (32) in the reference paper (Calibration of fringe projection profilometry: A comparative review)
% Eq. (18) in the reference paper（GPU-assisted high-resolution, real-time 3-D shape measurement）
% The matrix A is constructed using the internal and external parameters of the camera and projector (Ac, Ap) and the adjusted pixel coordinates
A = [Ac(1,1) - Ac(3,1) * uc, Ac(1,2) - Ac(3,2) * uc, Ac(1,3) - Ac(3,3) * uc;
Ac(2,1) - Ac(3,1) * vc, Ac(2,2) - Ac(3,2) * vc, Ac(2,3) - Ac(3,3) * vc;
Ap(1,1) - Ap(3,1) * up, Ap(1,2) - Ap(3,2) * up, Ap(1,3) - Ap(3,3) * up];

% Construct vector b using the translation portion of the external parameters of the camera and projector and the adjusted pixel coordinates
b = [Ac(3,4) * uc - Ac(1,4);
Ac(3,4) * vc - Ac(2,4);
Ap(3,4) * up - Ap(1,4)];
% Solve linear equations by matrix inversion and multiplication, obtain world coordinates (X, Y, Z)
XYZ_w = inv(A) * b;
% The resulting world coordinates are assigned to the corresponding elements in the Xws, Yws, and Zws matrices
Xws(y, x) = XYZ_w(1); 
Yws(y, x) = XYZ_w(2); 
Zws(y, x) = XYZ_w(3);
end
end
end

%% step4: Point cloud display

% Merge the three-dimensional coordinate data from the Xws, Yws, and Zws matrices into a matrix named xyzPoints
xyzPoints(:, 1) = Xws(:);
xyzPoints(:, 2) = Yws(:);
xyzPoints(:, 3) = Zws(:);
% A point cloud object, ptCloud, is created using the pointCloud function
% Contains the three-dimensional coordinate information from xyzPoints.
ptCloud = pointCloud(xyzPoints);
% Calculate the coordinate range of the point cloud in the X, Y, and Z directions
xlimits = [min(Xws(:)), max(Xws(:))];
ylimits = [min(Yws(:)), max(Yws(:))];
zlimits = ptCloud.ZLimits;
% A point cloud player object, player, is created using the pcplayer function and the initial axis range is set
player = pcplayer(xlimits,ylimits,zlimits);
% Sets the axis label of the point cloud player object
xlabel(player.Axes,'X (mm)');
ylabel(player.Axes,'Y (mm)');
zlabel(player.Axes,'Z (mm)');
% Use the view function to display point cloud data in the point cloud player
view(player,ptCloud);


