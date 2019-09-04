%% Create a 2D radial FID sequence and export for execution
% 
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
% 
% 1. Create slice selective RF pulse for imaging.
% 2. Create readout gradients in both directions 
% 3. Loop through number of  projections
% 4. Write the sequence to an open file format suitable for execution on a
%    scanner.
%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
Pulseq = uigetdir('Pick your Pulseq source code directory');
addpath(genpath('.'));
addpath(genpath(Pulseq));

system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s');

%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);

%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3;  % Field of view in (m)
Nx=256; Ny=256; % Matrix size
sliceThickness=5e-3; % in (m)
TR = 20e-3; % Repeition time in (s)
TE = 5e-3;%minimum 
dx = fov/Nx;
dy  =dx;
% radp = get_radkparams(dx,dy,fov);
%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences. 
%
% First, a slice selective RF pulse (and corresponding slice gradient) can
% be generated using the |makeSincPulse| function.
%
flip=15*pi/180;
[rf, gz] = mr.makeSincPulse(flip,system,'Duration',1.5e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);

%% Gradients
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$
% 
% Therefore the area of the readout gradient is $n\Delta k$.
deltak=1/fov;
kWidth = Nx*deltak;

readoutTime = 6.4e-3;
gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);

%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.
preTime=8e-4; %Need to figure this one out later!
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',1e-3);
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*preTime);

%% Calculate timing
delayTE=TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf)/2);
%delayTE=TE - mr.calcDuration(gzReph) -(mr.calcDuration(rf)/2)-preTime-gz.fallTime-gx.riseTime;% PP
% delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
%     - mr.calcDuration(gx) - mr.calcDuration(gzSpoil) - delayTE;
delayTR=TR - mr.calcDuration(gzReph) - mr.calcDuration(rf) ...
    - mr.calcDuration(gx) - mr.calcDuration(gzSpoil);
delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTR);

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
Np =256; %Cartesian phase encodes - play around with this number later for evals
Ns = ceil(pi*Np); % Number of spokes required for 360 based on Cartesian Phase encodes
% theta = linspace(0,360,radp.Ns); %This will be replaced with golden angle
%dtheta = 360/Ns;
%% Determining the angle
disp('---------')
disp('1. Normal.')
disp('2. Golden angle.')
disp('3. Tiny golden angle.');
disp('---------')
use_case=input ('Enter the type of angle: ');
switch use_case
    
    case 1
            dtheta = 360/Ns;
            rtheta=deg2rad(dtheta); % Radians
            X = ['Normal angle','=',num2str(dtheta),' deg'];
            disp(X);
    
    case 2
            N=1; %Determines the golden angle (111.24 deg)
            tau=(1+sqrt(5))/2;%  [Refer :Stefan Wundrak et.al., "Golden Ratio Sparse MRI Using Tiny Golden Angles", MRM-2016]
            phiN=pi/(tau+N-1); % [Refer :Stefan Wundrak et.al., "Golden Ratio Sparse MRI Using Tiny Golden Angles", MRM-2016
            phiNd=rad2deg(phiN);% Degree
            rtheta=phiN;
            dtheta=phiNd;
            X = ['Golden angle','=',num2str(phiNd),' deg'];
            disp(X);
            
    case 3
            N=7; % Tiny golden angle (23.6 deg)
            tau=(1+sqrt(5))/2;
            phiN=pi/(tau+N-1); % Radians
            phiNd=rad2deg(phiN);% Degree
            rtheta=phiN;
            dtheta=phiNd;
            X = ['Tiny Golden angle','=',num2str(phiNd),' deg.'];
            disp(X);
end
theta2 = 0:rtheta: (Ns-1)*rtheta; % For debugging
theta=mod(theta2,2*pi); % all angles are in Radians
% theta = get_goldenangle(Ns); % For better temporally resolved data acquisition - not meaningful for FID acq
% Np = length(theta); %equivalent to Ns
ktraj = zeros(Ns, adc.numSamples);
%%
figure(1001); axis([-1000 1000 -1000 1000]);grid on;
pause(10);
for np=1:Ns
            seq.addBlock(rf,gz);
            seq.addBlock(gzReph);
            kWidth_projx = kWidth.*cos(theta(np));
            kWidth_projy = kWidth.*sin(theta(np));

            gx = mr.makeTrapezoid('x',system,'FlatArea',kWidth_projx,'FlatTime',readoutTime);
            gy = mr.makeTrapezoid('y',system,'FlatArea',kWidth_projy,'FlatTime',readoutTime);
            ktraj(np,:) = get_ktraj(gx,gy,adc,1);
            
            %disp(atan2d(gy.flatArea, gx.flatArea));
            seq.addBlock(delay1);
            seq.addBlock(gx,gy,adc);
            seq.addBlock(gzSpoil);
            seq.addBlock(delay2);

end

%% Display the sequence 
ktraj = ktraj./max(abs(ktraj(:)))./2;
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2); 
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);
figure(1002);
seq.plot('TimeRange',[0 2*TR]);
%seq.plot();
fname = ['Rad2D_FID', num2str(Ns),'_',num2str(sliceThickness),'ktraj','_',num2str(dtheta),'deg'];
%uisave('ktrajs', fname );
save([fname,'.mat'], 'ktrajs' );
%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
fname = [fname, '.seq'];
seq.write(fname);

%%
% % Display the first few lines of the output file
% s=fileread('Rad2D_FID.seq');
% disp(s(1:309))
