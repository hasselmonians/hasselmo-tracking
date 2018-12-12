%LKHavParabolaFtestDiffMatrixOrderOct9dEXCELLENTcorrCol.m - This gives
%correct output with the flipping of columns and rows.

%LKHavParabolaFtestDiffMatrixOrderOct9.m - Test a different order to see
%how to write it in paper.

%LKHavParabolaFallFuncsForFIGUREXCELLENTwithNOTES.m - This adds notes about
%the function of the network.
%FASTER just plots the full array once for each trajectory (not
%sequentially within the prediction loop).

%LKHatWorksForSineCosStraightLineFixedBugEXCELLENTsaveThis.m - Got it to
%work for straight lines or sine waves in either dimension (problem was
%that 3x3 for Wy was being computed AFTER the 2xvX2 option.  Now it works for
%straight lines!

%This script plots a trajectory based on functions of x and y.  It uses the
%corner representation to extra the velocity between sequential points, and
%then the acceleration between sequential velocities. Then it creates a
%position, velocity, acceleration array to compute the dynamics of the
%network and predict the equation. Actual function is shown with red
%circles. Predicted function is shown with blue asterisks. The network can
%accurately predict sine, cosine, circle, straight lines, and quadratic
%polynomials. It cannot do cubic.  Exponentials work but only if the
%starting values are large enough (it runs into rounding errors because the
%matrix is being raised to ^tt).
close all
clear all
figure('Position',[10 600 400 1200]);
fig2D=gcf; axvX2D=gca; colormap(gray); 

%Code for selecting subplots
numplots=1;
h=zeros(1,numplots);
for k=1:numplots %number of subplots
    h(k)=subplot(1,numplots,k);
end

tic
TT=100; %Number of steps/elements for each trajectory
step=0.1 %Step size in time
time=[step:step:step*TT]'; %Time scale set by step size and duration

limit=6; %Step number for starting predictions. Can't be lower than 6, but 
%works with higher numbers (more empty, unpredicted circles)


axes(h(1));

%%MAIN LOOP RUNS THROUGH DIFFERENT TRAJECTORY SHAPES
for func=1:8; %1=parabola, 2=sine, 3=circle (max=8)

    %CREATE THE ARRAYS FOR THE INTERNAL AGENT
xv=zeros(TT,6);
predx=zeros(TT,3);
predy=zeros(TT,3);
corxva=zeros(TT,3);
coryva=zeros(TT,3);
corner=ones(3,3); %Corner will be the "object" that follows the trajectory

%%THIS SECTION CREATES THE INPUT FUNCTIONS THAT ARE BEING PREDICTED.
if func==1
    for tt=1:TT
        xv(tt,1)=10*time(tt); %cos(time(tt)); %
        xv(tt,4)=-40 + 100*time(tt) - 10*time(tt).^2; %sin(time(tt));
    end
elseif func==2
    xv(:,1)=10*time; %
    xv(:,4)=100*sin(time);
elseif func==3
    xv(:,1)=50+50*cos(time); %
    xv(:,4)=50*sin(time);
elseif func==4
    xv(:,1)=1*time.^2; %
    xv(:,4)=60*cos(time);
elseif func==5
    xv(:,1)=10*time; %
    xv(:,4)=100*exp(time/20).*sin(time);
elseif func==6;
    xv(:,1)=50+50*exp(-time/10).*sin(time);
    xv(:,4)=-20 + 100*exp(-time/10).*cos(time);
elseif func==7;
    xv(:,1)=10*time; %
    xv(:,4)=4*time;
elseif func==8
    xv(:,1)=50+40*cos(time); %
    xv(:,4)=200-30*time;
end
%END CREATION OF TRAJECTORIES


%PREDICT THE TRAJECTORIES

for tt=1:TT-1;
%UPDATE CORNER that reflects three features on the object based on current
%position,velocity xv.
oldcorner=corner;
corner(1:2,1:3)=[xv(tt+1,[1 4])+[.01 0]; xv(tt+1,[1 4]); xv(tt+1,[1 4])+[0 .01]]';

%COMPUTE affine matrix invW.cornerW with translation=velocity from INVERSE of 
%two consecutive corners (oldcorner and corner). Can only start at tt>1.
if tt>1 && tt<limit+1 %VELOCITY
     invW(tt).cornerW=(corner*inv(oldcorner)); %This gives translation velocity between sequential locations
%         ['x and y velocity'] 
%         invW(tt).cornerW
end
%COMPUTE affine matrix invW.cornerWW with translation=acceleration from 
%two consecutive invW.cornerW at tt and tt-1. This gives acceleration.
if tt>2 && tt<limit+1 %ACCELERATION
    invW(tt).cornerWW=(invW(tt).cornerW*inv(invW(tt-1).cornerW)); %This gives acceleration
%     ['x and y acceleration'] 
%     invW(tt).cornerWW
    
%Put extracted velocity and acceleration into a vector array into temporary arrays corxva 
%that contain position, velocity, acceleration. Need to do this because can't use 1:3 in struct    
corxva(tt,:)=[xv(tt,1) invW(tt).cornerW(1,3) invW(tt).cornerWW(1,3)];
coryva(tt,:)=[xv(tt,4) invW(tt).cornerW(2,3) invW(tt).cornerWW(2,3)];
end

%Now create the sequential x,v,a arrays. Note that this takes the elements
%from corxva because MATLAB can't copy sets of values from struct invW.
if tt>3 && tt<limit+1;
    xvX1=ones(3,3);
    xvX2=ones(3,3);
    xvY1=ones(3,3);
    xvY2=ones(3,3);
    xvX1(1:2,1:3)=corxva(tt-3:tt-1,1:2)'; %make affine array from the x and vel components
    xvX2(1:2,1:3)=corxva(tt-2:tt,1:2)';
    xvY1(1:2,1:3)=coryva(tt-3:tt-1,1:2)'; %make affine array from the y and vel components
    xvY2(1:2,1:3)=coryva(tt-2:tt,1:2)';

%COMPUTE the inverse PREDICTION matrix from the sequential created x,v,a array xvX1 and xvX2.
    invW(tt).Wx=(xvX2*inv(xvX1)); %3x3 prediction matrix for x dynamics (if curved)
    if abs(corxva(tt,3))<0.00001;  %This is to use 2x2 matrix if linear trajectory (no acceleration)
        invW(tt).Wx=eye(3,3);
        invW(tt).Wx(1:2,1:2)=(xvX2(1:2,1:2)*inv(xvX1(1:2,1:2))); %X dimension 2x2
    end
 %['x pos,xv, xa and Wx transition with air resistance']
 %[xvY1 invW(tt).Wx]
 
    invW(tt).Wy=(xvY2*inv(xvY1)); %3x3 prediction matrix for y dynamics (if curved)
    invW(tt).WyDiff=xvY2*inv(xvY1); %Y dimension 2x2
 
    if abs(coryva(tt,3))<0.00001;  %This is to use 2x2 matrix if linear trajectory (no acceleration)
        invW(tt).Wy=eye(3,3);
        invW(tt).Wy(1:2,1:2)=(xvY2(1:2,1:2)*inv(xvY1(1:2,1:2))); %Y dimension 2x2
        
        
    end
%  ['y pos, yv, ya and Wy transition']
%  [xvX1 invW(tt).Wy]
end %if t>3
end %for tt initial

for tt=1:TT;
%ADD PREDICTION CAPABILITIES based on derived matrix
graphdist=240;
if tt>limit %Wy is only available after tt
    predx(tt,:)=(invW(limit).Wx^(tt-(limit))*[corxva(limit,1:2) 1]')';
    predy(tt,:)=(invW(limit).Wy^(tt-(limit))*[coryva(limit,1:2) 1]')';
%        invW(limit).Wx %^(tt-(limit))
%        invW(limit).Wy %^(tt-(limit))
end



%end %for tt
if (tt<TT && func==1) || (func>1 && mod(tt,4)==1);
plot(xv(:,1),xv(:,4)-graphdist*func,'o','MarkerSize',8,'Color',[0.5 0.5 0.5]); hold on
plot(predx(:,1),predy(:,1)-graphdist*func,'k*'); hold on; 
xlim([-10 110]); ylim([-2100 100]);
drawnow
end

end %for tt


invW(limit).Wx
invW(limit).Wy
invW(limit).WyDiff

[inv(xvY1)  xvY2] 
[(inv(xvY1)*xvY2)  (inv(xvY1)*xvY2)']
[(xvY2'*inv(xvY1'))  ]
%pause

end %for func=1:7
toc
