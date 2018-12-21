%LVarSeparatePlotForPredEXCELLENTshowsPredTraj.m - This generates predicted
%trajectories and plots them on the right as the network progresses with
%larger ranges.  Many of the predictions are very bad, but some look good.
%However, looks like it might need more interpolation.

%LVamFocusOnRangesEtakeFromvvnextEXCELLENT.m - This seems to be effectively
%computing the acceleration (usually 1 because of flipped coordinates -
%will need to account for the negative value of the coordinates).

%LVamFocusOnRangesDnowDoAccel.m - now take consecutive velocities and
%compute acceleration. 

%LVamFocusOnRangeCveryGOODgetVelocity.m - This version collects points from
%the range, turns it into a corner/initcorner and computes the velocity.

%LVamFocusOnRangesB.m - This implements agents with different starting
%points on the trajectory that then expand their range to find next points.

%Rather than starting wide, might as well start with ALL input, and then
%grow search outward from one of those points.

%LVahMakeRangeMatricesForSearchAndMatch.m - Need to make range matrices
%that search large and then smaller ranges (search large and then narrow in
%for specific selection of three adjacent elements).  Then compute
%prediction matrix and evaluate accuracy of prediction matrix.  Try to
%program in as general a manner as possible.  Search is a mask template
%that is generated similar to a prediction matrix.
%Need to use some predictive subtraction, once one "agent" accounts for
%some aspect of the stimulus, need to block out the parts that are
%accounted for. For example, have a repeating texture on top of a smooth
%gradient. Should account for them separately.

%Need library of functions: 1. search, 2. match, 3. move, 4. compute
%inverse, 5. generate prediction, etc.  Then need an agent builder matrix
%that takes this library of functions and creates an agent - must choose
%functions in a meaningful order - do as linear rather than as a matrix?
%But some vector buffers would need to be accepped at multiple time points.
%Search should be large scope of search finding components. Wide field will
%give inaccurate estimate. Smaller search will be more accurate.
%Matching should involve phase range selection.

%%Can't just choose random points.  Need triplets.  Can't just choose random
%triplets, need adjacent triplets. Should search an X range to find
%adjacent x points to compute y dynamics.  INPUT vectors should be discete
%individual vectors that should come from a range.  Output vectors are specific.  
%SEARCH vector is a range.

%LVaeGetsTripletsEXCELLENT.m - this version extracts triplets of the x,y
%indices where space==1.  Note that to get ALL the trajectory elements,
%need to use 10*rand because rand has lots of overlap.

%LVacGravityAgain.m - Does nice plot of gravity trajectory and creates it
%as an array "space" that is plotted as an image.

%LVaaAttractors.m
%Idea that can try to create a higher level agent that looks at start and
%end points of lower level agents and appembles them in an intelligent
%manner (selecting only relevant sub-agents).  Are they indexed by relevant
%components?
%Idea is to have the indices allowing attractor spreading.

%Figure
close all
clear all
figure('Position',[10 400 1400 800]);
fig2D=gcf; axpvay2D=gca; colormap(gray); 

%Code for selecting subplots
numplots=3;
h=zeros(1,numplots);
for k=1:numplots %number of subplots
    h(k)=subplot(1,numplots,k);
end

%CREATE THE PATTERN (in a way that I think might be relevant to
%understanding)
TT=80;
pvay=zeros(3,TT); %vectors of locations pvay=pos,vel, acc, y dimension
Wvy=[1 1 1; 0 1 1; 0 0 1]; %Dynamical gravitation matrix
ground=2; %Level of ground (to avoid problems when it is zero or negative)
pvay(:,1)=[ground 20 -0.8]'; %Starting point pvay for trajectory pos=ground, vel=10, acc=-0.8
pvax=zeros(3,TT); %vectors of locations position, velocity, acceleration, x
Wvx=[1 1 0; 0 1 1; 0 0 1]; %Gravity matrix (row1,col2 is vel effect on pos
pvax(:,1)=[1 1 0]';  %Starting vector pos=1, vel=1, acc=0.

for tt=1:TT-1
    pvay(:,tt+1)=Wvy*pvay(:,tt); %Update x position 
    if pvay(1,tt+1)<ground; pvay(1,tt+1)=ground; pvay(2,tt+1)=-0.5*pvay(2,tt+1); end
    pvax(:,tt+1)=Wvx*pvax(:,tt);
end

%Create image separately so can compute max
maxspvy=ceil(max(pvay(1,:))); %Computes max pvay to determine plot range.
maxspvx=ceil(max(pvax(1,:))); %Computes max pvax value for plot range
space=zeros(maxspvy+10,maxspvx+10); %The overall space of the plot necessary for image pixelation.
for tt=1:TT
space(maxspvy+10-round(pvay(1,tt)),round(pvax(1,tt)))=1; %Set values in pixel image to 1 where
%trajectory is present.  eed to subtract from maxspvy to match plot
end

axes(h(1));
plot(pvax(1,:),pvay(1,:)); hold on;
xlim([1 maxspvx+10]); ylim([0 maxspvy+10]);

realspace=space; %Define realspace as the actual trajectory.
showspace=space; %Define show space as the space being shown in image.

axes(h(2));
image(64-64*0.2*showspace);

%%%%%%%%%%%%%%%%%%
%CREATE THE AGENTS TO UNDERSTAND THE PATTERN
%Use RANGE to find initial points, then narrow down to small range.

%%%%%%%
%PROPAGATE triplets and create agents that predict on basis of triplets.

%THIS FINDS THE INITIAL STARTING POINTS
[cx,cy]=find(space==1); %Find the spots where space==1
sense=[cx cy] %Create an column of x,y indices for those spots

AA=20; %length(sense) %Number of selected points.  length(sense) gives all
tripinx=ceil(rand(AA,1)*maxspvx) %Find the indexes for AA number of points

triplet=[sense(tripinx,:) ones(AA,1)]; %Create a triplet of values for each of the selcted points
tripadj=[maxspvy+10-triplet(:,1) triplet(:,2) ones(AA,1)] %%Tripajusted = Define the location of these points in image coordinates maxspvy-triplet

axes(h(1));
plot(tripadj(:,2),tripadj(:,1),'r*','MarkerSize',10);

WW=2; %Number of matrices in agent
VV=1; %Number of vectors in agent.
MM=10; %Number of steps of increasing range selection.
initrange=1.2;

axes(h(2));
%Make agents that expand search.
for mm=1:MM %Number of expansion search steps
    axes(h(2));
    
    for aa=1:AA; %Number of agents
        if mm==1;
            ag(aa).searchspace=zeros(maxspvy+10,maxspvx+10); %Create the plot of search space
            ag(aa).predspace=ag(aa).searchspace; %Create the plot of prediction space
            allpredspace=zeros(maxspvy+10,maxspvx+10);

        end
    %for ww=WW; %Number of agent matrix components that interact within agent
    %vector and Matrix 1 is expansion
        ag(aa).W(1).vvinit=triplet(aa,:); %Set the INITIAL STARTING POINT vvinit
        showspace(ag(aa).W(1).vvinit(1,1),ag(aa).W(1).vvinit(1,2))=0.2; %Mark vvinit in space
        range=1*mm; %initrange^mm; %Define Range of exploration that expands with internal agent time.
        ag(aa).W(1).Wrg=range*eye(4,4); %Range expansion matrix
        rx=maxspvx/50;
        ry=maxspvy/50; %Note that these will be used in opposite way x and y for image
        ag(aa).W(1).vvstrg=[ry -ry rx -rx]'; %vector for expansion of range x and y opposite for image
        ag(aa).W(1).vvrg=ag(aa).W(1).vvinit +  ag(aa).W(1).Wrg*ag(aa).W(1).vvstrg; %Search range
        
        for xx=1:maxspvy+10; %Needed to go larger to avoid losing points, not sure why
            for yy=1:maxspvx+10; %coordinates x and y are opposite due to image?
               if xx<ag(aa).W(1).vvrg(1,1) && xx>ag(aa).W(1).vvrg(2,1) && yy<ag(aa).W(1).vvrg(3,2) && yy>ag(aa).W(1).vvrg(4,2);
                   ag(aa).searchspace(xx,yy)=1;
                   showspace(xx,yy)=0.5;
               end
            end
        end
        
        %showspace
        
%         for zz=1:TT
%         space(maxspvy+10-ceil(pvay(1,zz)),ceil(pvax(1,zz)))=1;
%         end
        %showspace=showspace(1:maxspvy+10,1:maxspvx+10)+realspace;
        showspace(ag(aa).W(1).vvinit(1,1),ag(aa).W(1).vvinit(1,2))=0.2; %Repeat this so visible
        
        ag(aa).testspace=ag(aa).searchspace(1:maxspvy+10,1:maxspvx+10).*realspace;
        [tcx tcy]=find(ag(aa).testspace==1); %Find next spot within search space
        ag(aa).W(2).vvnext=[tcx tcy ones(length(tcy),1)]; %This is the next found position in the search space
        
        ag(aa).W(2).Wxva=zeros(3,3);
        ag(aa).W(2).Wdyny=zeros(3,3);
        %Wcorner=[1 0 0.5; 0 1 0.5; 0 0 1];
        %ag(aa).W(2).initcorner=Wcorner*ag(aa).W(1).vvinit';
        
        %IF THE number of points in test space are nonzero, then compute
        %corner to get velocity and acceleration AND then Wdyn
        if length(tcy)>0;
            for cc=1:length(tcy);
                %CREATE THE "corner" for computing affine matrix for
                %velocity and acceleration
            %ag(aa).W(2).corner(cc)=[ag(aa).W(1).vvinit'+[0.5 0 0]'  ag(aa).W(1).vvinit'  ag(aa).W(1).vvinit'+[0 0.5 0]'];
            ag(aa).W(2).corner(cc).c=[ag(aa).W(2).vvnext(cc,:)'+[0.5 0 0]'  ag(aa).W(2).vvnext(cc,:)'  ag(aa).W(2).vvnext(cc,:)'+[0 0.5 0]'];
            
                if cc>1  %Compute the velocity affine transform.
                ag(aa).W(2).Wvel(cc).c=ag(aa).W(2).corner(cc).c*inv(ag(aa).W(2).corner(cc-1).c);

                ag(aa).W(2).corner(cc-1).c   %Print the corners and the velocity affine matrix              
                ag(aa).W(2).corner(cc).c
                ag(aa).W(2).Wvel(cc).c
                ['above shows: corner(cc-1) corner(cc) and Wvel(cc)']
                end
                if cc>2  %Compute the acceleration affine transform

                    ag(aa).W(2).Wacc(cc).c=ag(aa).W(2).Wvel(cc).c*inv(ag(aa).W(2).Wvel(cc-1).c);
                    ag(aa).W(2).Wvel(cc-1).c %Print the velocity matrix and acceleration matrix
                    ag(aa).W(2).Wvel(cc).c
                    ag(aa).W(2).Wacc(cc).c
                    ['above: Wvel(cc-1) Wvel(cc) Wacc(cc)'] 
                
                    %if cc>2 && cc<6 - Create the arrays of pos, vel, acc.
                    %Wxva and Wyva
                    ag(aa).W(2).Wxva(1,cc-2)=ag(aa).W(2).corner(cc).c(1,2); %Take pos from corner
                    ag(aa).W(2).Wxva(2,cc-2)=ag(aa).W(2).Wvel(cc).c(1,3);
                    ag(aa).W(2).Wxva(3,cc-2)=ag(aa).W(2).Wacc(cc).c(1,3);
                    ag(aa).W(2).Wxva                     
                    ['above: ag(aa).W(2).Wxva']
                    ag(aa).W(2).Wyva(1,cc-2)=ag(aa).W(2).corner(cc).c(2,2); %Take pos from corner
                    ag(aa).W(2).Wyva(2,cc-2)=ag(aa).W(2).Wvel(cc).c(2,3);
                    ag(aa).W(2).Wyva(3,cc-2)=ag(aa).W(2).Wacc(cc).c(2,3); 
                    ag(aa).W(2).Wyva                     
                    ['above: ag(aa).W(2).Wyva']
                    
                    %Compute the dynamical matrix from the xva arrays. 
                    if cc>5; %ag(aa).W(2).Wxva(1,4)>0; %This tests that there are at least four columns active in Wxva
                        ag(aa).W(2).Wdynx=ag(aa).W(2).Wxva(:,2:4)*inv(ag(aa).W(2).Wxva(:,1:3));
                        ag(aa).W(2).Wdyny=ag(aa).W(2).Wyva(:,2:4)*inv(ag(aa).W(2).Wyva(:,1:3));
                        ag(aa).W(2).Wdynx
                        ag(aa).W(2).Wdyny
                      ['above: ag(aa).W(2).Wdynx and dyny'] 
                      PP=80;
                      ag(aa).W(2).vvpredx=zeros(3,PP);
                      ag(aa).W(2).vvpredy=zeros(3,PP);
                      for pp=1:PP
                        ag(aa).W(2).vvpredx(:,pp)=ag(aa).W(2).Wdynx^(pp)*ag(aa).W(2).Wxva(:,1);
                        ag(aa).W(2).vvpredy(:,pp)=[1 1 0; 0 1 0; 0 0 1]^(pp)*ag(aa).W(2).Wyva(:,1);
                        plotpredx=round(ag(aa).W(2).vvpredx(1,pp));
                        plotpredy=round(ag(aa).W(2).vvpredy(1,pp));
                        if plotpredx<maxspvy+10 && plotpredy<maxspvx+10 && plotpredx>0 && plotpredy>0;
                            ag(aa).predspace(plotpredx,plotpredy)=1.0; %[1:1:PP]
                        end
                      end
                      ag(aa).W(2).vvpredx
                      ag(aa).W(2).vvpredy
                      ['above is ag(aa).W(2).vvpredx']
                      image(64-64*(0.3*realspace + ag(aa).predspace(1:maxspvy+10,1:maxspvx+10))); %(1:maxspvy+10,1:maxspvx+10)));
                      %drawnow
                      aa
                      %pause(0.02);
                    end
                end %if cc>2    
            
            end %end for cc=1:length(tcy)
        end %if length(tcy>0)

                    ['length(tcy)   aa']   %Show the current length of the selected segment and agent number
                    [length(tcy)  aa]
                    %pause
        
        %ag(aa).W(2).Winv=
%         ag(aa).W(1).vvinit
%         ag(aa).W(1).vvrg
%         space
        image(64-64*(0.3*showspace)); %(1:maxspvy+10,1:maxspvx+10)));
        drawnow
%pause
    end %for aa
    for aa=1:AA
        allpredspace=allpredspace+ag(aa).predspace;
    end
    
    axes(h(3));
    image(64-64*(0.2*realspace+allpredspace));
    
    
end %for mm




%OLD SINE WAVE SEARCH PATTERN - Not really necessary
% %Start with a random output from "active" elements.
%  AA=3; %Number of agents
%  WW=1; %Number of matrices in each agent (search, predict, etc). Have a Wvy by Wvx array?
%  VV=3; %Number of vectors in each agent.  Start with input, search, output.
% % %Need RANGE matrices!!
% pp=200;
% ff=0.05; %This sets the frequency
% ampx=-.05; %This sets the amplitude
% ampy=0; %-.05
% LLx=0; %This changes line slope x
% LLy=0; %This changes line slope for y
% cx=26; %10;
% cy=13; %5;
% 
%  for aa=1:AA;
%      for ww=1:WW %Indicate with Ws, but Could have search matrix be 
%          %SET UP THE DYNAMICAL MATRIX
%         ag(aa).Ws(ww).W=zeros(6, 6);
%         ag(aa).Ws(ww).W(1:3,1:3)=[1 -ff 0; 1-ampx  1-ff  -cx; 0 0 1]; %; -0.1*cx 0 0*0.1*cx]; %0.9]; %; 0 0 0];
%         ag(aa).Ws(ww).W(4:6,4:6)=[1 1 0; 0 1 0; 0 0 1]; %[1 -ff 0; 1-ampy  1-ff  -cy; 0 0 1];
%         ag(aa).Ws(ww).vv=[cx+5 16 1 cy-3 .4 1]';
%      end
%  end
%  
%  axes(h(1));
%  for pp=1:pp
%    for aa=1:AA;
%      for ww=1:WW %Indicate with Ws, but Could have search matrix be 
%          ag(aa).Ws(ww).Vsx(:,pp+1)=(ag(aa).Ws(ww).W)^pp*ag(aa).Ws(ww).vv;; %Standard dynamical matrix.
%          plot(ag(aa).Ws(ww).Vsx(1,:),ag(aa).Ws(ww).Vsx(4,:),'g*','MarkerSize',14); %Search step matrix         
%      end
%    end
%  end
% 
