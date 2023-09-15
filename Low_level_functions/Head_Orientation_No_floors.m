function maze_view=Head_Orientation_No_floors(TrackingData)
 % Head_Orientation_No_floors determines the orientation of an animal in a chamber.
    %
    % Input:
    % - TrackingData: Structure containing the position and orientation data
    %
    % Output:
    % - maze_view: Matrix representing the coordinates in x,y,z of the head orientation projection on the maze walls

%Linear projection of rotated vectors (quiver) corresponding to orientation
%of the animal inside the chamber,
heaty=[];
heatx=[];

TrackingData.XPosition(TrackingData.XPosition>0.6)=nan;
TrackingData.YPosition(TrackingData.YPosition<0)=nan;
TrackingData.XPosition(TrackingData.XPosition<0)=nan; 
TrackingData.YPosition(TrackingData.YPosition>1.14)=nan;
TrackingData.ZPosition(TrackingData.ZPosition<0.2)=nan;
TrackingData.ZPosition(TrackingData.ZPosition>1.4)=nan;

%floor separation
floor=struct;
floor.A=TrackingData.ZPosition>0.2 & TrackingData.ZPosition<0.6;%first floor half rigth
floor.B=TrackingData.ZPosition>0.6 & TrackingData.ZPosition<1;%second floor half rigth
floor.C=TrackingData.ZPosition>1 ;%third floor half left


index=1:length(TrackingData.XPosition);

    WQuat= TrackingData.WQuat;
    iXQuat=TrackingData.iXQuat;
    jYQuat=TrackingData.jYQuat;
    kZQuat=TrackingData.kZQuat;
    
        R=quat2rotm([WQuat iXQuat jYQuat kZQuat]);%q = [w x y z]
   
    ux = reshape(R(1,1,:),[length(TrackingData.XPosition) 1]);
    vx = reshape(R(2,1,:),[length(TrackingData.XPosition) 1]);
    wx = reshape(R(3,1,:),[length(TrackingData.XPosition) 1]);
    uy = reshape(R(1,2,:),[length(TrackingData.XPosition) 1]);
    vy = reshape(R(2,2,:),[length(TrackingData.XPosition) 1]);
    wy = reshape(R(3,2,:),[length(TrackingData.XPosition) 1]);
    uz = reshape(R(1,3,:),[length(TrackingData.XPosition) 1]);
    vz = reshape(R(2,3,:),[length(TrackingData.XPosition) 1]);
    wz = reshape(R(3,3,:),[length(TrackingData.XPosition) 1]); 


clear WQuat iXQuat jYQuat kZQuat R i
%index=find(TrackingData.XPosition(social.A | social.B | social.C | social.D)); 

coefficients=[];
coefficientsz=[];

coefficients(:,1)=vy./uy; %m - slope 
coefficients(:,2)=TrackingData.YPosition(index)-(TrackingData.XPosition(index).*coefficients(:,1));% b point
x=[];
y=[];
x1=[];
x2=[];
x(:,1)=coefficients(:,2)>0&coefficients(:,2)<1.145; %evaluating y at x=0 (left wall)
x(:,2)=(0.6.*coefficients(:,1))+coefficients(:,2)>0&(0.6.*coefficients(:,1))+coefficients(:,2)<1.145; %evaluating y at x=0.6 (right wall)
x=x(:,1)==1 & x(:,2)==1;
x=double(x);
x(x==0)=nan;
x(h.UData>0&x>0)=0.6;
x(h.UData<0&x>0)=0;
y=(0.6.*coefficients(:,1))+coefficients(:,2);
y2=coefficients(:,2);
y(x==0)=y2(x==0);
x1=(-coefficients(:,2)./coefficients(:,1))>0&(-coefficients(:,2)./coefficients(:,1))<0.6; %evaluating x at y=0 (south wall)
x2=((1.145-coefficients(:,2))./coefficients(:,1))>0&((1.145-coefficients(:,2))./coefficients(:,1))<0.6;%evaluating x at y=1.145 (north wall)
x3=x1==1 & x2==1;
x1=(-coefficients(:,2)./coefficients(:,1));
x2=((1.145-coefficients(:,2))./coefficients(:,1));
x(h.VData<0&x3)=x1(h.VData<0&x3);
x(h.VData>0&x3)=x2(h.VData>0&x3);
y(h.VData<0&x3)=0;
y(h.VData>0&x3)=1.145;

y1=coefficients(:,2); %y values at x=0
y2=(0.6.*coefficients(:,1))+coefficients(:,2); %y values at x=0.6
y(isnan(x))=nan;

%vectors that cross two orthogonal walls L (diagonal)
x(isnan(x)&x2>0&x2<0.6&h.VData>0)=x2(isnan(x)&x2>0&x2<0.6&h.VData>0);
y(isnan(y)&x2>0&x2<0.6&h.VData>0)=1.145;


x(isnan(x)&x1>0&x1<0.6&h.VData<0)=x1(isnan(x)&x1>0&x1<0.6&h.VData<0);
y(isnan(y)&x1>0&x1<0.6&h.VData<0)=0;

y(isnan(y)&y2>0&y2<1.145&h.UData>0)=y2(isnan(x)&y2>0&y2<1.145&h.UData>0);
x(isnan(x)&y2>0&y2<1.145&h.UData>0)=0.6;

y(isnan(y)&y1>0&y1<1.145&h.UData<0)=y1(isnan(x)&y1>0&y1<1.145&h.UData<0);
x(isnan(x)&y1>0&y1<1.145&h.UData<0)=0;



heaty=y;
heatx=x;

% Z axis slope component

% x=h.UData;
% y=h.VData;
% z=h.WData;
coefficientsz=[];
heatz=[];
coefficientsz(:,1)=(h.WData)./(h.UData); %m - slope z/x
coefficientsz(:,2)=TrackingData.ZPosition(index)-(TrackingData.XPosition(index).*coefficientsz(:,1));% b point
coefficientsz(:,3)=(h.WData)./(h.VData); %m - slope z/y
coefficientsz(:,4)=TrackingData.ZPosition(index)-(TrackingData.YPosition(index).*coefficientsz(:,3));% b point
z1=coefficientsz(:,2);
z2=(coefficientsz(:,1).*0.6)+coefficientsz(:,2);
z3=coefficientsz(:,4);
z4=(coefficientsz(:,3).*1.145)+coefficientsz(:,4);
ztop=(z4>1.4&heaty==1.145) | (z3>1.4&heaty==0) | (z2>1.4&heatx==0.6) | (z1>1.4&heatx==0);
zbott=(z4<0&heaty==1.145) | (z3<0&heaty==0) | (z2<0&heatx==0.6) | (z1<0&heatx==0);
heatz(heatx==0.6&~zbott&~ztop)=z2(heatx==0.6&~zbott&~ztop);
heatz(heatx==0&~zbott&~ztop)=z1(heatx==0&~zbott&~ztop);
heatz(heaty==1.145&~zbott&~ztop)=z4(heaty==1.145&~zbott&~ztop);
heatz(heaty==0&~zbott&~ztop)=z3(heaty==0&~zbott&~ztop);
heatz(heatz==0)=nan;
heatz(ztop)=1.4;
heatz(zbott)=0;
%Projection to top and bottom
x1=(-coefficientsz(:,2)./coefficientsz(:,1));
x2=((1.4-coefficientsz(:,2))./coefficientsz(:,1));
heatx((heatx==0 | heatx==0.6)&(zbott))=x1((heatx==0 | heatx==0.6)&(zbott));
heatx((heatx==0 | heatx==0.6)&(ztop))=x2((heatx==0 | heatx==0.6)&(ztop));

y1=(-coefficientsz(:,4)./coefficientsz(:,3));
y2=((1.4-coefficientsz(:,4))./coefficientsz(:,3));
heaty((heaty==0 | heaty==1.145)&(zbott))=y1((heaty==0 | heaty==1.145)&(zbott));
heaty((heaty==0 | heaty==1.145)&(ztop))=y2((heaty==0 | heaty==1.145)&(ztop));
heatz=heatz';
maze_view=[];
maze_view(:,2)=heaty;
maze_view(:,1)=heatx;
maze_view=maze_view(1:length(heatz),:);
maze_view(:,3)=heatz;
maze_view(end+1:length(heaty),:)=nan;
% close all
end
