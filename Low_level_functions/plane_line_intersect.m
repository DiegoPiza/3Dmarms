function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.
 
% Edited by Diego B Piza Jan of 2020
%           Edits: works for array of multiple points, all the vectors need to
%           be the same length.
%          

if length(n)>3
    dim=2;
else
    dim=1;
end

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D =dot(n,u,dim);
N = -dot(n,w,dim);
check(1:length(D))=0;
 check(abs(D) < 10^-7 & N == 0)=2; % The segment is parallel to plane,% The segment lies in plane
 check(abs(D) < 10^-7 & ~(N == 0))=0;
%compute the intersection parameter
sI = N ./ D;
I = P0+ sI.*u;
I(abs(I)<1e-6)=0;
I=round(I,4,'significant');
check((sI < 0 | sI > 1))=3;%The intersection point  lies outside the segment, so there is no intersection
check(~(sI < 0 | sI > 1))=1;
check(isnan(sI))=nan;
check=check';
end
