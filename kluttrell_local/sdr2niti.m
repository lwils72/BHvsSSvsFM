function [n,t]=sdr2niti(strike,dip,rake)
% [n,t]=sdr2niti(strike,dip,rake)
%   input:
%     strike,dip,rake = 3 1-column vectors in Aki and Richards format
%   output:
%     n = normal vector components (3-column vector: nx,ny,nz)
%     t = traction vector components (3-column vector: tx,ty,tz)
%

%clear

  %strike=200;
  %dip=65;
  %rake=-120;
  
  nx=cosd(strike).*sind(dip);
  ny=-sind(strike).*sind(dip);
  nz=cosd(dip);
  
  tx=-sind(rake).*cosd(dip).*cosd(strike)+sind(strike).*cosd(rake);
  ty= sind(rake).*cosd(dip).*sind(strike)+cosd(strike).*cosd(rake);
  tz= sind(rake).*sind(dip);
  
  sx=sind(strike);
  sy=cosd(strike);
  sz=0;
  
  n=[nx,ny,nz];
  t=[tx,ty,tz];
  
  %figure(1),clf
  %plot3([0 nx],[0 ny],[0 nz],'.-'),hold on
  %plot3([0 sx],[0 sy],[0 sz],'.-k')
  %plot3([0 tx],[0 ty],[0 tz],'*-k')
  %plot3([-1,1],[0,0],[0,0],'k')
  %plot3([0,0],[-1,1],[0,0],'k')
  %plot3([0,0],[0,0],[-1,1],'k')
  %grid on,axis equal,xlabel('x'),ylabel('y'),zlabel('z'),view(20,30)
  %title(['\theta=',num2str(strike),', \delta=',num2str(dip),', \nu=',num2str(rake)])
    
  %[X,Y]=meshgrid([-1:.1:1],[-1:.1:1]);
  %Z=-(nx*X+ny*Y)/nz;
  %Z(find(abs(Z)>2))=NaN;
  %surface(X,Y,Z)
