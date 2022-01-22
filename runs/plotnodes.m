close all
clear all

%Low order grid
coords_lo=load('COORDS_LO.dat');
x=coords_lo(:,1);
y=coords_lo(:,2);
z=coords_lo(:,3);
ip=coords_lo(:,4)+1;


%Edges nodes
coords_ho=load('COORDS_HO_edges.dat');
xho_e=coords_ho(:,1);
yho_e=coords_ho(:,2);
zho_e=coords_ho(:,3);
ip_e=coords_ho(:,4)+1;
clear coords_ho;

%Faces nodes
coords_ho=load('COORDS_HO_faces.dat');
xho_f=coords_ho(:,1);
yho_f=coords_ho(:,2);
zho_f=coords_ho(:,3);
ip_f=coords_ho(:,4)+1;
clear coords_ho;

%Volume nodes
coords_ho=load('COORDS_HO_vol.dat');
xho_v=coords_ho(:,1);
yho_v=coords_ho(:,2);
zho_v=coords_ho(:,3);
ip_v=coords_ho(:,4)+1;
clear coords_ho;


dx    = 0.15; dy = 0.15; dz = 0.01; % displacement so the text

for i=1:length(x)
    
    scatter3(x(i), y(i), z(i), 140, 'Filled'); hold on
    
    ipstr   = ip(i);
    
    dx    = 0.015*max(x); dy = 0.015*max(y); 0.001*max(z); % displacement so the text
    tx = text(x(i)+dx, y(i)+dy, z(i)+dz, int2str(ipstr), 'FontSize', 12);
end
hold on
for i=1:length(xho_e)
    
    scatter3(xho_e(i), yho_e(i), zho_e(i), 120);
    
    ipstr   = ip_e(i);
    dx    = 0.015*max(xho_e); dy = 0.015*max(yho_e); 0.001*max(zho_e); % displacement so the text
    tx = text(xho_e(i)+dx, yho_e(i)+dy, zho_e(i)+dz, int2str(ipstr), 'FontSize', 12);
end
for i=1:length(xho_f)
    
    scatter3(xho_f(i), yho_f(i), zho_f(i), 120);
    
    ipstr   = ip_f(i);
    dx    = 0.015*max(xho_f); dy = 0.015*max(yho_f); 0.001*max(zho_f); % displacement so the text
    tx = text(xho_f(i)+dx, yho_f(i)+dy, zho_f(i)+dz, int2str(ipstr), 'FontSize', 12);
end

for i=1:length(xho_v)
    
    scatter3(xho_v(i), yho_v(i), zho_v(i), 120);
    
    ipstr   = ip_v(i);
    dx    = 0.015*max(xho_v); dy = 0.015*max(yho_v); 0.001*max(zho_v); % displacement so the text
    tx = text(xho_v(i)+dx, yho_v(i)+dy, zho_v(i)+dz, int2str(ipstr), 'FontSize', 12);
end

