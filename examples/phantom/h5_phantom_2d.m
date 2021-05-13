% STEP1: create a 2d phantom pattern (3 circles)
phantom_name = 'circles.h5';
% resolution 
dx = 1; % mm 
dy = 1; % mm

x = -100:dx:100; % x range
y = -80:dy:80; % y range
[posx posy] = meshgrid(x,y);
tissue_distrib = zeros(size(posx)); % background value 0.
tissue_distrib((posx.^2 + posy.^2)<25^2)= 1;   % radius 25, center at the origin, tissue index 1.
tissue_distrib(((posx - 50).^2+(posy - 50).^2)<15^2)= 2;   % radius 15, center at (50,50), tissue index 2.
tissue_distrib(((posx + 50).^2+(posy + 50).^2)<15^2)= 3;   % radius 15, center at (-50,-50),tissue index 3.
% plot the phantom
imagesc(tissue_distrib)
%---------------------------------------------------------------------
% STEP2: create the h5 phantom file.
if exist(phantom_name,'file')==2
    delete(phantom_name);
end

res=[dx dy]*1e-3; % unit in m.
h5create(phantom_name, '/phantom/resolution', size(res));
h5write(phantom_name, '/phantom/resolution', res);

tissue_distrib = uint16(tissue_distrib);
h5create(phantom_name, '/phantom/tissue_dist', size(tissue_distrib), 'Datatype', 'uint16');
h5write(phantom_name, '/phantom/tissue_dist', tissue_distrib);