function drawMito(pos,length,height,rotation)

x0 = pos(1);
y0 = pos(2);
rotationVector = ([cos(rotation), -sin(rotation); sin(rotation), cos(rotation)]);
xExtrema = ([-length/2, length/2, length/2, -length/2]);
yExtrema = ([-height/2, -height/2, height/2, height/2]);

for i = 1:4
    rotateExtrema(:,i) = rotationVector * [xExtrema(i); yExtrema(i)];
end

x_lower_left = x0 + rotateExtrema(1,1);
x_lower_right = x0 + rotateExtrema(1,2);
x_upper_right = x0 + rotateExtrema(1,3);
x_upper_left = x0 + rotateExtrema(1,4);
y_lower_left = y0 + rotateExtrema(2,1);
y_lower_right = y0 + rotateExtrema(2,2);
y_upper_right = y0 + rotateExtrema(2,3);
y_upper_left = y0 + rotateExtrema(2,4);

x_coor = [x_lower_left x_lower_right x_upper_right x_upper_left];
y_coor = [y_lower_left y_lower_right y_upper_right y_upper_left];

patch('Vertices',[x_coor; y_coor]','Faces',[1 2 3 4],'Edgecolor','k','Facecolor','k');
end