%Perhitungan waktu
function tcal = t_cal(x1,y1,x0,y0,vp)
lx = length(x0);
for i = 1 : lx
    tcal(i) = sqrt((x1-x0(i)).^2+(y1-y0(i)).^2)./vp;
end
end