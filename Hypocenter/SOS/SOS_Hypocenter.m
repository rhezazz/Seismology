%Penentuan hiposenter dengan metode inversi non-linear menggunakan
%algoritma mSOS (modified Symbiotic Organism Search)
%Mohammad Rheza Zamani
clear all;
clc;
%Earthquake station position
x = [20 50 40 10];
y = [10 25 50 40];
z = [0.5 0.5 0.5 0.5];
%Synthetic Model
xs = 40;
ys = 30;
zs = -20;
vp = 7;
%Data sintetik
tobs = t_cal(xs,ys,zs,x,y,z,vp);
%Definisi Ruang Model
npop = 100; 
nitr = 100; 
x_min = -100;
x_max = 100;
y_min = -100;
y_max = 100;
z_min = -100;
z_max = 100;

%Membuat Video
v = VideoWriter('Penentuan Hiposenter SOS.avi');
open(v);
%Membuat model acak
for ipop = 1 : npop
    model(ipop,1) = x_min + rand*(x_max-x_min);
    model(ipop,2) = y_min + rand*(y_max-y_min);
    model(ipop,3) = z_min + rand*(z_max-z_min);
    tcal(ipop,:) = t_cal(model(ipop,1),model(ipop,2),model(ipop,3),x,y,z,vp);
    err(ipop) = misfit(tobs(1,:),tcal(ipop,:));
end
%Proses inversi
for itr = 1 : nitr
    for i = 1 : npop
        idx = find(err ==min(err));
        model_best = model(idx(1),:);
        %Mutualisme
        j = randi(npop,1);
        k = randi(npop,1);
        if j==i || k==i
            j = randi(npop,1);
            k = randi(npop,1);
        end
        model_mut = [model(i,:);model(j,:)];
        mv_m =(model(i,:)+model(j,:))/2;
        bf = 1;
        for l = 1 : 2
            mod_mut(l,:) = model_mut(l,:) + rand*(model(k)-mv_m*bf);
            if mod_mut(l,1) < x_min
                mod_mut(l,1) = x_min;
            end
            if mod_mut(l,2) < y_min
                mod_mut(l,2) = y_min;
            end
            if mod_mut(l,3) < z_min
                mod_mut(l,3) = z_min;
            end
            if mod_mut(l,1) > x_max
                mod_mut(l,1) = x_max;
            end
            if mod_mut(l,2) > y_max
                mod_mut(l,2) = y_max;
            end
            if mod_mut(l,3) > z_max
                mod_mut(l,3) = z_max;
            end
        end
        %Hitung model untuk prosedur mutualisme
        for l = 1 : 2
            [tcal_mut] = t_cal(mod_mut(l,1),mod_mut(l,2),mod_mut(l,3),x,y,z,vp);
            err_mut(l) = misfit(tobs,tcal_mut);
            %Update model jika  nilai misfit lebih baik proses mutualisme
            if l == 1
                if err_mut(l)<err(i)
                    model(i,:) = mod_mut(l,:);
                    err(i) = err_mut(l);
                end
            else
                if err_mut(l)<err(j)
                    model(j,:) = mod_mut(l,:);
                    err(j) = err_mut(l);
                end
            end
        end
        %Komensalisme
        j = randi(npop,1);
        if j == i
            j = randi(npop,1);
        end
        mod_com = model(i) +(0.4+0.9*rand)*(model_best-model(j));
             if mod_com(1) < x_min
                mod_com(1) = x_min;
            end
            if mod_com(2) < y_min
                mod_com(2) = y_min;
            end
            if mod_com(3) < z_min
                mod_com(3) = z_min;
            end
            if mod_com(1) > x_max
                mod_com(1) = x_max;
            end
            if mod_com(2) > y_max
                mod_com(2) = y_max;
            end
            if mod_com(3) > z_max
                mod_com(3) = z_max;
            end
        %Perhitungan misfit untuk prosedur komensalisme
        [tcal_com] = t_cal(mod_com(1),mod_com(2),mod_com(3),x,y,z,vp);
         err_com = misfit(tobs,tcal_com);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if err_com < err(i)
             model(i,:) = mod_com;
             err(i) = err_com;
         end
         %Parasitisme
         j = randi(npop,1);
         if j == i 
             j = randi(npop,1);
        end
         mod_par = model(i,:);
         p1 = randi(3,1)
         if p1 ==1
            mod_par(i,1) = x_min + rand*(x_max-x_min);
         elseif p1==2
            mod_par(i,2) = y_min + rand*(y_max-y_min);
         elseif p1==3
            mod_par(i,3) = z_min + rand*(z_max-z_min);
         end
         %Perhitungan misfit untuk tahap parasitisme
         [tcal_par] = t_cal(mod_par(1),mod_par(2),mod_par(3),x,y,z,vp);
         err_par = misfit(tobs,tcal_par);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if err_par < err(i)
             model(j,:) = mod_par(1,:);
             err(j) = err_par;
         end
    end
    %Update model terbaik untuk setiap iterasi
    Emin = 100;
    for ipop = 1 : npop
        Emin = err(ipop);
        model_baru = model(ipop,:);
    end
    %Nilai misfit terbaik
    Egen(itr)=Emin;
%Ploting model
hold on
plot3(model_baru(1),model_baru(2),-model_baru(3),'bo','MarkerFaceColor','r','MarkerSize',7)
plot3(x,y,z,'rv','MarkerFaceColor','b','MarkerSize',10)
plot3(xs,ys,zs,'go','MarkerFaceColor','k','MarkerSize',7) 
xlim([x_min x_max])
ylim([y_min y_max])
zlim([z_min z_max])
grid on
xlabel('x (m)','FontWeight','bold')
ylabel('y (m)','FontWeight','bold')
zlabel('z (m)','FontWeight','bold')
title(['Penentuan Hiposenter || Hasil Inversi ke-.',num2str(itr)],'FontWeight','bold')
subtitle(['x = ',num2str(model_baru(1)),' ; y = ',num2str(model_baru(2)),' ; z = ',num2str(-model_baru(3)),' || ERMS = ',num2str(Egen(itr))],'FontWeight','bold')
legend('Hiposenter Hasil inversi','Stasiun','Hiposenter Teori')
set(gcf, 'Position', get(0, 'Screensize'));
view(3)
pause(0.05)
frame = getframe(gcf);
writeVideo(v,frame);
if itr ~= nitr
clf
end
end
close(v);
%Plot grafik misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on

saveas(figure(2),'Grafik Misfit Hiposenter SOS')