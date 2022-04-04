%Penentuan Episenter dengan metode inversi non-linear menggunakan
%Simulated Annealing Algorithm
%Mohammad Rheza Zamani

clear all;
clc;
%Earthquake station position
x = [20 50 40 10];
y = [10 25 50 40];
%Synthetic Model
xs = 40;
ys = 30;
model_sin = [xs ys];
vp = 7;
%Data sintetik
tobs = t_cal(xs,ys,x,y,vp);
%Definisi Ruang Model
nitr = 500; 
x_min = 1;
x_max = 50;
y_min = 1;
y_max = 50;
T = 5;
dec = 0.05;
%Model awal
%model1(1,1) = x_min + rand*(x_max-x_min);
%model1(1,2) = y_min + rand*(y_max-y_min);
model1(1,1) = 20;
model1(1,2) = 20;
tcal1(1,:) = t_cal(model1(1),model1(2),x,y,vp)
E1 = misfit(tobs(1,:),tcal1(1,:));

for itr = 1 : nitr
    model2(1,1) = x_min + rand*(x_max-x_min);
    model2(1,2) = y_min + rand*(y_max-y_min);
    tcal2(1,:) = t_cal(model2(1),model2(2),x,y,vp);
    E2 = misfit(tobs(1,:),tcal2(1,:));
    delta_E = E2 - E1;
    if delta_E<0
        model1 = model2;
        E1 = E2;
        if model1(1)<x_min && model1(2)<y_min
            model1 = [x_min ymin];
        elseif model1(1)>x_max && model1(2)>y_max
            model1 = [x_max y_max];
        end
    else
        P = exp(-delta_E/T); 
        if P>= rand
           model1 = model2;
           E1 = E2;
        if model1(1)<x_min && model1(2)<y_min
            model1 = [x_min ymin];
        elseif model1(1)>x_max && model1(2)>y_max
            model1 = [x_max y_max];
        end
    end
    end
   Egen(itr) = E1;
   T = T*(1-dec);
   Temperature(itr) = T;
%Ploting model
hold on
plot(model1(1),model1(2),'bo','MarkerFaceColor','r','MarkerSize',7)
plot(x,y,'rv','MarkerFaceColor','b','MarkerSize',10)
plot(xs,ys,'go','MarkerFaceColor','k','MarkerSize',7) 
xlim([x_min x_max*1.5])
ylim([y_min y_max*1.5])
grid on
xlabel('x (m)','FontWeight','bold')
ylabel('y (m)','FontWeight','bold')
title(['Penentuan Episenter || Hasil Inversi ke-.',num2str(itr)],'FontWeight','bold')
subtitle(['x = ',num2str(model1(1)),' ; y = ',num2str(model1(2)),' || ERMS = ',num2str(Egen(itr))],'FontWeight','bold')
legend('Episenter Hasil inversi','Stasiun','Episenter Teori')
set(gcf, 'Position', get(0, 'Screensize'));
pause(0.1)
if itr ~= nitr
clf
end   
end

%Plot grafik misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
grid on