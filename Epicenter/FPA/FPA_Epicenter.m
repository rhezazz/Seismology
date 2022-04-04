%Penentuan Episenter dengan metode inversi non-linear menggunakan
%algoritma (Flower Pollination Algorithm)
%Mohammad Rheza Zamani
%Reference : Yang, Xin-She. (2014). Nature-Inspired Optimization Algorithms. 
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
npop = 100; 
nitr = 500; 
x_min = 1;
x_max = 100;
y_min = 1;
y_max = 100;
prp = 0.8;

%Membuat Video
v = VideoWriter('Penentuan Episenter FA.avi');
open(v);

%Membuat model acak
for ipop = 1 : npop
    model(ipop,1) = x_min + rand*(x_max-x_min);
    model(ipop,2) = y_min + rand*(y_max-y_min);
    tcal(ipop,:) = t_cal(model(ipop,1),model(ipop,2),x,y,vp);
    err(ipop) = misfit(tobs(1,:),tcal(ipop,:));
end
 idx = find(err ==min(err));
model_best = model(idx(1),:);
%Proses Inversi
for iitr =  1 : nitr
    for i = 1 : npop
        if rand < prp
            n = 1;
            m = length(model_sin);
            beta = 1.5;
            [L]=levy(n,m,beta);%levy distribution
            npop_mod=model(i)+(0.1.*L.*(model(i)-model_best));
            if npop_mod(1) < x_min
                npop_mod(1) = x_min;
            end
            if npop_mod(2) < y_min
                npop_mod(2) = y_min;
            end
            if npop_mod(1) > x_max
                npop_mod(1) = x_max;
            end
            if npop_mod(2) > y_max
                npop_mod(2) = y_max;
            end
        else
            epsilon=rand;
            JK=randperm(npop);
            npop_mod=model(i)+epsilon*(model(JK(1),:)-model(JK(2),:));
            if npop_mod(1) < x_min
                npop_mod(1) = x_min;
            end
            if npop_mod(2) < y_min
                npop_mod(2) = y_min;
            end
            if npop_mod(1) > x_max
                npop_mod(1) = x_max;
            end
            if npop_mod(2) > y_max
                npop_mod(2) = y_max;
            end
        end
       model_new = npop_mod;
       tcal_new = t_cal(model_new(1),model_new(2),x,y,vp);
       E = misfit(tobs(1,:),tcal_new);
       %Update model
        if E < err(i)
             model(i,:) = model_new(1,:);
             err(i) = E;
         end
      end
          %Update model terbaik untuk setiap iterasi
    Emin = 100;
    for ipop = 1 : npop
        Emin = err(ipop);
        model_best = model(ipop,:);
    end
    %Nilai misfit terbaik
    Egen(iitr)=Emin;
    
 %Ploting model
hold on
plot(model_best(1),model_best(2),'bo','MarkerFaceColor','r','MarkerSize',7)
plot(x,y,'rv','MarkerFaceColor','b','MarkerSize',10)
plot(xs,ys,'go','MarkerFaceColor','k','MarkerSize',7) 
xlim([x_min x_max])
ylim([y_min y_max])
grid on
xlabel('x (m)','FontWeight','bold')
ylabel('y (m)','FontWeight','bold')
title(['Penentuan Episenter || Hasil Inversi ke-.',num2str(iitr)],'FontWeight','bold')
subtitle(['x = ',num2str(model_best(1)),' ; y = ',num2str(model_best(2)),' || ERMS = ',num2str(Egen(iitr))],'FontWeight','bold')
legend('Episenter Hasil inversi','Stasiun','Episenter Teori')
set(gcf, 'Position', get(0, 'Screensize'));
pause(0.05)
frame = getframe(gcf);
writeVideo(v,frame);
if iitr ~= nitr
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


%Levy Function
function [z] = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);

    sigma_u = (num/den)^(1/beta);

    u = normrnd(0,sigma_u^2,n,m); 
    
    v = normrnd(0,1,n,m);

    z = u./(abs(v).^(1/beta));
end