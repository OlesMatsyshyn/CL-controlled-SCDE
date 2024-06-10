newcolors = [60 180 75
             0 130 200 
             245 130 48
             145 30 180
             240 50 230
             210 245 60
             230 25 75
             0 128 128
             220 190 255
             170 110 40
             255 250 200
             128 0 0
             170 255 195
             128 128 0
             255 215 180
             0 0 128
             128 128 128
             255 255 255
             0 0 0]/255;
colororder(newcolors)

DAarr = zeros(5,100);
DBarr = zeros(5,100);
J_arr = zeros(5,100);
eta = zeros(1,6);
DE = zeros(1,6); 
for i = 0 : 5
    if mod(i,2) == 0
        index = num2str(fix(i/2));
    else
        index = [num2str(fix(i/2)) '5'];
    end
    name = ['SCDEres00' index '.mat'];
    load(name)
    DAarr(i+1,:) = DAq;
    DBarr(i+1,:) = DBq;
    J_arr(i+1,:) = Jxq;
    Jxq = real(Jxq);
    eta(i+1) = abs(abs(max(Jxq))-abs(min(Jxq)))/(abs(max(Jxq))+abs(min(Jxq)));
    DE(i+1) = DeltaE;
end
 
%%

hold on
grid off
x0=300;
y0=150;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'TickLength',[0.015 0.015]);
subplot(2,2,1)
hold on
grid off
for i = 1 :2: 5
    plot(qx,DAarr(i,:),'Linewidth',3,'Color',newcolors(i,:));
end
leg = legend('$\Delta_E/[\rm eV] = 0$','$\Delta_E/[\rm eV] = 5\cdot 10^{-3} $','$\Delta_E/[\rm eV] = 10\cdot 10^{-3}$');
set(leg,'Interpreter','latex','location','northwest');
set(leg,'FontSize',20);
ax = gca;
set(gca,'XAxisLocation', 'bottom', 'YAxisLocation', 'left');

yrule = ax.YAxis;
xrule = ax.XAxis;
yrule.FontSize = 20;
xrule.FontSize = 20;

tit = title('$\Delta_A$');
set(tit,'Interpreter','latex');
set(tit,'FontSize',40);

yticks((0:0.4:1.2)/100)
yticklabels(string(((0:0.4:1.2)/100)))
xticks((-0.5:0.25:0.5))
xticklabels(string((-0.5:0.25:0.5)))
ylim([0 1.3/100])

ylabel('[eV]','interpreter','latex','Fontsize',30)
xlabel('$q_x a$','interpreter','latex','Fontsize',30)

%%
subplot(2,2,2)
hold on
grid off
x0=300;
y0=150;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'TickLength',[0.015 0.015]);

hold on
grid off
for i = 1 :2: 5
    plot(qx,DBarr(i,:),'Linewidth',3,'Color',newcolors(i,:));
end
leg = legend('$\Delta_E/[\rm eV] = 0$','$\Delta_E/[\rm eV] = 5\cdot 10^{-3} $','$\Delta_E/[\rm eV] = 10\cdot 10^{-3}$');
set(leg,'Interpreter','latex','location','northwest');
set(leg,'FontSize',20);
ax = gca;
set(gca,'XAxisLocation', 'bottom', 'YAxisLocation', 'left');

yrule = ax.YAxis;
xrule = ax.XAxis;
yrule.FontSize = 20;
xrule.FontSize = 20;

tit = title('$\Delta_B$');
set(tit,'Interpreter','latex');
set(tit,'FontSize',40);

yticks((0:0.4:1.2)/100)
yticklabels(string(((0:0.4:1.2)/100)))
xticks((-0.5:0.25:0.5))
xticklabels(string((-0.5:0.25:0.5)))
ylim([0 1.3/100])

ylabel('[eV]','interpreter','latex','Fontsize',30)
xlabel('$q_x a$','interpreter','latex','Fontsize',30)
%%
subplot(2,2,3)
hold on
grid off
x0=300;
y0=150;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'TickLength',[0.015 0.015]);

hold on
grid off
for i = 1 :2: 5
    plot(qx,real(J_arr(i,:)),'Linewidth',3,'Color',newcolors(i,:));
end
leg = legend('$\Delta_E/[\rm eV] = 0$','$\Delta_E/[\rm eV] = 5\cdot 10^{-3} $','$\Delta_E/[\rm eV] = 10\cdot 10^{-3}$');
set(leg,'Interpreter','latex','location','northwest');
set(leg,'FontSize',20);
ax = gca;
set(gca,'XAxisLocation', 'bottom', 'YAxisLocation', 'left');

yrule = ax.YAxis;
xrule = ax.XAxis;
yrule.FontSize = 20;
xrule.FontSize = 20;

tit = title('$J_x$');
set(tit,'Interpreter','latex');
set(tit,'FontSize',40);

yticks((-2.5:1:2.5)/1000)
yticklabels(string(((-2.5:1:2.5)/1000)))
xticks((-0.5:0.25:0.5))
xticklabels(string((-0.5:0.25:0.5)))
ylim([-2.5 2.5]/1000)

ylabel('[eV]','interpreter','latex','Fontsize',30)
xlabel('$q_x a$','interpreter','latex','Fontsize',30)

%%
subplot(2,2,4)
grid off
x0=300;
y0=150;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'TickLength',[0.015 0.015]);
plot(DE, eta,'ko','Linewidth',3);
ylabel('$\eta$','interpreter','latex','Fontsize',30)
xlabel('$\Delta_E$','interpreter','latex','Fontsize',30)

ax = gca;
set(gca,'XAxisLocation', 'bottom', 'YAxisLocation', 'left');

yrule = ax.YAxis;
xrule = ax.XAxis;
yrule.FontSize = 20;
xrule.FontSize = 20;

yticks((0:0.05:0.2))
yticklabels(string(((0:0.05:0.2))))
xticks(((0:0.4:1.2)/100))
xticklabels(string(((0:0.4:1.2)/100)))
ylim([0 0.2])