 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
%% Band structure computation
SimpleBZ_con

sizeT = 20;
Tarr = linspace(0.000000001,0.02,sizeT);

DAq = Tarr*0+0.01; % seed value
DBq = Tarr*0+0.008; % seed value
Jq  = Tarr*0;
dkxdky = (1/nx) * (1/ny) * abs(a1_star_strained(1)  * a2_star_strained(2)...
                              -a1_star_strained(2)  * a2_star_strained(1))/(2*pi)^2;
if dkxdky == 0
    cprintf('err','The BZ is crazy, fix the strain...\n');
    return
end

tic
cprintf('hyper','Starting the canlculation... Hold on... It"s not so simple...\n')% XX component (just a nice dislay of the statement)
for it = 1:sizeT
    Temp = Tarr(it);% The system temperature
    disp(['Temperature array index ', num2str(it), ' started'])
    fFD  = @(x) 1/2-tanh(0.5*x/Temp)/2;
    iter = 0;             % number of iterations from this q completed

    discrA = 1;           % just to inititalize while
    discrB = 1;

    energy = zeros (4, nx+1, ny+1);
    DAnn   = zeros (4, nx+1, ny+1);
    DBnn   = zeros (4, nx+1, ny+1);
    while (discrA > 10^(-2) && abs(DAq(it)) > 10^(-5))||(discrB > 10^(-2) && abs(DBq(it)) > 10^(-5))  % convergence cryterium
        iter = iter + 1;  % update iteration index
        for i = 0 : nx
            for j = 0 : ny
                k = (i/nx) * a1_star_strained + (j/ny) * a2_star_strained;
                [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),0,DAq(it),DBq(it)));
                energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
                for n = 1 : nbands
                    DAnn(n,i+1,j+1)  =  Vec(:,n)'*dDAH_BdG*Vec(:,n);
                    DBnn(n,i+1,j+1)  =  Vec(:,n)'*dDBH_BdG*Vec(:,n); 
                end
            end
        end  
        DAmem = DAq(it);                 % memorise the D before the update
        DBmem = DBq(it);
        % Updating the gap
        DAq(it) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DAnn)));
        DBq(it) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DBnn)));
         
        discrA = abs((DAq(it)-DAmem))/abs((DAq(it)+DAmem)); % discrepancy indicator
        discrB = abs((DBq(it)-DBmem))/abs((DBq(it)+DBmem)); % discrepancy indicator
        discr = discrA+discrB; 
        disp(['it = ', num2str(iter), ' fin, 1000*DA = ', num2str(fix(DAq(it)*100000)/100),...
              ', 1000*DB = ', num2str(fix(DBq(it)*100000)/100),', discrepancy*1000 = ',num2str((discr*10^3))]);      
    end
    if abs(DAq(it)) < 10^(-5)
        DAq(it) = 0;
    end
    if abs(DBq(it)) < 10^(-5)
        DBq(it) = 0;
    end
    
end

plot(Tarr,DAq,Tarr,DBq)
xlabel('T')
legend('D_A','D_B')
% Comments:
% Convergence with the discr 10^(-2) and 10^(-3) was tested
% The difference is insignificant, but -2 is much faster
%