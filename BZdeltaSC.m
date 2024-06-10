 % This is the main file, it runs all the code
clc       % clean the command window
clear all % delete all 
%% Band structure computation
SimpleBZ

nqx = 50;
qx = linspace(-0.4,0.4,nqx);

DAq = qx*0+0.012; % seed value
DBq = qx*0+0.008; % seed value
Jxq = qx*0;
Jyq = qx*0;
% Ok value 1.1869e-05
dkxdky = (1/nx) * (1/ny) * abs(a1_star_strained(1)  * a2_star_strained(2)...
                              -a1_star_strained(2)  * a2_star_strained(1))/(2*pi)^2;
if dkxdky == 0
    cprintf('err','The BZ is crazy, fix the strain...\n');
    return
end

filename = 'progress.txt';
fid = fopen( filename, 'wt' );
fprintf( fid,'%d\n',0);
fclose(fid);

filenameR = 'results.txt';
fidR = fopen( filenameR, 'wt' );
fprintf(fidR,'q DA DB Jx Jy \n\n');
fclose(fidR);


cprintf('hyper','Starting the canlculation... Hold on... It"s not so simple...\n')% XX component (just a nice dislay of the statement)
l = 1;
for iqx = 1:nqx
    iter = 0;             % number of iterations from this q completed
    discrA = 1;           % just to inititalize while
    discrB = 1;
    energy = zeros (4, nx+1, ny+1);
    DAnn   = zeros (4, nx+1, ny+1);
    DBnn   = zeros (4, nx+1, ny+1);
    while (discrA > 10^(-3) && abs(DAq(iqx)) > 10^(-5))||(discrB > 10^(-3) && abs(DBq(iqx)) > 10^(-5))
        iter = iter + 1;  % update iteration index
        % Do not count the boundary twise! start from 0 and end at nx-1,
        % or from 1 and end at nx
        for i = 1 : nx
            for j = 1 : ny
                k = (i/nx) * a1_star_strained + (j/ny) * a2_star_strained;
                [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),qx(iqx),DAq(iqx),DBq(iqx)));
                energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
                for n = 1 : nbands
                    DAnn(n,i+1,j+1)  =  Vec(:,n)'*dDAH_BdG*Vec(:,n);
                    DBnn(n,i+1,j+1)  =  Vec(:,n)'*dDBH_BdG*Vec(:,n); 
                end
            end
        end

        DAmem = DAq(iqx);                 % memorise the D before the update
        DBmem = DBq(iqx);
        % Updating the gap
        DAq(iqx) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DAnn)));
        DBq(iqx) = - 0.5 * U * dkxdky * sum(sum(sum(fFD(energy).*DBnn)));
         
        discrA = abs((DAq(iqx)-DAmem))/abs((DAq(iqx)+DAmem)); % discrepancy indicator
        discrB = abs((DBq(iqx)-DBmem))/abs((DBq(iqx)+DBmem)); % discrepancy indicator   ', discrepancy = ',num2str(fix(discr*10^6)/10^6)]);      
    end
    if DAq(iqx) < 10^(-5)
       DAq(iqx) = 0;
    end
    if DBq(iqx) < 10^(-5)
       DBq(iqx) = 0;
    end
    for i = 0 : l*nx
        for j = 0 : l*ny
            k = (i/(l*nx)) * a1_star_strained + (j/(l*ny)) * a2_star_strained;
            [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),qx(iqx),DAq(iqx),DBq(iqx)));
            for n = 1 : nbands
                Jxq(iqx) = Jxq(iqx) + fFD(Val(n)).*Vec(:,n)'*dxH_BdG(k(1),k(2),qx(iqx))*Vec(:,n);
                %Jyq(iqx) = Jyq(iqx) + fFD(Val(n)).*Vec(:,n)'*dyH_BdG(k(1),k(2),qx(iqx))*Vec(:,n);
            end
        end
    end
    p = update_progress(filename);
    disp([num2str(fix(100*p/nqx)),'% finished'])

    fileID = fopen(filenameR,'a+');
    A = [qx(iqx); DAq(iqx); DBq(iqx); Jxq(iqx)* dkxdky/l; Jyq(iqx)* dkxdky/l];
    fprintf(fileID,'%12.8f %12.8f %12.8 %12.8 %12.8f\n\n',A);
    fclose(fileID);
end
Jxq = Jxq * dkxdky/l;

save('SCDEres0025.mat','qx','Jxq',"DAq","DBq","epsilon","th","t0","mu","U","Delta","DeltaE","Temp")

figure
plot(qx,DAq,qx,DBq)
xlabel('q_x')
legend('\Delta_A','\Delta_B')

figure
plot(qx,real(Jxq))
legend('J_x')
xlabel('q_x')
ylabel('J_q')



