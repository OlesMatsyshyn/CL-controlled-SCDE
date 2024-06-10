SimpleBZ_con % initialization

%% Band structure computation
% Compute the band structure and collect the results
% Structure of array below as follows

% Electronic Bands        
band   = zeros (2, 2, nx+1, ny+1);
energy = zeros (2, nx+1, ny+1);

kx = linspace(0,1,nx+1);
ky = linspace(0,1,ny+1);
[KX,KY] = meshgrid(ky,kx);

for i = 0 : nx
    for j = 0 : ny
        k = Kp + (i/nx) * a1_star_strained + (j/ny) * a2_star_strained;
        [Vec,Val] = eigenshuffle(H_eff(k(1),k(2)));
        band(:,:,i+1,j+1) = Vec;       % column eigenvectors at which kx, at which ky, ;
        energy(:,i+1,j+1) = Val(:);    % band, at which kx, at which ky;
    end
end

nbands = 2;
unnormalizationfactorBand = zeros(1,nbands);
unorthogonalityfactorBand = zeros(1,nbands);
uneigenalityfctorBand = zeros(1,nbands);

for n = 1 : nbands
    for i = 0 : nx
        for j = 0 : ny
            k = Kp + (i/nx) * a1_star_strained + (j/ny) * a2_star_strained;
            unnormalizationfactorBand(n) = unnormalizationfactorBand(n) + abs(1 - band(:,n,i+1,j+1)'*band(:,n,i+1,j+1));
            for m = 1 : nbands
                if m == n
                    continue;
                end
                unorthogonalityfactorBand(n) = unorthogonalityfactorBand(n) + abs(band(:,n,i+1,j+1)'*band(:,m,i+1,j+1));
            end    
            uneigenalityfctorBand(n) = uneigenalityfctorBand(n) + sum(abs(H_eff(k(1),k(2))*band(:,n,i+1,j+1) - energy(n,i+1,j+1)*band(:,n,i+1,j+1)));
        end
    end
end
% Full Basis
%band(:,1,54,78)*band(:,1,54,78)'+band(:,2,54,78)*band(:,2,54,78)'+band(:,3,54,78)*band(:,3,54,78)'+band(:,4,54,78)*band(:,4,54,78)';
if (sum(unnormalizationfactorBand) + sum(unorthogonalityfactorBand) + sum(uneigenalityfctorBand)> 10^(-5))
    cprintf('err','Too large numerical error...\n');
    return
end
%% Display f
fmapr = zeros(nx+1,ny+1);
fmapi = zeros(nx+1,ny+1);
for i = 0 : nx
    for j = 0 : ny
        k = Kp + (i/nx) * a1_star_strained + (j/ny) * a2_star_strained;
        fmapr(i+1,j+1) = abs(f_strained(k(1),k(2)));
        fmapi(i+1,j+1) = imag(f_strained(k(1),k(2)));
    end
end
%figure 
%imagesc(fmapr)
%% Spectrum 3D and the FS
if seeBD == 1
    figure
    %set(gcf,'Position',[810 310 800 800])
    % Plot the Band structure at mu = 0
    % figure
    surf(KX,KY,reshape(energy(1,:,:),nx+1,ny+1),'EdgeColor','none')
    hold on
    surf(KX,KY,reshape(energy(2,:,:),nx+1,ny+1),'EdgeColor','none')
    surf(KX,KY,mu+0*reshape(energy(1,:,:),nx+1,ny+1),'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[0, 0, 0])
    % axis equal
    axis off
    zlim([min(min(energy(2,:,:))) max(max(energy(1,:,:)))])
    D_eff = min(min(energy(1,:,:)))-max(max(energy(2,:,:)));
    xlabel('$\vec{a}^*_{1,strained}$','Interpreter','latex')
    ylabel('$\vec{a}^*_{2,strained}$','Interpreter','latex')
    figure
    set(gcf,'Position',[10 310 800 800])
    energyslice = zeros(nx+1,ny+1,2);
    for i = 1 : 2
        enl = energy(i,:, :);
        enl(enl<mu-0.001) = NaN;
        enl(enl>mu+0.001) = NaN;
        energyslice(:,:,i) = double(~isnan(enl));
    end
    imagesc(kx,ky,sum(energyslice,3))
    xlabel('$\vec{a}^*_{1,strained}$','Interpreter','latex')
    ylabel('$\vec{a}^*_{2,strained}$','Interpreter','latex')
    axis equal
end
% Periodicity check
% Band
p_check = norm(reshape(energy(1,1,:)-energy(1,end,:),ny+1,1)); % Looks good!
% Hamiltonian
pHx = 0;
pHy = 0;
% over x
for j = 0 : ny
        k = Kp + (0/nx) * a1_star_strained + (j/ny) * a2_star_strained;
        pHx = pHx + sum(sum(abs(H_eff(k(1),k(2)))));
        k = Kp + (nx/nx) * a1_star_strained + (j/ny) * a2_star_strained;
        pHx = pHx - sum(sum(abs(H_eff(k(1),k(2)))));
end
% over y
for i = 0 : nx
        k = Kp + (i/nx) * a1_star_strained + (0/ny) * a2_star_strained;
        pHy = pHy + sum(sum(abs(H_eff(k(1),k(2)))));
        k = Kp + (i/nx) * a1_star_strained + (ny/ny) * a2_star_strained;
        pHy = pHy - sum(sum(abs(H_eff(k(1),k(2)))));
end
if (pHx + pHy > 10^(-5))
    cprintf('err','Periodicity is violated...\n');
    return
end