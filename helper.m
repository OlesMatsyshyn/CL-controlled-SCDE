for iqx = 1:nqx
    for i = 0 : l*nx
        for j = 0 : l*ny
            k = (i/(l*nx)) * a1_star_strained + (j/(l*ny)) * a2_star_strained;
            [Vec,Val] = eigenshuffle(H_BdG(k(1),k(2),qx(iqx),DAq(iqx),DBq(iqx)));
            for n = 1 : nbands
                Jxq(iqx) = Jxq(iqx) + fFD(Val(n)).*Vec(:,n)'*dxH_BdG(k(1),k(2),qx(iqx))*Vec(:,n);
                Jyq(iqx) = Jyq(iqx) + fFD(Val(n)).*Vec(:,n)'*dyH_BdG(k(1),k(2),qx(iqx))*Vec(:,n);
            end
        end
    end
end