

%%
set(0,'DefaultFigureWindowStyle','docked')
clearvars

%% simulate alpha(lambda) = Eobs(lambda)/Eth with dl = 120 m



V0=300e-3;
l=linspace(20,4000,4000);
dl2=120;dl30=30;
i=1;
for lambda=l
    try
        kk=2*pi./lambda;
        V=@(x) V0*cos(kk.*x);
        Eth=@(x) V0*kk.*sin(kk.*x);
        Eobs=@(x,dl1) -(V(x+dl1)-V(x))/dl1;
        x=linspace(0,10*lambda,1000);
        
        alpha120(i)=max(Eobs(x,dl2))./max(Eth(x));
        alpha30(i)=max(Eobs(x,dl30))./max(Eth(x));
        
        [~,x1] = findpeaks(Eobs(x,dl2));
        [~,x2] = findpeaks(Eth(x));
        x1 = x1(2);
        x2 = x2(abs(x(x2)-x(x1)) == min(abs(x(x2)-x(x1))));
        dx = x(x2)-x(x1);
        phase120(i) = kk*dx;
        
        [~,x12] = findpeaks(Eobs(x,dl30));
        [~,x22] = findpeaks(Eth(x));
        x12 = x12(2);
        x22 = x22(abs(x(x22)-x(x12)) == min(abs(x(x22)-x(x12))));
        dx2 = x(x22)-x(x12);
        phase30(i) = kk*dx2;
    catch
        phase120(i) = nan;
        phase30(i) = nan;
    end
    i=i+1;
end
Ix = (find(diff(phase120)>pi/2));
j=1;
L=length(Ix);
if mod(L,2)~=0;Ix = [1 Ix];end
for i = L:-1:1
    if j == 1
        phase120(Ix(end-1):Ix(end)) = phase120(Ix(end-1):Ix(end)) + pi;
        j=0;
        Ix(end) = [];
    else
        j=1;
        Ix(end) = [];
    end
end
clear Ix j L
Ix = (find(diff(phase30)>pi/2));
j=1;
L=length(Ix);
if mod(L,2)~=0;Ix = [1 Ix];end
for i = L:-1:1
    if j == 1
        phase30(Ix(end-1):Ix(end)) = phase30(Ix(end-1):Ix(end)) + pi;
        j=0;
        Ix(end) = [];
    else
        j=1;
        Ix(end) = [];
    end
end

phase120(diff(phase120)>pi/2) = nan;
phase30(diff(phase30)>pi/2) = nan;
%% Generate electric and magnetic field data to check if the code works


B0=1;numf = 100; numk = 100;
freq1=logspace(2,3.55,numf);
wn=2*pi*freq1;
vph=90e3;Lmda=vph./freq1';

Thb=0:180;
Phid=0:360;
k_real = ones(length(Thb),length(Phid),numf,3)*nan;
k_measured = ones(length(Thb),length(Phid),numf,3)*nan;
k_maxvar_measured = ones(length(Thb),length(Phid),numf,3)*nan;
k_maxvar_120 = ones(length(Thb),length(Phid),numf,3)*nan;
eig_maxvar = ones(length(Thb),length(Phid),numf,3)*nan;
freq_wavelet = ones(length(Thb),length(Phid),numf)*nan;
ii = 1;mm=1;
for thb=Thb
    jj = 1;
    
    for phid=Phid
        %% Create synthetic data
        s=1e4;fs=8192;
        t=linspace(0,s/fs,s);
        B=B0*[sind(thb)*cosd(phid) sind(thb)*sind(phid) cosd(thb)];
        B=repmat(B,[s 1]);
        Bprobes=irf.ts_vec_xyz(t,B);
        b=mean(Bprobes.data)/norm(mean(Bprobes.data));
        
        
        
        k=2*pi*freq1'./vph.*[sind(thb)*cosd(phid) sind(thb)*sind(phid) cosd(thb)];k=k';
        k_real(ii,jj,1:length(freq1),1:3) = k';
        k_80th = rotz(45)*k;
        V0=300e-3;dl13=sqrt(2)*60;dl56=28.15;dl12=120;
        V=@(x,y,z,kx,ky,kz,w,t) V0*cos(kx.*x+ky.*y+kz.*z-w.*t);
        
        
        EThx=@(x,y,z,kx,ky,kz,w,t) V0.*kx.*sin(kx.*x+ky.*y+kz.*z-w.*t);
        EThy=@(x,y,z,kx,ky,kz,w,t) V0.*ky.*sin(kx.*x+ky.*y+kz.*z-w.*t);
        EThz=@(x,y,z,kx,ky,kz,w,t) V0.*kz.*sin(kx.*x+ky.*y+kz.*z-w.*t);
        
        EObx=@(dl,kx,ky,kz,w,t) - (V(dl/2,0,0,kx,ky,kz,w,t)-V(-dl/2,0,0,kx,ky,kz,w,t))/dl;
        EOby=@(dl,kx,ky,kz,w,t) - (V(0,dl/2,0,kx,ky,kz,w,t)-V(0,-dl/2,0,kx,ky,kz,w,t))/dl;
        EObz=@(dl,kx,ky,kz,w,t) - (V(0,0,dl/2,kx,ky,kz,w,t)-V(0,0,-dl/2,kx,ky,kz,w,t))/dl;
        
        for i=1:length(t)
            
            V1(i)=sum(V(dl12/2,0,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            V2(i)=sum(V(-dl12/2,0,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            V3(i)=sum(V(0,dl12/2,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            V4(i)=sum(V(0,-dl12/2,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            V5(i)=sum(V(0,0,dl56/2,k(1,:),k(2,:),k(3,:),wn,t(i)));
            V6(i)=sum(V(0,0,-dl56/2,k(1,:),k(2,:),k(3,:),wn,t(i)));
            
            
            ETh(i,1)=sum(EThx(0,0,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            ETh(i,2)=sum(EThy(0,0,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            ETh(i,3)=sum(EThz(0,0,0,k(1,:),k(2,:),k(3,:),wn,t(i)));
            
            
            EObs(i,1)=sum(EObx(dl12,k(1,:),k(2,:),k(3,:),wn,t(i)));
            EObs(i,2)=sum(EOby(dl12,k(1,:),k(2,:),k(3,:),wn,t(i)));
            EObs(i,3)=sum(EObz(dl56,k(1,:),k(2,:),k(3,:),wn,t(i)));
        end
        
        
        ETh(1:ceil(s*0.1),:)=[];
        EObs(1:ceil(s*0.1),:)=[];
        c_eval('V?(1:ceil(s*0.1))=[];',1:6);
        t(1:ceil(s*0.1))=[];
        V1=irf.ts_scalar(t,V1);
        V2=irf.ts_scalar(t,V2);
        V3=irf.ts_scalar(t,V3);
        V4=irf.ts_scalar(t,V4);
        V5=irf.ts_scalar(t,V5);
        V6=irf.ts_scalar(t,V6);
        ETTh = irf.ts_vec_xyz(t,ETh);
        
        SCpot = irf.ts_scalar(t,[V1.data V2.data V3.data V4.data V5.data V6.data]);
        T1 = [V1.time(1) V1.time(end)];
        
        %% Calculate different E fields
        
        out = mms.fk_powerspec_E80(SCpot,T1,'w0',4*5.36,'numf',numf,'numk',numk,'f',[100 4000],'return_fields',1,'correct_timeshifts',0);
        
        f = out{1}.f;
        freq_wavelet(ii,jj,1:numf) = f;
        k_80_x = out{1}.k_x;
        pow_80_x = out{1}.pow_x;
        k_80_x = [k_80_x+6*min(k_80_x) k_80_x+4*min(k_80_x) k_80_x+2*min(k_80_x) k_80_x k_80_x + 2*max(k_80_x) k_80_x + 4*max(k_80_x) k_80_x + 6*max(k_80_x)];
        pow_80_x = [pow_80_x;pow_80_x;pow_80_x;pow_80_x;pow_80_x;pow_80_x;pow_80_x];
        
        k_80_y = out{2}.k_y;
        pow_80_y = out{2}.pow_y;
        k_80_y = [k_80_y+6*min(k_80_y) k_80_y+4*min(k_80_y) k_80_y+2*min(k_80_y) k_80_y k_80_y + 2*max(k_80_y) k_80_y + 4*max(k_80_y) k_80_y + 6*max(k_80_y)];
        pow_80_y = [pow_80_y;pow_80_y;pow_80_y;pow_80_y;pow_80_y;pow_80_y;pow_80_y];
        
        E120 = out{3}.E120;
        E120 = irf_filt(E120,100,0,[],3);
        
        
        
        
        %% get frequency range
        %%%%%include a forloop over frequencies with F as loop variable
        Mxp_x=max(max(pow_80_x));
        Mxp_y=max(max(pow_80_y));
        
        
        thresh=1e-5;%% if signal has power less than this threshhold (whichis equal to digitization noise) then don't calculate
        if Mxp_x<=thresh || Mxp_y<=thresh
            a120x=ones(1,numf);a120y=ones(1,numf);k_sp=ones(numf,2)*nan;k2=ones(numf,3)*nan;k22=ones(numf,3)*nan;kmaxvar=ones(numf,3)*nan;kmaxvar2=ones(numf,3)*nan;L1=ones(numf,3)*nan;
            Lam=ones(1,numf)*nan;Lam2=ones(1,numf)*nan;Lam_sp=ones(1,numf)*nan;thkb=ones(1,numf)*nan;thkb2=ones(1,numf)*nan;
            p120x=ones(1,numf)*nan;p120y=ones(1,numf)*nan;
            ix=ones(1,numf)*nan;iy=ones(1,numf)*nan;
            k_sp80=ones(numf,2)*nan;Lamx_80r=ones(1,numf)*nan;Lamy_80r=ones(1,numf)*nan;
            
        else
            
            
            
            noise_lev=10^-5;
            
            flim(1,1:2) = [200 3300];
            klimx(1,1:2) =  [min(k_80_x) max(k_80_x)];
            klimy(1,1:2) =  [min(k_80_y) max(k_80_y)];
            %%
            a120x=ones(1,numf)*nan;a120y=ones(1,numf)*nan;
            p120x=ones(1,numf)*nan;p120y=ones(1,numf)*nan;
            ix=ones(1,numf)*nan;iy=ones(1,numf)*nan;
            k_sp=ones(numf,2)*nan;k2=ones(numf,3)*nan;k22=ones(numf,3)*nan;kmaxvar=ones(numf,3)*nan;kmaxvar2=ones(numf,3)*nan;L1=ones(numf,3)*nan;
            Lam=ones(1,numf)*nan;Lam2=ones(1,numf)*nan;Lam_sp=ones(1,numf)*nan;thkb=ones(1,numf)*nan;thkb2=ones(1,numf)*nan;
            k_sp80=ones(numf,2)*nan;Lamx_80r=ones(1,numf)*nan;Lamy_80r=ones(1,numf)*nan;
            
            
            for Fidx=1:numf
                
                %%calculate lambda for each case
                
                
                con = find(f(Fidx)>=flim(:,1) & f(Fidx) <=flim(:,2),1);
                
                if isempty(con)
                else
                    
                    con1x = (k_80_x>= min(klimx(con,:)) & k_80_x<= max(klimx(con,:)));
                    con1y = (k_80_y>= min(klimy(con,:)) & k_80_y<= max(klimy(con,:)));
                    pow_80_xtemp = pow_80_x(:,Fidx)';
                    pow_80_ytemp = pow_80_y(:,Fidx)';
                    
                    itempx=find((pow_80_xtemp ==max(pow_80_xtemp) & (con1x)));
                    k80xtemp = k_80_x(itempx);
                    ftemp = interp1(k_80th(1,:),freq1,k80xtemp);
                    ftemp = abs(ftemp-f(Fidx));
                    itempx = itempx(ftemp==min(ftemp));
                    clear ftemp
                    if isempty(itempx);continue;end
                    ix(Fidx)=itempx;
                    Lamx_80=2*pi/k_80_x(itempx);
                    Fm_80_x=f(Fidx);
                    Vp_80_x=Fm_80_x*Lamx_80;
                    
                    itempy=find((pow_80_ytemp ==max(pow_80_ytemp) & (con1y)));
                    k80ytemp = k_80_y(itempy);
                    ftemp = interp1(k_80th(2,:),freq1,k80ytemp);
                    ftemp = abs(ftemp-f(Fidx));
                    itempy = itempy(ftemp==min(ftemp));
                    clear ftemp
                    if isempty(itempy);continue;end
                    iy(Fidx)=itempy;
                    Lamy_80=2*pi/k_80_y(itempy);
                    Fm_80_y=f(Fidx);
                    Vp_80_y=Fm_80_y*Lamy_80;
                    
                    
                    
                    
                    %%make k vectors
                    
                    k_80=2*pi./[Lamx_80 Lamy_80];L80=2*pi/norm(k_80);
                    kn_80=k_80/norm(k_80);
                    
                    Rz=rotz(-45);Rz=Rz(1:2,1:2);%% put k in E120 coordinate instead of E80
                    k_80r=Rz*k_80';kn_80r=k_80r/norm(k_80r);
                    Lamx_80r(Fidx)=2*pi/k_80r(1);Lamy_80r(Fidx)=2*pi/k_80r(2);
                    
                    k_sp80(Fidx,1:2)=k_80;
                    k_sp(Fidx,1:2)=k_80r;
                    a120x(Fidx)=interp1(l,alpha120,abs(Lamx_80r(Fidx)));
                    p120x(Fidx)=interp1(l,phase120,abs(Lamx_80r(Fidx)));
                    
                    a120y(Fidx)=interp1(l,alpha120,abs(Lamy_80r(Fidx)));
                    p120y(Fidx)=interp1(l,phase120,abs(Lamy_80r(Fidx)));
                    
                    if isnan(a120x(Fidx)); a120x(Fidx)=1;p120x(Fidx) = 0;end
                    
                    if pow_80_x(itempx,Fidx)<noise_lev ; a120x(Fidx)=1;p120x(Fidx) = 0;end
                    
                    if isnan(a120y(Fidx)); a120y(Fidx)=1;p120y(Fidx) = 0;end
                    if pow_80_y(itempy,Fidx)<noise_lev ; a120y(Fidx)=1;p120y(Fidx) = 0;end
                    
                    
                    %%if power less than threshold remove measurement
                    if pow_80_x(itempx,Fidx)<noise_lev; k_sp(Fidx,1) = nan;k_sp80(Fidx,1) = nan;ix(Fidx)=nan;end
                    if pow_80_y(itempy,Fidx)<noise_lev; k_sp(Fidx,2) = nan;k_sp80(Fidx,2) = nan;iy(Fidx)=nan;end
                    
                    %%if wavelength along both directions is infinite then remove measurement
                    if abs(k_sp(Fidx,1)) <=0.0018 && abs(k_sp(Fidx,2)) <= 0.0018; k_sp(Fidx,1:2) = nan;k_sp80(Fidx,1:2) = nan;end
                    
                    clear itempx jtempx itempy jtempy
                end
            end
        end
        p120x(isnan(a120x)) = 0;
        p120y(isnan(a120y)) = 0;
        a120x(isnan(a120x)) = 1;
        a120y(isnan(a120y)) = 1;
        %%to account for when alpha~0, the 0.01 comes from trial and error
        a120x(a120x<0.01) = 0.01;
        a120y(a120y<0.01) = 0.01;
        %% perfomr fourier transform on E120 -correct amplitudes then apply ifft
        Fs=1/(E120.time(2)-E120.time(1));
        Len=length(E120.time);n=2^nextpow2(Len);
        F=0:(Fs/n):(Fs/2-Fs/n);
        X=fft(E120.x.data,n);
        Y=fft(E120.y.data,n);
        ff=flip(f);
        a120xf=flip(a120x);a120yf=flip(a120y);
        p120xf=flip(p120x);p120yf=flip(p120y);
        df=diff(ff);
        
        for i=1:length(ff)-1
            
            
            
            itemp1=find(F>=ff(i) & F<ff(i)+df(i));
            itemp2=n-itemp1+2;
            
            X(itemp1)=X(itemp1)/a120xf(i);
            X(itemp2)=X(itemp2)/a120xf(i);
            Y(itemp1)=Y(itemp1)/a120yf(i);
            Y(itemp2)=Y(itemp2)/a120yf(i);
            clear itemp1 itemp2
            
        end
        E120xc=ifft(X);E120yc=ifft(Y);
        E120c=E120;
        E120c.data(:,1)=(E120xc(1:Len));
        E120c.data(:,2)=(E120yc(1:Len));
        Ec=E120c;
        
        %% get thkb and lambda from Ec
        
        fs = Ec.time(2)-Ec.time(1);fs=1/fs;
        
        I = find(~isnan(k_sp(:,1)) & ~isnan(k_sp(:,2)));
        for i =  I'
            
            fmt=f(i)*1.01;
            if f(i)*1.01>=fs/2; fmt=0.999*fs/2;end
            p=0.99;
            Etemp = irf_filt(Ec,f(i)*p,fmt,[],3);%% Corrected E field in frequency interval of interest
            Etemp2 = irf_filt(E120,f(i)*p,fmt,[],3);%% Uncorrected E field in frequency interval of interest
            [~,l1,v] = irf_minvar(Etemp);%% apply minvar on E
            [~,l2,v2] = irf_minvar(Etemp2);%% apply minvar on E
            
            
            
            %%%getting thkb
            kmaxvar(i,:) = v(1,:);
            kmaxvar2(i,:) = v2(1,:);
            
            dp = interp1(l,phase120,abs(2*pi./k_sp(i,:)));
            kmaxvar(i,dp>pi) = -kmaxvar(i,dp>pi);
            
            L1(i,:) = l1;
            K1 = kmaxvar(i,1:2)/norm(kmaxvar(i,1:2));%%spin plane component of k_maxvar
            K2 = k_sp(i,:)/norm(k_sp(i,:));
            
            
            if abs(acosd(dot(K1,K2))) > 90 %% if K1 is antiparallel to K2 flip kmaxvar
                kmaxvar = -kmaxvar;
                kmaxvar2 = -kmaxvar2;
            end
            
            
            th = acosd(kmaxvar(i,3)/norm(kmaxvar(i,:)));%%angle between k_maxvar and z
            k2(i,:) = [k_sp(i,:) norm(k_sp(i,:))/sind(th)*cosd(th)];%%get full k from k_sp assuming no attenuation in z direction
            
            
        end
        k_measured(ii,jj,I,1:3) = k2(I,:);
        k_maxvar_measured(ii,jj,I,1:3) = kmaxvar(I,:);
        k_maxvar_120(ii,jj,I,1:3) = kmaxvar2(I,:);
        eig_maxvar(ii,jj,I,1:3) = L1(I,:);
        clearvars -except  l alpha120 alpha30 phase120 phase30 B0 numf numk freq1 wn vph Lmda Thb Phid k_real k_measured k_maxvar_measured k_maxvar_120 eig_maxvar freq_wavelet ii jj thb phid mm
        jj = jj + 1;
        mm = mm + 1;
        disp(mm/(length(Thb)*length(Phid))*100);
        
    end
    ii = ii + 1;
    
end

save('simulation_data.mat','k_real','k_measured','k_maxvar_measured','k_maxvar_120','eig_maxvar','freq_wavelet','freq1','Thb','Phid','Lmda','vph','l','alpha30','alpha120','phase120','phase30')