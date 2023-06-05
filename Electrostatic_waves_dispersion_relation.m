
set(0,'DefaultFigureWindowStyle','docked')
clearvars

ic=1;%%spacecraft number

%% set times and load data
%%time of the waveburst of interest
T1=irf.tint('2021-01-10T17:27:14.48/2021-01-10T17:27:14.63');

%%extended time around the waveburst
Tint  = T1 +[-1 1];



c_eval('Egse=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
c_eval('Bgse=mms.get_data(''B_gse_brst_l2'',Tint,?);',ic);
c_eval('defatt?=mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',ic);
c_eval('defatt?.zdec=mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',ic);
c_eval('Vi = mms.get_data(''Vi_gse_fpi_brst_l2'',Tint+[-60 60],?);',ic);
c_eval('Ne = mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ic);
c_eval('Te = mms.get_data(''Te_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('Ti = mms.get_data(''Ti_gse_fpi_brst_l2'',Tint,?);',ic);

dir=mms.db_list_files('mms1_edp_brst_l2_scpot',T1);dir_name=dir.name;
dir=[dir.path '/' dir_name];
d=dataobj(dir);
d=d.data;
V=d.mms1_edp_dcv_brst_l2;V=V.data;V_2=V;
Time=d.mms1_edp_epoch_brst_l2;Time=Time.data;
SCpot=TSeries(EpochUnix(Time),V);
zphase = mms.db_get_variable('mms1_ancillary_defatt','zphase',Tint);


%% Getting the Zphase in the probes coordinate system (PCS)

zphase=irf.ts_scalar(zphase.time,zphase.zphase);
zphase=zphase.tlim(T1+[-1 1]);


zphase=zphase.resample(Egse)+30;%%resample and make the angle so the inverse rotation align the x direction with probe 1
zphase=zphase.tlim(T1);



%% puting the fields in DSL coordinate system
c_eval('Bdsl = mms_dsl2gse(Bgse,defatt?,-1);', ic);
c_eval('Edsl = mms_dsl2gse(Egse,defatt?,-1);', ic);
c_eval('Vidsl = mms_dsl2gse(Vi,defatt?,-1);', ic);
Bdsl=Bdsl.resample(Edsl);Vidsl=Vidsl.resample(Edsl);


%% Spin plane interferometry
numf = 100; numk = 100;
out = mms.fk_powerspec_E80(SCpot,T1,'w0',4*5.36,'numf',numf,'numk',numk,'f',[100 4000],'return_fields',1);

f = out{1}.f;

k_80_x = out{1}.k_x;
pow_80_x = out{1}.pow_x;
%%extend the domain in wavenumber to mitigate aliasing
k_80_x = [k_80_x+4*min(k_80_x) k_80_x+2*min(k_80_x) k_80_x k_80_x + 2*max(k_80_x) k_80_x + 4*max(k_80_x)];
pow_80_x = [pow_80_x;pow_80_x;pow_80_x;pow_80_x;pow_80_x];

k_80_y = out{2}.k_y;
pow_80_y = out{2}.pow_y;
%%extend the domain in wavenumber to mitigate aliasing
k_80_y = [k_80_y+4*min(k_80_y) k_80_y+2*min(k_80_y) k_80_y k_80_y + 2*max(k_80_y) k_80_y + 4*max(k_80_y)];
pow_80_y = [pow_80_y;pow_80_y;pow_80_y;pow_80_y;pow_80_y];

E120 = out{3}.E120;
E120 = irf_filt(E120,100,0,[],3);
%% Rotate B and V to PCB

Bprobes=Bdsl.tlim(T1);
Viprobes = Vidsl.tlim(T1);
tic
Zphse=mean(zphase.data);
for i=1:length(Bprobes.data)
    Bdsltemp=Bprobes.data(i,:);
    Vidsltemp=Viprobes.data(i,:);
    Bdslt(1:3,i)=rotz(-Zphse)*Bdsltemp';
    Vidslt(1:3,i)=rotz(-Zphse)*Vidsltemp';
    clear Bdsltemp
end

toc
Bprobes.data=Bdslt';
Viprobes.data=Vidslt';
b=mean(Bprobes.data)/norm(mean(Bprobes.data));%% background B in PCB
clear Bdslt Vidslt





%% simulate alpha(lambda) = Eobs(lambda)/Eth and \Delta \phi with dl = 120 m and 30 m



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
%%account for discontinuitites in phase difference and wrap it so it would
%%be from 0 to 2pi
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

%% get frequency and wavenumber ranges of dispersion relations in thespin plane
%%%%%include a forloop over frequencies with F as loop variable
Mxp_x=max(max(pow_80_x));
Mxp_y=max(max(pow_80_y));


thresh=1e-5;%% if signal has power less than this threshhold (which 
%is equal to digitization noise) then don't calculate
if Mxp_x<=thresh || Mxp_y<=thresh
    a120x=ones(1,numf);a120y=ones(1,numf);k_sp=ones(numf,2)*nan;k2=ones(numf,3)*nan;k22=ones(numf,3)*nan;
    Lam=ones(1,numf)*nan;Lam2=ones(1,numf)*nan;Lam_sp=ones(1,numf)*nan;thkb=ones(1,numf)*nan;thkb2=ones(1,numf)*nan;
else
    
    Vidshp=mean(Viprobes.data);
    Vidsh0 = dot(mean(Viprobes.data),b)*b; %%assuming k||b
    VE80 = rotz(45)*Vidsh0';
    
    noise_lev=10^-5.5;
    
    fh=figure;
    h80x=subplot('Position',[0.07 0.25 0.4 0.4]);
    h80y=subplot('Position',[0.55 0.25 0.4 0.4]);
    surf(h80x,k_80_x,f,log10(pow_80_x)')
    shading(h80x,'interp')
    view(h80x,0,90)
    axis(h80x,'tight')
    xlabel(h80x,'k 1/m')
    ylabel(h80x,'f Hz')
    title(h80x,['E80x' newline 'Vdsh_{80x} = ' num2str(round(VE80(1))) ' km/s'])
    colorbar(h80x)
    c=caxis(h80x);
    caxis(h80x,[log10(noise_lev) c(2)])
    hold(h80x,'on')
    
    surf(h80y,k_80_y,f,log10(pow_80_y)')
    shading(h80y,'interp')
    view(h80y,0,90)
    axis(h80y,'tight')
    xlabel(h80y,'k 1/m')
    ylabel(h80y,'')
    yticklabels(h80y,'')
    title(h80y,['E80y' newline 'Vdsh_{80y} = ' num2str(round(VE80(2))) ' km/s'])
    colorbar(h80y)
    caxis(h80y,c)
    caxis(h80y,[log10(noise_lev) c(2)])
    hold(h80y,'on')
    colormap jet
    drawnow
    
    
    %%
    nn = input('How many dispersion relations do you see?');
    fidx = [];kidx_x = [];kidx_y = [];
    a = annotation(fh,'textbox',[0.3 0.5 0.3 0.3],'String',['Choose f limits for DR # ' num2str(1) ],'FitBoxToText','on','FontSize',30);
    for ii = 1:nn
        if ii ~= 1
            a.String = ['Choose f limits for DR # ' num2str(ii)];
        end
        gtemp = ginput(2);
        flim(ii,1:2) = gtemp(:,2);
        clear gtemp
        a.String = ['Choose kx limits for DR # ' num2str(ii)];
        gtemp = ginput(2);
        klimx(ii,1:2) = gtemp(:,1);
        a.String = ['Choose ky limits for DR # ' num2str(ii)];
        gtemp = ginput(2);
        klimy(ii,1:2) = gtemp(:,1);
        
    end
%% Calculate the spin plane dispersion relation
    a120x=ones(1,numf)*nan;a120y=ones(1,numf)*nan;
    p120x=ones(1,numf)*nan;p120y=ones(1,numf)*nan;
    ix=ones(1,numf)*nan;iy=ones(1,numf)*nan;
    k_sp=ones(numf,2)*nan;k2=ones(numf,3)*nan;k22=ones(numf,3)*nan;
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
            
            itempx=find((pow_80_xtemp ==max(pow_80_xtemp) & (con1x)),1);
            if isempty(itempx);continue;end
            ix(Fidx)=itempx;
            Lamx_80=2*pi/k_80_x(itempx);
            Fm_80_x=f(Fidx);
            Vp_80_x=Fm_80_x*Lamx_80;
            
            itempy=find((pow_80_ytemp ==max(pow_80_ytemp) & (con1y)),1);
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

%% plot detected maxima to fk plots
IX = ix(~isnan(ix));IY = iy(~isnan(iy));
fx = f(~isnan(ix));fy = f(~isnan(iy));
Px = ones(size(ix))*nan;Py = ones(size(iy))*nan;
for i = 1:length(IX)
    
    plot3(h80x,k_80_x(IX(i)),fx(i),log10(pow_80_x(IX(i),(f==fx(i))))','*k')
    Px((ix == IX(i))) = (pow_80_x(IX(i),(f==fx(i))));
    
end
for i = 1:length(IY)
    plot3(h80y,k_80_y(IY(i)),fy(i),log10(pow_80_y(IY(i),(f==fy(i))))','*k')
    Py((iy == IY(i))) = (pow_80_y(IY(i),(f==fy(i))));
end

Pt = log10(Px+Py);
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
for i = I'
    
    fmt=f(i)*1.01;
    if f(i)*1.01>=fs/2; fmt=0.999*fs/2;end
    Etemp = irf_filt(Ec,f(i)*0.99,fmt,[],3);%% Corrected E field in frequency interval of interest
    Etemp2 = irf_filt(E120,f(i)*0.99,fmt,[],3);%% Uncorrected E field in frequency interval of interest
    
    [~,l1,v] = irf_minvar(Etemp);%% apply minvar on corrected E
    [~,l2,v2] = irf_minvar(Etemp2);%% apply minvar on uncorrected E
    
    
    %%%getting thkb
    kmaxvar = v(1,:);
    kmaxvar2(i,:) = v2(1,:);
    
    %%deal with phase shift effect on maximum variance direction. If phase
    %%shift due to short wavelength effects is larger than pi the electric
    %%field flips direction, so we flip it again here to account for it.
    dp = interp1(l,phase120,abs(2*pi./k_sp(i,:)));
    kmaxvar(dp>pi) = -kmaxvar(dp>pi);
    
    K1 = kmaxvar(1:2)/norm(kmaxvar(1:2));%%spin plane component of k_maxvar
    K2 = k_sp(i,:)/norm(k_sp(i,:));%%spin plane component from interferometry
    
    if abs(acosd(dot(K1,K2))) > 90 %% if K1 is antiparallel to K2 flip kmaxvar i.e remove ambiguity in minimum variance anallysis
        kmaxvar = -kmaxvar;
        kmaxvar2(i,:) = -kmaxvar2(i,:);
    end
    
    
    th = acosd(kmaxvar(3)/norm(kmaxvar));%%angle between k_maxvar and z
    k2(i,:) = [k_sp(i,:) norm(k_sp(i,:))/sind(th)*cosd(th)];%%get full k from k_sp assuming no attenuation in z direction
    thkb(i) = acosd(dot(b,k2(i,:)/norm(k2(i,:))));
     
    
    %%%getting lambda
    
    Lam(i) = 2*pi/norm(k2(i,:));
    
    
    
end


k2=k2(I,:);
Lam = Lam(I);



%% Compare with simulation
load('simulation_data.mat')
%
knorm = k2./vecnorm(k2')';
theta_kz = acosd(knorm(:,3));

phi_kz = atan2d(knorm(:,2)./vecnorm(knorm(:,1:2)')',knorm(:,1)./vecnorm(knorm(:,1:2)')');
phi_kz(phi_kz<0) = 360 + phi_kz(phi_kz<0) ;

TTkr = acosd(dot(k_real,k_measured,4)./(vecnorm(k_real,2,4).*vecnorm(k_measured,2,4)));%%angle between real and measured direction of propagation 
TTkr(TTkr>90) = 180 - TTkr(TTkr>90);
clear Bb
Bb(1,1,1,1:3) = b;%%background B field 
Bb = (repmat(Bb,size(k_real,1),size(k_real,2),size(k_real,3)));
TTkb = acosd(dot(k_real,Bb,4)./(vecnorm(k_real,2,4).*vecnorm(Bb,2,4)));%% angle between background B field and real k


TTkz120 = acosd(k_maxvar_120(:,:,:,3)./(vecnorm(k_maxvar_120,2,4)));%% polar angle of the maximum variance direction of the uncorrecte E field 
TTkz120(TTkz120>90) = 180 - TTkz120(TTkz120>90);


Phkx120 = atan2d(k_maxvar_120(:,:,:,2)./(vecnorm(k_maxvar_120(:,:,:,1:2),2,4)),k_maxvar_120(:,:,:,1)./(vecnorm(k_maxvar_120(:,:,:,1:2),2,4)));%%azimuthal angle of the maximum variance direction of the uncorrected E field
Phkx120 = abs(Phkx120);%%for E120 maximum variance I put the azimuthal angle in the first quadrant because of the pi ambiguity
Phkx120(Phkx120>90) = 180 - Phkx120(Phkx120>90);

TTkzm = acosd(k_measured(:,:,:,3)./(vecnorm(k_measured,2,4))); %% polar angle of the measured wave vector


Phkxm= atan2d(k_measured(:,:,:,2)./(vecnorm(k_measured(:,:,:,1:2),2,4)),k_measured(:,:,:,1)./(vecnorm(k_measured(:,:,:,1:2),2,4)));
Phkxm(Phkxm<0) = 360 + Phkxm(Phkxm<0);


L_real = vph./freq_wavelet;%%get wavelength at frequency of the wavelet transform
L_measured = 2*pi./vecnorm(k_measured,2,4);
L_ratio = L_real./L_measured;

ii=zeros(size(I))*nan;jj=ii;kk=ii;
ii2=zeros(size(I))*nan;jj2=ii2;kk2=ii2;
m=1;
for i = 1:length(I)
    try
        Ltemp = Lam(i);
        thtemp = theta_kz(i);
        phitemp = phi_kz(i);
        
        th120temp = acosd(kmaxvar2(I(i),3));th120temp(th120temp>90) =180 - th120temp(th120temp>90);
        phi120temp = atan2d(kmaxvar2(I(i),2)./norm(kmaxvar2(I(i),1:2)),kmaxvar2(I(i),1)./norm(kmaxvar2(I(i),1:2)));
        
        phi120temp = abs(phi120temp);%%for E120 maximum variance I put the azimuthal angle in the first quadrant because of the pi ambiguity
        phi120temp(phi120temp>90) = 180 - phi120temp(phi120temp>90);
        
        
        
        dphitemp = abs(Phkxm - phitemp);    dthtemp = abs(TTkzm-thtemp);dltemp = abs(L_measured - Ltemp);
        dth120temp = abs(TTkz120 - th120temp); dphi120temp = abs(Phkx120 - phi120temp);
        dphi120temp(dphi120temp>90) = 180-dphi120temp(dphi120temp>90);
        
        mdphi = 5; mdth = 5; mdl = 5; mdth120 = 10; mdphi120 = 10;
        Ii = find(dphitemp < mdphi & dthtemp < mdth & dltemp < mdl & dth120temp < mdth120 & dphi120temp < mdphi120);
        
        
        
        dphitemp2 = dphitemp(Ii)/mdphi;
        dthtemp2 = dthtemp(Ii)/mdth;
        dltemp2 = dltemp(Ii)/mdl;
        dth120temp2 = dth120temp(Ii)/mdth120;
        dphi120temp2 = dphi120temp(Ii)/mdphi120;
        c=1;lim = 1e-2;Ii2=[];
        while c
            Ii2 = find(dphitemp2<=lim & dthtemp2<=lim & dltemp2<=lim & dth120temp2<=lim & dphi120temp2<=lim,1);
            lim = lim +1e-2;
            if ~isempty(Ii2) || lim>1
                c = 0;
            end
        end
        
        
        [iit2,jjt2,kkt2] = ind2sub(size(TTkzm),Ii((Ii2)));
        [iit,jjt,kkt] = ind2sub(size(TTkzm),Ii);
        ii(m,1:length(iit)) = iit;
        jj(m,1:length(iit)) = jjt;
        kk(m,1:length(iit)) = kkt;
        ii2(m,1:length(iit2)) = iit2;
        jj2(m,1:length(iit2)) = jjt2;
        kk2(m,1:length(iit2)) = kkt2;
        m=m+1;
        
        
        
        clear Ltemp thtemp phitemp th120temp phi120temp dphitemp dthtemp dltemp dth120temp dphi120temp Ii dphitemp2...
            dthtemp2 dltemp2 dth120temp2 dphi120temp2 Ii2
    catch e
        ii(m,:) = nan;
        jj(m,:) = nan;
        kk(m,:) = nan;
        m=m+1;
        
    end
end

%%
%%%variables for the solution that minimizes all quantities in equation 6 in the manuscript simultaneously
tkr2 = ones(size(ii,1),1)*nan;
tkb2 = ones(size(ii,1),1)*nan;
Lrat2 = ones(size(ii,1),1)*nan;
Lreal2 = ones(size(ii,1),1)*nan;
FF = ones(size(ii))*nan;
%%%variables for the all solutions that satisfy the condition in equation 6
%%%in the manuscript
tkr = ones(size(ii))*nan;
tkb = ones(size(ii))*nan;
Lm = ones(size(ii))*nan;
Lrat = ones(size(ii))*nan;
Lreal = ones(size(ii))*nan;
kr = ones([size(ii) 3])*nan;

for i = 1:size(ii,1)
    try
        
        tkr2(i) = TTkr(ii2(i),jj2(i),kk2(i));
        tkb2(i) = TTkb(ii2(i),jj2(i),kk2(i));
        Lrat2(i) = L_ratio(ii2(i),jj2(i),kk2(i));
        Lreal2(i) = L_real(ii2(i),jj2(i),kk2(i));
    catch
        tkr2(i) = nan;
        tkb2(i) = nan;
        Lrat2(i) = nan;
        Lreal2(i) = nan;
    end
    
    for j = 1:size(ii,2)
        try
            FF(i,j) = freq_wavelet(ii(i,j),jj(i,j),kk(i,j));
            tkr(i,j) = TTkr(ii(i,j),jj(i,j),kk(i,j));
            tkb(i,j) = TTkb(ii(i,j),jj(i,j),kk(i,j));
            Lm(i,j) = L_measured(ii(i,j),jj(i,j),kk(i,j));
            Lrat(i,j) = L_ratio(ii(i,j),jj(i,j),kk(i,j));
            Lreal(i,j) = L_real(ii(i,j),jj(i,j),kk(i,j));
            kr(i,j,:) = k_real(ii(i,j),jj(i,j),kk(i,j),:)./norm(squeeze(k_real(ii(i,j),jj(i,j),kk(i,j),:)));
        catch
            FF(i,j) = nan;
            tkr(i,j) = nan;
            tkb(i,j) = nan;
            Lm(i,j) = nan;
            Lrat(i,j) = nan;
            Lreal(i,j) = nan;
            kr(i,j,:) = nan;
        end
    end
end


%% plot
w120=irf_wavelet(E120,'nf',100,'returnpower',1,'f',[100 4000],'wavelet_width',4*5.36);
p = w120.p; p = p{1}+p{2}+p{3};
w120.p = p;
noise_lev=10^-5;fs=20;


h= figure;
set(h,'Units','normalized','Position',[0.1 0 0.6 1])

h1 = subplot('Position',[0.1 0.86 0.72 0.11]);
h5 = subplot('Position',[0.1 0.75 0.72 0.11]);
h2 = subplot('Position',[0.1 0.64 0.72 0.11]);

h3 = subplot('Position',[0.18 0.35 0.23 0.23]);
h4 = subplot('Position',[0.51 0.35 0.23 0.23]);

h6 = subplot('Position',[0.18 0.05 0.23 0.23]);
h7 = subplot('Position',[0.51 0.05 0.23 0.23]);

irf_plot(h1,E120.z,'r')
hold(h1,'on')
irf_plot(h1,E120.y,'b')
irf_plot(h1,E120.x,'k')
h1.XTickLabel=[];
text(h1,0.7,0.84,'E_x','FontSize',fs,'Units','normalized','BackgroundColor','w','color','k');
text(h1,0.75,0.84,'E_y','FontSize',fs,'Units','normalized','BackgroundColor','w','color','b');
text(h1,0.8,0.84,'E_z','FontSize',fs,'Units','normalized','BackgroundColor','w','color','r');
text(h1,0.7,0.2,'Uncorrected','FontSize',fs,'Units','normalized','BackgroundColor','w','color','k');



irf_plot(h5,E120c.x,'k')
hold(h5,'on')
irf_plot(h5,E120c.y,'b')
irf_plot(h5,E120c.z,'r')
h5.XTickLabel=[];
text(h5,0.7,0.84,'E_x','FontSize',fs,'Units','normalized','BackgroundColor','w','color','k');
text(h5,0.75,0.84,'E_y','FontSize',fs,'Units','normalized','BackgroundColor','w','color','b');
text(h5,0.8,0.84,'E_z','FontSize',fs,'Units','normalized','BackgroundColor','w','color','r');
text(h5,0.7,0.2,'Corrected','FontSize',fs,'Units','normalized','BackgroundColor','w','color','k');

irf_spectrogram(h2,w120)
h2c = colorbar(h2);
h2c.Label.String = ['log10 E^2' newline '(mV/m)^2 Hz^{-1}'];
h2c.Label.FontSize = fs;
h2c.Position = [0.85 0.64 0.01 0.11];
Tint = [E120.time(1) E120.time(end)];
irf_zoom([h1 h2 h5],'x',Tint);
irf_plot_axis_align([h1 h2 h5]);
ylabel(h1,'E (mV/m)','FontSize',fs)
ylabel(h5,'E (mV/m)','FontSize',fs)
ylabel(h2,'f (kHz)','FontSize',fs)
h2.XLabel.FontSize= fs;
set(h2,'YScale','log')


surf(h3,k_80_x,f,log10(pow_80_x)')
shading(h3,'interp')
view(h3,0,90)
hold(h3,'on')

for i = 1:length(IX)
    
    plot3(h3,k_80_x(IX(i)),fx(i),log10(pow_80_x(IX(i),(f==fx(i))))','*k')
    
end
axis(h3,'tight')
xlabel(h3,'k_{x80} (m^{-1})','FontSize',fs)
ylabel(h3,'f (Hz)','FontSize',fs)



surf(h4,k_80_y,f,log10(pow_80_y)')
shading(h4,'interp')
hold(h4,'on')
for i = 1:length(IY)
    plot3(h4,k_80_y(IY(i)),fy(i),log10(pow_80_y(IY(i),(f==fy(i))))','*k')
end
view(h4,0,90)
axis(h4,'tight')
xlabel(h4,'k_{y80} (m^{-1})','FontSize',fs)
ylabel(h4,'')
yticklabels(h4,'')
h4c = colorbar(h4);
c=caxis(h4);
caxis(h4,[log10(noise_lev) c(2)])
caxis(h2,[log10(noise_lev) c(2)])
caxis(h3,[log10(noise_lev) c(2)])

h4c.Position = [0.75 0.35 0.01 0.23];
h4c.Label.String = ['log10 E^2' newline '(mV/m)^2 Hz^{-1}'];
h4c.Label.FontSize = fs;
colormap(gcf,'jet')
box([h1 h2 h3 h4 h5],'on')




text(h1,0.9,0.85,'(a)','FontSize',fs,'Units','normalized','BackgroundColor','w');
text(h5,0.9,0.85,'(b)','FontSize',fs,'Units','normalized','BackgroundColor','w');
text(h2,0.9,0.85,'(c)','FontSize',fs,'Units','normalized','BackgroundColor','w');
text(h3,0.9,0.9,'(d)','FontSize',fs,'Units','normalized','BackgroundColor','w');
text(h4,0.9,0.9,'(e)','FontSize',fs,'Units','normalized','BackgroundColor','w');

[iix,jjx] = find(Lreal<45);
Lreal(iix,jjx) = nan;tkb(iix,jjx) = nan;
kr(iix,jjx,:) = nan;
mL = nanmedian(Lreal');
sL = nanstd(Lreal');
mtkb = nanmedian(tkb');
stkb = nanstd(tkb');
mkr = squeeze(nanmedian(kr,2));
mkr = mkr./vecnorm(mkr')';
%


for i = 1:nn
    Itemp = find(f>=min(flim(i,:)) & f<= max(flim(i,:)) & ~isnan(k_sp(:,1)) & ~isnan(k_sp(:,2)));
    
end





slope = [];ft = fittype('a/x');

for i = 1:nn
    I2 = find(f(I)>=min(flim(i,:)) & f(I)<= max(flim(i,:)) & ~isnan(k_sp(I,1)) & ~isnan(k_sp(I,2)));
    plot(h6,tkb((I2),:),f(I(I2)),'*b')
    hold(h6,'on')
    errorbar(h6,mtkb(I2),f(I(I2)),zeros(size(I2)),zeros(size(I2)),stkb(I2),stkb(I2),'*-k','linewidth',2)
    
    
    
    p0 = plot(h7,Lreal((I2),:),f(I(I2)),'*b');
    hold(h7,'on')
    ixt = find(~isnan(mL(I2)));
    p1 = errorbar(h7,mL(I2),f(I(I2)),zeros(size(I2)),zeros(size(I2)),sL(I2),sL(I2),'*-k','linewidth',3);
    fittemp=fit(mL(I2(ixt))',f(I(I2(ixt))),ft);
    slope = fittemp.a;
    clear fittemp
    if i ==1
        p2 = plot(h7,slope./f(I(I2)),f(I(I2)),'r','linewidth',2);
        
        text(h7,0.5,mean(f(I(I2)))/max(f(I)),['f = (' num2str(round(slope/1000)) ' km/s)/\lambda '],'Units','normalized','FontSize',fs,'color','r')
    else
        p3 = plot(h7,slope./f(I(I2)),f(I(I2)),'g','linewidth',2);
        text(h7,0.5,mean(f(I(I2)))/(1.1*max(f(I))),['f = (' num2str(round(slope/1000)) ' km/s)/\lambda '],'Units','normalized','FontSize',fs,'color','g')
    end
end
grid([h6 h7],'on')
xlabel(h6,'\theta_{kB}^\circ','FontSize',fs)
xlabel(h7,'\lambda (m)','FontSize',fs)
h7.YTickLabel=[];
ylabel(h6,'f (Hz)','FontSize',fs)
text(h6,0.85,0.85,'(f)','FontSize',fs,'Units','normalized')
text(h7,0.85,0.85,'(g)','FontSize',fs,'Units','normalized')
if nn>1
    leg = legend([p0(1) p1 p2 p3],'All solutions','median \pm \sigma of all solutions','f_1 = V_{ph1}/\lambda_1','f_2 = V_{ph2}/\lambda_2','interpreter','tex');
else
    leg = legend([p0(1) p1 p2],'All solutions','median \pm \sigma of all solutions','f = V_{ph}/\lambda','interpreter','tex');
end
xlim(h6,[0 180])
set(leg,'Position',[0.75 0.1 0.17 0.1])