function [dat] = SBUTPAS(fl,fh, fs, ap, as)
%　斎藤正徳 1978 「漸化式ディジタル・フィルタの自動設計」
%  Butterworth型 Bandpass Filter 関数
%  入力値
%  FL:低周波側カットオフ周波数,FH:高周波側カットオフ周波数,
%  FS:ストップバンド周波数
%  AP:MAX Attendation in pass band , 
%  AS:min
%
%  出力値


% h = zeros(4, 1);
hp = pi/2.0;
wl = min(abs(fl), abs(fh))*pi;
wh = max(abs(fl), abs(fh))*pi;
ws = abs(fs)*pi;

if wl ~= 0 && wl ~= wh && wh < hp && ws ~= 0 ...
               && ws < hp && (ws-wl)*(ws-wh) > 0

    clh = 1./(cos(wl)*cos(wh));
    dp = sin(wh-wl)*clh;
    ww = tan(wl)*tan(wh);
    ts = tan(ws);
    ds = abs(ts-ww/ts);
    PA = min(abs(ap), abs(as));
    SA = max(abs(ap), abs(as));
    if PA == 0
        PA = 0.5;
    end
    if SA == 0
        SA = 5.0;
    end
    n  = max(2, round(abs(log(PA/SA)/log(dp/ds))));
    cc = exp(log(PA*SA)/double(n))/(dp*ds);
    C  = sqrt(cc);
    ww = ww*cc;
    dp = hp/double(n);
    K =  floor(n / 2);
    m = K*2;
    L  = 0;
    g  = 1.0;
    fj = 1.0;
    h = zeros(1, 4 * K);

    for j = 1:K
        dj = complex(cos(dp*fj),sin(dp*fj))*0.5;
        fj = fj+2;
        cq = sqrt(dj^2+ww);
        R(1) = dj +cq;
        R(2) = dj -cq;
        g = g*cc;
        for l = 1:2
            RE = real(R(l))^2;
            RI = imag(R(l));
            A = 1./((C+RI)^2+RE);
            g = g*A;
            h(L+1)= 0;
            h(L+2)= -1;
            h(L+3)= 2.*((RI-C)*(RI+C)+RE)*A;
            h(L+4)= ((RI-C)^2+RE)*A;
            L = L+4;
        end
    end
    gn = g;
    if n ~= m
        m = m + 1;
        WPC = cc*cos(wh-wl)*clh;
        WMC = -cc*cos(wh+wl)*clh;
        A = 1./(WPC + C);
        gn = g*C*A;
        h(L+1) = 0;
        h(L+2) = -1;
        h(L+3) = 2.*WMC*A;
        h(L+4) = (WPC-C)*A;
    end
else
    n = 0;
    m = 0;
    gn = 0.0;
end


dat.h = h;
dat.m = m;
dat.gn = gn;
dat.n = n;

dat.fl=num2str(fl);
dat.fh=num2str(fh);
dat.fs=num2str(fs);
dat.ap=num2str(ap);
dat.as=num2str(as);
dat.tp = 1;








