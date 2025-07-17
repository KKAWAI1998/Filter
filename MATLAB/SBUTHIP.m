function [Fil] = SBUTHIP(fp, fs, ap, as)
%　斎藤正徳 1978 「漸化式ディジタル・フィルタの自動設計」
%  Butterworth型 Highpass Filter 関数
%  
%  

Fil.tp = 3;
Fil.fp=fp;
Fil.fs=fs;
Fil.ap=ap;
Fil.as=as;

h = zeros(4, 1);
hp = pi/2.0;

wp = max(abs(fp), abs(fs))*pi;
ws = min(abs(fp), abs(fs))*pi;

if ws ~= 0 && ws ~= wp && wp < hp
    tp = tan(wp);
    ts = tan(ws);
    pa = min(abs(ap), abs(as));
    sa = max(abs(ap), abs(as));
    if pa == 0
        pa = 0.5;
    end
    if sa == 0
        sa = 5.0;
    end
    n = max(2, fix(abs(log(sa/pa)/log(tp/ts))+0.5));
    cc = exp(log(pa*sa)/double(n))*(tp*ts);
    c = sqrt(cc);

    dp = hp/double(n);
    m = fix(n/2);
    k = m*4;
    g = 1.0;
    fj = 1.0;
    c2 = -2.0*(1.0 - c)*(1.0 + c);

    

    for j = 1 : 4 : k
        sj = cos(dp*fj)^2;
        tj = sin(dp*fj);
        fj = fj + 2;
        a = 1.0/((c + tj)^2 + sj);
        g = g*a;
        h(j) = -2.0;
        h(j+1) = 1.0;
        h(j+2) = c2*a;
        h(j+3) = ((c-tj)^2 + sj)*a;
    end
    gn = g;
    if mod(n, 2) ~= 0
        m = m + 1;
        gn = g/(c+1.0);
        h(k+1) = -1.0;
        h(k+2) = 0.0;
        h(k+3) = (c - 1.0)/(c+1.0);
        h(k+4) = 0.0;
    end
else
    n = 0;
    m = 0;
    gn = 0.0;
end


Fil.h = h;
Fil.m = m;
Fil.gn = gn;
Fil.n = n;
