function [dat] = SRECRES(h, m, gn,fm,df, n)

if m > 0 && n > 0
    K = 4*m;
    C = gn ^2;
    CS = cos(fm*pi*2);
    SN = sin(fm*pi*2);
    DC = -2.*sin(df * pi)^2;
    DS = sin(df*pi*2);
    G = zeros(1, n);
    P = zeros(1, n);

    for i = 1:n
        GG = C;
        PP = 0;
        for j = 1:4:K
            A = (h(j+1)+1.0)*CS+h(j);
            B = (h(j+1)-1.0)*SN;
            GG = GG * (A^2 + B^2);
            if GG > 0
                PP = PP + atan2(B,A);
                A = (h(j+3)+1.0)*CS+h(j+2);
                B = (h(j+3)-1.0)*SN;
                GG = GG /(A^2 + B^2);
                PP = PP - atan2(B,A);
                if abs(PP) > pi
                    PP = PP - r_sign(2*pi,PP);
                end
            end
        end
        G(i) = GG;
        P(i) = PP;
        W = CS;
        CS = CS*DC-SN*DS+CS;
        SN = W*DS+SN*DC+SN;
        W = (1.0-CS^2-SN^2)*0.5;
        CS = CS*W+CS;
        SN = SN*W+SN;
    end
end

dat.G = G;
dat.P = P;
end

function result = r_sign(a, b)
% Example r_sign function definition
% This function returns the value of a with the sign of b
if b >= 0
    result = abs(a);
else
    result = -abs(a);
end
end