function table_write_tf(sysd,wi,wf,step)

pw=sysd;

dado=[];

for ww=wi:step:wf
    w1=1i*ww;
    m1=abs(evalfr(pw,w1));
    m1db=20*log10(abs(evalfr(pw,w1)));
    f1=angle(evalfr(pw,w1))*180/pi;
    if f1>0,
        f1=f1-360;
    end
    dado=[dado; ww m1 m1db f1];
end


[l1,c1]=size(dado);

fprintf('       w        Amplitude       dB        Fase\n')
fprintf('                                        \n')
for xx=1:l1,
    fprintf('  %8.2f     %8.4f     %8.4f  %8.4f\n',dado(xx,:));
end


