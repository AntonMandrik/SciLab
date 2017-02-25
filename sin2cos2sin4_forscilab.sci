clear
nmtr=5;

l0=20;
A=2.2;
g=9.81;
s=poly(0,'s');
den=[1 0.1 g./l0];
nom=[1];
nomw=poly(nom,'s','c');
denw=poly(den,'s','c');
i=sqrt(-1);
pi=3.1415;
ci=1;

W=syslin('c',nomw,denw);
ci=1;
wst=0.05;
wmax=0.2;
G1p=cell(1,7);

//ñîçäàåì ìàòðèöó íà äèàãîíàëè êîòîðîé ñòîÿò ïåðåäàòî÷íûå ôóíêöèè îò nw
//îñòàëüíûå ýëåìåíòû ìàòðèöû - íîëþ.
clear den1;
for wn=1:1:nmtr+1
den1(wn,3)=den(1).*(wn-1).*(wn-1);
den1(wn,2)=-den(2).*(wn-1);
den1(wn,1)=den(3);
den1c(wn)=poly(den1(wn,:),'s','c');
end

clear den2;
for wn=1:1:nmtr+1
den2(wn,3)=den(1).*(wn-1).*(wn-1);
den2(wn,2)=den(2).*(wn-1);
den2(wn,1)=den(3);
den2c(wn)=poly(den2(wn,:),'s','c');
end

for Wwi=0:1:nmtr
Ww(nmtr-Wwi+1,nmtr-Wwi+1)=syslin('c',[1],den1c(Wwi+1));
Ww(nmtr+Wwi+1,nmtr+Wwi+1)=syslin('c',[1],den2c(Wwi+1));
end
//%ñîçäàåì ìàòðèöó ñ ýëåìåíòàìè -w äî ñðåäíåé ñòðîêè è w ïîñëå ñðåäíåé ñòðîêè
for wi=1:1:(nmtr*2+1)
Wneg(wi,wi)=wi-nmtr-1;
end
//ñîçäàåì ìàòðèöó ñ ýëåìåíòàìè -w äî ñðåäíåé ñòðîêè è w ïîñëå ñðåäíåé ñòðîêè


f11=0:0.1:6.3; 

for ppp=[2 4 6 8 10]
for wi=1:1:round((wmax/wst)) 
 
for fi=1:1:length(f11) 
    f=f11(fi);
   
    
Wsin2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wsin2wf(ni-2,ni)=(A./2).*exp(i.*0.5.*pi-i.*f);
Wsin2wf(ni,ni-2)=-(A./2).*exp(i.*0.5.*pi+i.*f);
end
Wsin2wf=-g.*(l0.^-2).*Wsin2wf;

f1=f+0.5.*pi./2;
Wcos2wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=3:1:(nmtr*2+1)
Wcos2wf(ni-2,ni)=(A./2).*exp(i.*0.5.*pi-i.*f1);
Wcos2wf(ni,ni-2)=-(A./2).*exp(i.*0.5.*pi+i.*f1);
end
Wcos2wf=4*(l0^-1)*(-1)*wi*wi*Wneg*Wcos2wf*i;

f2=2.*f;
Wsin4wf=zeros(nmtr*2+1,nmtr*2+1);
for ni=5:1:(nmtr*2+1)
Wsin4wf(ni-4,ni)=(A./2).*exp(i.*0.5.*pi-i.*f2);
Wsin4wf(ni,ni-4)=-(A./2).*exp(i.*0.5.*pi+i.*f2);
end
Wsin4wf=-2.*(l0.^-2)*(-1)*wi*wi*(A.^1)*Wneg*Wsin4wf*i;

Wa=Wsin2wf+Wcos2wf+Wsin4wf;

for gi=1:1:nmtr*2+1
    for gk=1:1:nmtr*2+1 
      Www(gi,gk)=repfreq(Ww(gi,gk),wi.*wst);
    end
end

G=-Www*Wa;

mag(wi,ppp,fi)=abs(((G^ppp)(nmtr+2,nmtr+2))^(1/ppp));


end
end

    for wi=1:1:round((wmax./wst))
    Mf(wi,ppp)=f11(1);
    for fi=2:1:length(f11)
    
        if (mag(wi,ppp,1)<mag(wi,ppp,fi))
        mag(wi,ppp,1)=mag(wi,ppp,fi);
        Mf(wi,ppp)=f11(fi);
        end
    
    end
    end

ppp

//wi=1:1:round((wmax./wst));
//plot(wi*wst,mag(:,ppp,1),'color',[250*0.8*(1.1+sin(ci))/1.7 250*0.8*(1.1+sin(ci+4*pi/3))/17 250*0.8*(1.1+sin(ci+2.*pi./3))/1.7]/255)
//ci=ci+(6.28./10);
end



