i.loxodrom.dist <-
function(x1, y1, x2, y2, epsilon=0.0001){
dis<-numeric(length(x1))
rerde<-6368
deltax<-abs(x2*pi/180-x1*pi/180)
deltay<-abs(y2*pi/180-y1*pi/180)
tga<-deltax/(log(tan(pi/4+y2*pi/360))-log(tan(pi/4+y1*pi/360))) 

dis[abs(x1-x2)<epsilon&abs(y1-y2)<epsilon]<-0
dis[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]<-abs(cos(y1[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]*pi/180)*deltax[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)])
dis[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]<-abs(deltay[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]/cos((pi-atan(tga[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]))))
dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]<--deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]))
dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]<-abs(deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)])))
dis[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)])))
dis[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)])))
dis*rerde
}

