############## DUI
mul_iso<-c(971,427,136,95,348,1025,300,286)
mgsg<-c(90,51,13,6,18,94,16,13)
mg<-c(1,4,0,0,2,3,5,1)
sg<-c(9,2,2,5,5,11,5,4)

all<-rbind( mul_iso,mgsg,mg,sg)
colnames(all)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')
all
 
res_all<-chisq.test(all,correct=F)
res_all$expected


mgsg_1<-rbind( mul_iso,mgsg)
colnames(mgsg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')

res_mgsg<-chisq.test(mgsg_1,correct=F)
res_mgsg$expected

mg_1<-rbind( mul_iso,mg)
colnames(mg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')


res_mg<-chisq.test(mg_1,correct=F)
res_mg$expected


sg_1<-rbind( mul_iso,sg)
colnames(sg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')

res_sg<-chisq.test(sg_1,correct=F)
res_sg$expected


################# DEMG
mul_iso<-c(971,427,136,95,348,1025,300,286)
mgsg<-c(170,73,17,16,51,194,44,34)
mg<-c(7,1,0,0,2,10,1,1)
sg<-c(15,10,3,2,5,32,6,3)

all<-rbind( mul_iso,mgsg,mg,sg)
colnames(all)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')
all

res_all<-chisq.test(all,correct=F)
res_all$expected
res_all

mgsg_1<-rbind( mul_iso,mgsg)
colnames(mgsg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')

res_mgsg<-chisq.test(mgsg_1,correct=F)
res_mgsg$expected
res_mgsg

mg_1<-rbind( mul_iso,mg)
colnames(mg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')


res_mg<-chisq.test(mg_1,correct=F)
res_mg$expected
res_mg

sg_1<-rbind( mul_iso,sg)
colnames(sg_1)<-c('ES','MES','MEE','IR','ATTS','ATSS','A5','A3')

res_sg<-chisq.test(sg_1,correct=F)
res_sg$expected
res_sg
