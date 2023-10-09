function dy = SINDyODEsine(t,y,param)

yPool = poolDataSine(y',length(y),param.polyorder,param.usesine,param.sineorder);
dy = (yPool*param.Xi)';