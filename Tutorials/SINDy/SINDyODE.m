function dy = SINDyODE(t,y,param)

yPool = poolData(y',length(y),param.polyorder);
dy = (yPool*param.Xi)';