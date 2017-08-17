function f,B,x
 
  return,B[0]*x + B[1]
end

function quad_func,const,x1,x2,m1,m2
  return, (m1*x1)+(m2*x2)+const
end

function try, x
  return, x
end

;function MY_OODR, x, y, sx=sx,sy=sy,beta0=[1.,2.]
;  python.import,"scipy.odr"
;  x=[1,2,3,4]
;  y=x
;  sx=stddev(x)
;  sy=stddev(y)
;  beta0=[1.,2.]
;
;  mydata=python.Realdata(python.wrap(list(X1,Y,sx,sy)))
;  mybeta=python.wrap(List(beta0))
;  mymodel=python.Model('quad_func')
;
;  py_odr=python.import('scipy.odr')
;  myodr=py_odr.odr
;  myresult=myodr(mymodel,beta=[1.,2.],mydata)
;  return, myresult
;end

pro try_ODR

; TSI_data= read_latis_data('nrl2_observational_composite_tsi_v02r01', '2003-06-15', '2003-07-30', query='replace_missing(-9999999)')
;  SSD_data= read_latis_data('nrl2_sunspot_darkening_v02r01', '2003-06-15','2003-07-30')
;  MGII_data= read_latis_data('nrl2_facular_brightening_v02r01', '2003-06-15','2003-07-30')
;  TIME_data=TSI_data.Time
;  UNCERTAINTY_data=TSI_data.TSI_UNCERTAINTY
;
;  TSI_filter=TSI_data.TSI
;  SSD_Filter=SSD_data.SSD
;  MGII_filter=MGII_data.MGII
;
;  convert_Time=make_array(size=size(TIME_data))
;  n=size(TIME_data,/N_Elements)-1
;  for j=0,n do convert_Time[j]=iso_date2mjdn(TIME_data[j])
;  for j=0,n do convert_Time[j]=convert_Time[j]
;  
;  i=where(TSI_filter NE -9999999)
;  TSI= TSI_filter[i]
;  UNCERTAINTY=UNCERTAINTY_data[i]
;  SSD=SSD_filter[i]
;  MGII=MGII_filter[i]
;  final_Time=float(convert_Time[i])
;
;  TSI_quiet=1360.45
;  SSD_quiet=0
;  MGII_quiet=0.1502
;  
;  Delta_TSI=TSI - TSI_quiet
;  Delta_SSD=SSD - SSD_quiet
;  Delta_MGII=MGII- MGII_quiet
;  
;  X1=[Delta_SSD]
;  X2=[Delta_MGII]
;  Y=[Delta_TSI]
;  X=[TRANSPOSE(X1), TRANSPOSE(X2)]
;  x_err=20% of the relative validation
;  y_err=0.
;  
;  
;  result = REGRESS(X, Y, SIGMA=sigma, CONST=const,MEASURE_ERRORS=UNCERTAINTY, yfit=yfit)
  stop
  python.import,"scipy.odr"
  x=[1,2,3,4]
  y=x
  sx=stddev(x)
  sy=stddev(y)
  beta0=[1.,2.]

  mydata=python.Realdata(python.wrap(list(x,y,sx,sy)))
  mybeta=python.wrap(List(beta0))
  ananomousfunc=lambda(quad_func:quad_func)
  mymodel=python.Model(f([1,1],1))
  myodr=python.ODR(python.Model('f'),beta0,x,y)

  
  py_odr=python.import('scipy.odr')
  myodr=py_odr.ODR
  myresult=myodr([mymodel,beta0,mydata])
  
  
end





