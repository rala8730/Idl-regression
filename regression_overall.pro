function mul_regression_line,x1,x2,m1,m2,b
  y=(m1*x1)+(m2*x2)+b
  data_x_y=[[x1],[x2],[y]]
  return,data_x_y
end

Function GetDebSSD, data1, data2, data3
  i=0
  SSD=0
  time=''
  foreach element1, data2 do begin
    foreach element2, data1 do begin
      if element1 eq element2 then begin
        SSD=[SSD,data3[i]]
        time=[Time,element1]
      endif
    endforeach
    i=i+1
  endforeach
  SSD=SSD[1:*]
  Time=Time[1:*]
  data={Time:Time,SSD:SSD}
  return,data
end

function getmodel,TSI_filter, SSD_Filter, MGII_Filter, time_data, UNCERTAINTY_data

convert_Time=make_array(size=size(TIME_data))
n=size(TIME_data,/N_Elements)-1
for j=0,n do convert_Time[j]=iso_date2mjdn(TIME_data[j])
for j=0,n do convert_Time[j]=convert_Time[j]

i=where(TSI_filter NE -9999999)
TSI= TSI_filter[i]
UNCERTAINTY=UNCERTAINTY_data[i]
SSD=SSD_filter[i]
MGII=MGII_filter[i]
final_Time=float(convert_Time[i])

TSI_quiet=1360.45
SSD_quiet=0
MGII_quiet=0.1502

Delta_TSI=TSI - TSI_quiet
Delta_SSD=SSD - SSD_quiet
Delta_MGII=MGII- MGII_quiet

X1=[Delta_SSD]
X2=[Delta_MGII]
Y=[Delta_TSI]
X=[TRANSPOSE(X1), TRANSPOSE(X2)]

result = REGRESS(X, Y, SIGMA=sigma, CONST=const,MEASURE_ERRORS=UNCERTAINTY, yfit=yfit)
;PRINT, 'Constant: ', const
;PRINT, 'Coefficients: ', result[*]
;PRINT, 'Standard errors: ', sigma ;obs-real error for tsi

;p=plot(yfit)
;p1=plot(const+result[0,0]*X1+result[0,1]*X2,'r',/over)

regline_x1_y=mul_regression_line(x1,x2,result[0],result[1],const)
plt_TSI_point=plot(final_Time,Delta_TSI,color='red',XTITLE="TIME",YTITLE="$\Delta$ TSI",TITLE='$\Delta$ TSI vs Time',FONT_SIZE='7.3',MARGIN=0.099,Layout=[2,2,1])
t1=text(53000,-3,"Constant:"+string(const),/data)
t1=text(53000,-2.5,"Coefficients:"+string(result[0] + result[1]),/data)
t1=text(53000,-2,"Standard errors:"+string(sigma),/data)
plt_TSI_prime_point=plot(final_Time,regline_x1_y[*,2],color='blue',XTITLE="TIME",YTITLE="Model TSI",TITLE='TSI_prime vs TIME',FONT_SIZE='7.3',MARGIN=0.099,layout=[2,2,3],/current)
residual=plot(final_time,(Delta_TSI-regline_x1_y[*,2]),color='green',TITLE='Residual',XTITLE="TIME",YTITLE="MeasureTSI- ModelTSI",FONT_SIZE='7.3',MARGIN=0.099,layout=[2,2,2],/current)
standard_deviation=plot(final_time,(Delta_TSI-regline_x1_y[*,2])^2,color='purple',Title='Standard Deviation',XTITLE="TIME",YTITLE="(MeasureTSI- ModelTSI)^2",FONT_SIZE='7.3',MARGIN=0.099,layout=[2,2,4],/current)
Correlation=correlate([X1,X2],Y)
print, correlation,"corr"
stop
end 

function sorcetsi
  ;Getting data from lisird/latis and filtering with date and adding -9999999 for missing value to filter the index
  TSI_data=read_latis_data('sorce_tsi_24hr_l3','2003-02-25','2016-12-31',query='replace_missing(-9999999)')
  SSD_data= read_latis_data('nrl2_sunspot_darkening_v02r01', '2003-02-25','2016-12-31')
  MGII_data=read_latis_data('nrl2_facular_brightening_v02r01', '2003-02-25','2016-12-31')
  time_data=mjd2iso_date(jd2mjd(TSI_data.time))
  UNCERTAINTY_data=TSI_data.MEASUREMENT_UNCERTAINTY_1AU

  TSI_filter=TSI_data.TSI_1AU
  SSD_filter=SSD_data.SSD
  MGII_filter=MGII_data.MGII

  TSI_24hr=getmodel(TSI_filter,SSD_filter, MGII_Filter, time_data, UNCERTAINTY_data)
end

function debrecen
  TSI_data= read_latis_data('nrl2_observational_composite_tsi_v02r01', '2003-02-25','2016-12-31', query='replace_missing(-9999999)')
  MGII_data= read_latis_data('nrl2_facular_brightening_v02r01', '2003-06-15','2016-12-31')

  TIME_data=TSI_data.Time
  UNCERTAINTY_data=TSI_data.TSI_UNCERTAINTY

  TSI_filter=TSI_data.TSI
  MGII_filter=MGII_data.MGII

  file='~/Desktop/nrlssi/data/debrecen_sunspot_blocking_corrected_heliographic_1975-01-01_2016-12-31.txt'
  RESTORE, 'myPlotTemplate.sav',/verb
  Data=read_ascii(file,Template=plotTemplate)

  data1=TSI_data.Time
  data2=Data.Field1
  data3=Data.Field3

  Deb_SSD=GetDebSSD(data1,data2,data3)
  modeldeb=getmodel(TSI_filter,Deb_SSD.SSD, MGII_Filter,Deb_SSD.Time, UNCERTAINTY_data)
end

function observational
  TSI_data= read_latis_data('nrl2_observational_composite_tsi_v02r01', '2003-06-15', '2003-07-30', query='replace_missing(-9999999)')
  SSD_data= read_latis_data('nrl2_sunspot_darkening_v02r01', '2003-06-15','2003-07-30')
  MGII_data= read_latis_data('nrl2_facular_brightening_v02r01', '2003-06-15','2003-07-30')
  TIME_data=TSI_data.Time
  UNCERTAINTY_data=TSI_data.TSI_UNCERTAINTY

  TSI_filter=TSI_data.TSI
  SSD_Filter=SSD_data.SSD
  MGII_filter=MGII_data.MGII
  obs=getmodel(TSI_filter,SSD_Filter, MGII_Filter, time_data, UNCERTAINTY_data)
end

pro regression_overall
  
debrecen=debrecen()
sorcetsi=sorcetsi()
;observational=observational()
end