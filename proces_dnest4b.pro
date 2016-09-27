pro proces_DNest4b

cd ,'/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy/'
dir = '/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy'

;
info_filename = dir + '/'   + 'posterior_sample.txt'
temp_sample = READ_ASCII(info_filename)
posterior_sample = temp_sample.field00001

n = 14    ; number of parameters 
total_sample =posterior_sample[0:14,*]  

titles = ['x','y','mag','re','nser', 'axrat','ang','box','mag-bar', 'rout','a','b','q-bar','ang-bar','box-bar']


nil = calc_stats(total_sample, titles)


index = [2,3,4,5,6]

;# obtain subset without images     [0xcen, 1ycen, 2mag, 3re, 4nser, 5axrat, 6ang, 7box, mag8, rout9,a10, b11, axrat12, ang13, box14] 
titels = ['x','y','mag','re','nser', 'axrat','ang','box','mag-bar', 'rout','a','b','q-bar','ang-bar','box-bar']
;range=[[0,100],[0,100],[0,40], [0,100],[0,10],[0,1],[0,180],[-1,1],[0,40],[0,1],[0,40],[0,40],[0,1],[0,180], [-1,1] ] 
range=[[0,100],[0,100],[15,20], [10,25],[4,8],[0,1],[0,180],[-1,1],[15,20],[0.4,1.2],[-2,9],[-2,3],[0,1],[0,180], [-1,1] ]
range=[[0,100],[0,100],[17.,19.], [15.,19],[5.,8],[0.1,0.5],[0,90],[-1,1],[17.,19.],[0.4,1.2],[-2,9],[-2,3],[0,1],[0,180], [-1,1] ]

sub_sample = posterior_sample[index,*]
sub_titles = titles[index]
subrange= range[*,index]
; Set up the display window.
;cgDisplay, 9*n, 6*n
;!P.Multi = [0, n+1, n+1]
;!X.OMargin = [0, 25]
;!Y.OMargin = [0, 8]



image_name = strcompress(dir + '/' + 'triangle_disc_IDL'  + '.eps' , /remove_all)
n = n_elements(index)
;pos = cgLayout([n, n], OXMargin=[n, n], OYMargin=[n, n], XGap=5, YGap=8)
;pos = cgLayout([n, n], OXMargin=[n, n], OYMargin=[n, n], XGap=4, YGap=7)
pos = cgLayout([n, n],  OXMargin=[2, 2], OYMargin=[2, 2], XGap=3, YGap=3)


cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
cgDisplay, 1200, 1200

;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
loadct, 33


z =n
;pos = cgLayout([z,z])
pos = cgLayout([z,z], OXMargin=[4, 3], OYMargin=[4, 3], XGap=4, YGap=4)

FOR j=0,z-1  DO BEGIN
  for i= 0L, j do begin
    if (i eq  j ) then begin
      cghistoplot, sub_sample[i,*], NoErase=j NE 0, Position=pos[*,i+z*j],xtitle=sub_titles[i],Title=None, charsize =1, xrange=[subrange[0,i],subrange[1,i]], xticks=3
    endif else begin
      cgplot, sub_sample[i,*], sub_sample[j,*], NoErase=j NE 0, Position=pos[*,i+z*j],xtitle=sub_titles[i],ytitle=sub_titles[j],psym=2, charsize =1,  xrange=[subrange[0,i],subrange[1,i]],yrange=[subrange[0,j],subrange[1,j]] , xticks=3
    endelse
  endfor
ENDFOR



; Clear the multiple plot variable.
!P.Multi = 0
  
   
;plot, output_data[*,0]
cgPS_Close


cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
image_name = strcompress(dir + '/' + 'triangle_disc_IDL'  + '.png' , /remove_all)
STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
spawn, strinG


angle = 90

n = n_elements(titels)
for i = 0, n-1 do begin
     for j = (i+1), n-1 do begin
            
          if (i ne j) then begin
                 
                 
                  image_name = strcompress(dir + '/' + 'scatterplot_'+   titels[i] + '_ '+ titels[j]  +'.eps' , /remove_all)       
                  cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
                 ; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
                  cgDisplay, 1200, 1200      
                  ;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
                  loadct, 33
                  cgplot,  posterior_sample[i,*], posterior_sample[j,*], xtitle=titles[i],ytitle=titles[j], psym=3 

                  cgPS_Close
                  
                  
                  cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
                  image_name = strcompress(dir + '/' + 'scatterplot_'+   titels[i] + '_ '+ titels[j]  +'.png' , /remove_all)                  
                  STRING =  'convert '+ image_name +' -rotate -90'  + ' ' +  image_name
                  spawn, strinG
                   


         endif
   
  
  endfor
  
  

endfor





index = [8,9,10,11, 12,13]
sub_sample = posterior_sample[index,*]
sub_titles = titles[index]
subrange= range[*,index]
; Set up the display window.
;cgDisplay, 9*n, 6*n
;!P.Multi = [0, n+1, n+1]
;!X.OMargin = [0, 25]
;!Y.OMargin = [0, 8]



image_name = strcompress(dir + '/' + 'triangle_bar_IDL'  + '.eps' , /remove_all)
n = n_elements(index)
;pos = cgLayout([n, n], OXMargin=[n, n], OYMargin=[n, n], XGap=5, YGap=8)
;pos = cgLayout([n, n], OXMargin=[n, n], OYMargin=[n, n], XGap=4, YGap=7)
pos = cgLayout([n, n],  OXMargin=[2, 2], OYMargin=[2, 2], XGap=3, YGap=3)


cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
cgDisplay, 1200, 1200

;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
loadct, 33


z =n
;pos = cgLayout([z,z])
pos = cgLayout([z,z], OXMargin=[4, 3], OYMargin=[4, 3], XGap=4, YGap=4)

FOR j=0,z-1  DO BEGIN
  for i= 0L, j do begin
    if (i eq  j ) then begin
      cghistoplot, sub_sample[i,*], NoErase=j NE 0, Position=pos[*,i+z*j],xtitle=sub_titles[i],Title=None, charsize =1, xrange=[subrange[0,i],subrange[1,i]],xticks=3
    endif else begin
      cgplot, sub_sample[i,*], sub_sample[j,*], NoErase=j NE 0, Position=pos[*,i+z*j],xtitle=sub_titles[i],ytitle=sub_titles[j],psym=3, charsize =1,  xrange=[subrange[0,i],subrange[1,i]],yrange=[subrange[0,j],subrange[1,j]],xticks=3
      ;cgplot, sub_sample[i,*], sub_sample[j,*], NoErase=j NE 0, Position=pos[*,i+z*j],xtitle=sub_titles[i],ytitle=sub_titles[j],psym=2, charsize =1,  xrange=[subrange[0,i],subrange[1,i]],yrange=[subrange[0,j],subrange[1,j]],XTickFormat='exponent', YTickFormat='exponent'
    endelse
  endfor
ENDFOR



; Clear the multiple plot variable.
!P.Multi = 0


;plot, output_data[*,0]
cgPS_Close


cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
image_name = strcompress(dir + '/' + 'triangle_bar_IDL'  + '.png' , /remove_all)
STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
spawn, strinG





stop 

end



function calc_stats, data, titels
                           
;confidence_rate 95
;
 

s = size(data)
ndim = s[1]
len = s[2]

for i = 0L, ndim-1 do begin
   h = data[i,*]           ; select one row
   m = median(h)           ; calculate median
   index = sort(h)         ; create index
   h = h[index]
   left_index =    (len/100)*2.5
   right_index =    (len/100)*(100-2.5)
   ;print,titles[i], 'median',m, 'credible interval [',h[left_index],',',h[right_index],']'  
   ; print,titles[i], 'median',m ,' + ',h[right_index]-m, ' - ',m-h[left_index]
   
   print,FORMAT='(A0, "$ ",F7.3, "^{+", F5.3,"}_{-",F5.3,"}$ &")',titels[i], m,h[right_index]-m,m-h[left_index]
   ;PRINT, FORMAT='("You are ", I0, " years old, ", A0)', age, name

endfor
return, 0

end


FUNCTION Exponent, axis, index, number

  ; A special case.
  IF number EQ 0 THEN RETURN, '0'
  
  ; Assuming multiples of 10 with format.
  ex = String(number, Format='(e8.0)')
  pt = StrPos(ex, '.')
  
  first = StrMid(ex, 0, pt)
  sign = StrMid(ex, pt+2, 1)
  thisExponent = StrMid(ex, pt+3)
  
  ; Shave off leading zero in exponent
  WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)
  
  ; Fix for sign and missing zero problem.
  IF (Long(thisExponent) EQ 0) THEN BEGIN
    sign = ''
    thisExponent = '0'
  ENDIF
  
  ; Make the exponent a superscript.
  IF sign EQ '-' THEN BEGIN
    RETURN, first + 'x10!U' + sign + thisExponent + '!N'
  ENDIF ELSE BEGIN
    RETURN, first + 'x10!U' + thisExponent + '!N'
  ENDELSE
  
END