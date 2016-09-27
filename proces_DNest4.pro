pro proces_DNest4

cd ,'/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy/'
dir = '/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy'

;
info_filename = dir + '/'   + 'posterior_sample.txt'
temp_sample = READ_ASCII(info_filename)
posterior_sample = temp_sample.field00001

n = 14    ; number of parameters 


index = [2,3,4,5,6]

;# obtain subset without images     [0xcen, 1ycen, 2mag, 3re, 4nser, 5axrat, 6ang, 7box, mag8, rout9,a10, b11, axrat12, ang13, box14] 
titles = ['x','y','mag','re','nser', 'axrat','ang','box','mag-bar', 'rout','a','b','q-bar','ang-bar','box-bar'] 
sub_sample = posterior_sample[index,*]
sub_titles = titles[index]
; Set up the display window.
;cgDisplay, 9*n, 6*n
;!P.Multi = [0, n+1, n+1]
;!X.OMargin = [0, 25]
;!Y.OMargin = [0, 8]



image_name = strcompress(dir + '/' + 'triangle_'  + '.eps' , /remove_all)



n = n_elements(index)
pos = cgLayout([n, n], OXMargin=[n, n], OYMargin=[n, n], XGap=5, YGap=8)



cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
cgDisplay, 1200, 1200

;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
loadct, 33


z =n
pos = cgLayout([z,z])

FOR j=0,z-1  DO BEGIN
  for i= 0L, j do begin
    if (i eq  j ) then begin
      cghistoplot, sub_sample[i,*], NoErase=j NE 0, Position=pos[*,i+z*j], Title='Plot ' + StrTrim(j+1,2)
    endif else begin
      cgplot, sub_sample[i,*], sub_sample[j,*], NoErase=j NE 0, Position=pos[*,i+z*j], Title='Plot ' + StrTrim(j+1,2)
    endelse
  endfor
  
  
ENDFOR


    
    ;    cgText, 0.5, 0.925, /Normal, 'Example Plot Layout', Alignment=0.5, Charsize=cgDefCharsize()*1.25

;; !P.MULTI = [0, n, n]            ; new page,  1 line, 2 rows
;  for i=0, 4 do begin
;               cghistoplot, sub_sample[i,*] , position = pos[i,*]                                      
;                 
;   endfor
;

 

 


;Note that the first element, !P.MULTI (0), is set to zero to cause the next plot to begin a new page. 
;To make four plots per page with two columns and two rows, use the following statement:

;!P.MULTI = [1, 2, 2]
;cgplot, sub_sample[0,*],sub_sample[1,*]
;cgplot, sub_sample[1,*],sub_sample[2,*]



; Clear the multiple plot variable.
!P.Multi = 0
  
   
;plot, output_data[*,0]
cgPS_Close


cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
image_name = strcompress(dir + '/' + 'triangle_'  + '.png' , /remove_all)
STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
spawn, strinG


stop 

end