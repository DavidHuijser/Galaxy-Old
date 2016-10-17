pro test_DNest4


  DEFSYSV, '!FALSE', 0
  DEFSYSV, '!TRUE', 1
  
  
  ; load data
  cd ,'/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy/Data'
  dir = '/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy/Data'
  info_filename = dir + '/'   + 'fitim.fits'
  
  fits_image = readfits(info_filename)
  fits_image = transpose(fits_image)
  
  
  info_filename = dir + '/'   + 'psfim.fits'
  fits_psf = readfits(info_filename)

  info_filename = dir + '/'   + 'sigma.fits'
  fits_sigma = readfits(info_filename)

  
  cd ,'/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy/'
  dir = '/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy'
  
;  
;   info_filename = dir + '/'   + 'file_2.txt'
;  temp_data = READ_ASCII(info_filename)
;  file_data = temp_data.field001
;  file_data =  transpose(file_data)
;  
;  print,total(file_data)
;  
;  info_filename = dir + '/'   + 'file_5.txt'
;  temp_data = READ_ASCII(info_filename)
;  file_data = temp_data.field001
;  file_data =  transpose(file_data)
;print,total(file_data)
  
;  
;  temp = readfits(dir+'temp.fits')
;  
;  psf = readfits(dir+'psfim.fits')
;
;  info_filename = dir + '/'   + 'Output_results4.txt'
;  temp_data = READ_ASCII(info_filename)
;  output_data = transpose(temp_data.field01)
;  ; output_data =  transpose(output_data )
;  
;  info_filename = dir + '/'   + 'neW_image.txt'
;  temp_data = READ_ASCII(info_filename)
;  new = temp_data.field001
;  new =  transpose(new)
;
;  info_filename = dir + '/'   + 'Final.txt'
;  temp_data = READ_ASCII(info_filename)
;  final = temp_data.field001
;  final =  transpose(final)
;
;  info_filename = dir + '/'   + 'Final_blurred.txt'
;  temp_data = READ_ASCII(info_filename)
;  blur = temp_data.field001
;  blur =  transpose(blur)
;
; 
;;  ; read in test from aarons code
;;  info_filename = dir + '/'   + 'aarons.txt'
;;  temp_data = READ_ASCII(info_filename)
;;  aarons = temp_data.field001
;;  aarons =  transpose(aarons)
;
;





;
;  chi_square = totAl((temp - conv_data)^2); 
;  
;
;
;   print,chi_square
;   helpme = size(temp)
;   metadata = [helpme[1], helpme[2]]
;   helpme = size(output_data)
;   sample_size = helpme[1]
;   ; recreate image
;   
;    
   start = 1
   sample_size = 234
   
   for index = start, sample_size-1 do begin
   
   
     info_filename = strcompress(dir + '/' + 'file_' + string(index) + '.txt' , /remove_all)
     temp_data = READ_ASCII(info_filename)
     final = transpose(temp_data.field001)
     c_residual = (final - fits_image )
     
     image_name1 = strcompress(dir + '/' + 'multi_img_' + string(index) + '.eps' , /remove_all)
     
     !Y.OMargin = [4.5, -0.5]  ; onder , boven
     !X.OMargin = [10, 2]  ; links, rechts
     !P.Multi=[0,3,1]
     ; !P.Font = 0 
     ;!p.position=[0.1, 0.1, 0.95, 0.95]        ; font things
     cgPS_Open, image_name1, /ENCAPSULATED, TT_FONT='Times Italic'
     cgDisplay, 1200, 600
     
     ;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
     loadct, 33
     ;pos = cgLayout([2,2], OXMargin=[3, 3], OYMargin=[3, 14], XGap=5, YGap=12)
     ;Position=pos[*,0]
    
     ; original image
          
     ; reconstructed image 
     cgImage, final, margin=0.7,  Background='white'
     cgText,  0.1,  -0.15, Alignment=0.5, 'model-image', Charsize=cgDefCharsize()*1.0
     
     str =  'max vale = ' + string(max(final))
     cgText,  0.1,  -0.25, Alignment=0.5,str, Charsize=cgDefCharsize()*1.0
     
;     cgImage, idl_residual, margin=0.7,  Background='white'
;     
     
     cgIMAGE,fits_image , margin=0.7,  Background='white' ; , XRange=[metadata[2],metadata[3]],yrange=[metadata[4],metadata[5]], /Axes, margin=0.5, Background='white'
     cgText,    3.7,  -0.45, Alignment=0.5, textoidl('Original Image' ), Charsize=cgDefCharsize()*1.0
     str =  'max vale = ' + string(max(fits_image))
     cgText,  3.7,  -0.85, Alignment=0.5,str, Charsize=cgDefCharsize()*1.0

     
     
;     cgIMAGE,final, margin=0.7,  Background='white'
;     cgText,  3.7,  -0.32, Alignment=0.5, 'C++-image', Charsize=cgDefCharsize()*1.2
;     ;   cgText,    0.0,  8.2, Alignment=0.5, textoidl('Model' ), Charsize=cgDefCharsize()*1.2
     
     cgIMAGE, c_residual, margin=0.7,  Background='white'
     cgText,  3.7,  -0.85, Alignment=0.5, 'Residual', Charsize=cgDefCharsize()*1.0
     ;   cgText,    0.0,  8.2, Alignment=0.5, textoidl('Model' ), Charsize=cgDefCharsize()*1.2
     
     cgPS_Close
     
     
     ; Create a PNG file with a width of 600 pixels.
     
     cgPS2Raster, image_name1, /PNG, Width=600, /DELETE_PS
     
 ;    image_name2 = strcompress(dir + '/' + 'multi_img_' + string(index+1) + '.png' , /remove_all)
    ; STRING =  'convert '+ image_name1  + image_name2
     ;STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
     ;image_name2 = strcompress(dir + '/' + 'multi_img_' + string(l+1) + '.png' , /remove_all)
     ;string =  'convert -density 300 '+ image_name + ' -rotate -0  ' + image_name2
 ;    spawn, strinG
     
     ;    spawn, 'rm *.eps'
     ;cgPS2Raster, image_name, /DELETE_PS
     !P.Multi = 0
     !p.font=-1
     !p.position=[0.0,0.0,0.0,0.0]
     !Y.OMargin = [0, 0]
     !X.OMargin = [0, 0]
     
     
     
   endfor
   
   stop 
   
 ; STRING =  'convert '+ image_name +' -rotate -90 ' + image_name   stop 
   
   
   zz= convolve(img, psf )
   r_e = parameters[5]  
   I_e = parameters[4]
   
   sersic_index = parameters[6];    // sersic-index
   q =  parameters[7];    // Axis Ratio (b/a)
   position_angle = parameters[8];   // position angle
   ;k = 2*sersic_index - 0.324  
   ;intflux =  2*!Pi*
  ; print, "integrated flux",intflux  
   

   print,"Log likelihood: " ,  -0.5*total((temp - new)^2) ,  output_data[sample_size-1, 0]
   
   
   
   stop 
 
  z = maxlocation(img)
  zz =  maxlocation(new)
  print, "      IDL-img      C++ IMG "  
  print ,z ,"           " , zz
  print ,total(img) ,"           " , total(new), "differece = ", total(abs(new - img))
  
  

  
;  
;  Maxi =  1136.88
;  Mini = -67.6093
;  magzp = 30
;
;  
;  
;  
;  print, (alog10(Maxi) - 0.4*30) / (-0.4)
;  print, (alog10(Mini) - 0.4*30) / (-0.4)
  

stop 


;tab = READFITS( dir+'fitim.fits', header)    ; gives an array of 157 x 171  with a total flux of   1.07620e+06
                                     ; 0x0 is the lower left corner,  157 wide and 171 high.
                                     ; 
                                     ;
                                     ;  
                                     
                                     
;psf = READFITS( dir+'psfim.fits')
;pstamp = READFITS( dir+'pstamp.fits')


zeropoint = 30;
exp_time = 1
psf_m_tot  = 16.0509;              integrate  magnitude

f_tot = exp_time*10^( (psf_m_tot - zeropoint)/ (-2.5))



;




min_x = 0
max_x = 25 
min_y =  133.38 -12.5
max_y =   133.38 + 12.5
;print, total(tab[min_x:max_x-1, min_y:max_y-1]   - psf*362600.01)
;sub = tab[min_x:max_x-1, min_y:max_y-1]
;tab[min_x:max_x-1, min_y:max_y-1]   = tab[min_x:max_x-1, min_y:max_y-1] - psf*332332
;tab[min_x:max_x-1, min_y:max_y-1]   = tab[min_x:max_x-1, min_y:max_y-1] - psf*379875
                                      
 x = size(taB)
 width = x[1]
 height = x[2]
 print, width, height                                    
 
totals = make_array(height, /double)
for i = 0L, height-1 do begin
  totals[i] = total(tab[*,i])  
endfor

;;original image read from FITS in tab
;;'altered image read from txt in tab2 


info_filename = dir + '/'   + 'Converted.txt'
temp_data = READ_ASCII(info_filename)
conv_data = temp_data.field001
conv_data =  transpose(conv_data)

info_filename = dir + '/'   + 'Circ.txt'
temp_data = READ_ASCII(info_filename)
circ_data = temp_data.field0001
circ_data =  transpose(circ_data)


;temp = readfits(dir+'temp.fits')
stop 

info_filename = dir + '/'   + 'Galaxy.txt'
temp_data = READ_ASCII(info_filename)
info_data = temp_data.field001
info_data =  transpose(info_data)

info_filename = dir + '/'   + 'TryOutImage2.txt'
temp_data2 = READ_ASCII(info_filename)
info_data2 = temp_data2.field001
info_data2 =  transpose(info_data2)

print, "original:", total(info_data)
print, "polar:", total(circ_data)
print, "converted polar:", total(conv_data)

print, total( tab[20:150, *] ), total( info_data[20:150,*]),  total( conv_data[20:150,*])


THX = indgen(130) + 20

;mdisp, tab[20:150, *] -  info_data[20:150,*]*2.24



;temp = readfits(dir+'ps.fits')

temp = readfits(dir+'fitim.fits',header)
new_tab = temp[50:150,25:125]
writefits,dir+'temp2.fits',new_tab, header


;temp = readfits(dir+'fitim.fits',header)




stop 



zeropoint = 30;
exp_time = 1
m_tot  = 16.0509;              integrate  magnitude
r_e = 22.2706 
n = 2                      ;sersic index  
kappa = 2*n - 0.331 
q = 0.8390
f_tot = exp_time*10^( (m_tot - zeropoint)/ (-2.5))

factors =   2.0*!pi*r_e^2*exp(kappa)*n*kappa^(-2*n)*gamma(2*n)*q
sigma_e = f_tot/ factors

print, sigma_e       ; 55.9523  for no-boxyness

c = 0
boxy = !pi*c/(4*beta(1/c,1+1/c) )
sigma_e = f_tot*boxy/ factors
print, sigma_e       ; 55.9523  for no-boxyness


stop
  end
  
  
 function create_new_model, parameters, metadata
   xsize = metadata[0]
   ysize = metadata[1]
   
   new_image =  make_Array(xsize, ysize, value=0.0, /double)
   
   x_pos = round(parameters[2])
   y_pos = round(parameters[3])
   I_e = parameters[4];    // Ie
   R_e = parameters[5];    // Re
   sersic_index = parameters[6];    // sersic-index
   q =  parameters[7];    // Axis Ratio (b/a)
   position_angle = parameters[8];   // position angle
   
    c= parameters[9]
    b_m= 2.0*sersic_index - 0.324;
    
    
   cos_theta = cos(position_angle)
   sin_Theta = sin(position_angle)
   for xindex = 0L, xsize-1 do begin
         for yindex = 0L, ysize-1 do begin
                 
                 x_prime =  (xindex - x_pos)*cos_theta + (yindex - y_pos)*sin_theta;
                 y_prime =  -(xindex - x_pos)*sin_theta + (yindex - y_pos)*cos_theta;               
                ; R =    sqrt(     x_prime^2.0 +  (y_prime/q)^(2.0)    ) + 0.00000001;
                  
                 R =    (    abs( x_prime)^(c+2.0) +  abs((y_prime)/q)^(c+2.0)    )^ (1.0/(c+2.0) );
                                 
                 e_factor =  q/sqrt(  q*q*cos_theta^2 + sin_theta^2) ;
                ; e_factor =1.000;       
                 factor =  (R/(e_factor*R_e))^(1.0/sersic_index) -1 ;
                 if (-b_m*factor gt  -10 ) then begin
                    new_image[xindex,yindex] =  I_e*exp(-b_m*factor)
                 endif
                                                                       
                   
         endfor    
    endfor
   return, new_image
    
 end
 
 function maxlocation, img
   helpme = size(img)
   index = where(max(img) eq img)
   
   xsize = helpme[1]  
   ysize = helpme[2]
   if n_elements(index eq 1) and (index ne -1) then begin
       xpos = index mod xsize 
       ypos = index / xsize
   endif
            
   return,  [xpos, ypos]
 
 end
 
 
