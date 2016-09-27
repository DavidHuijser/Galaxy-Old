; this work with Test_Sampler4, where sigma is added to the parameters 
pro Make_mock_galaxies

;ash cd   temp = readfits(dir+'fitim.fits',header)
;  new_tab = temp[50:150,25:125]
;  writefits,dir+'temp2.fits',new_tab, header

 abs_max_objects  = 10 
 abs_max_parameters = 15
 

  DEFSYSV, '!FALSE', 0
  DEFSYSV, '!TRUE', 1
  
  !P.Multi = 0
  
  cd , '/home/dhui890/Documents/Research/C/RJ_SS'
  dir = '/home/dhui890/Documents/Research/C/RJ_SS/'


  name = 'Final.txt'
  ;  info_filename = dir + '/'   + 'neW_image.txt'
  temp_data = READ_ASCII(name)
  new = temp_data.field001
  final =  transpose(new)
 
   name = 'mockdata.txt'   
    ;  info_filename = dir + '/'   + 'neW_image.txt'
    temp_data = READ_ASCII(name)
    new = temp_data.field001
    new =  transpose(new)
    
        
    info_filename = dir    + 'MockdataInfo.txt'
    temp_data = read_ascii(info_filename)
    output_data = transpose(temp_data.field001)
    
    sample_size = 1 
    
    object_types      =  make_Array(abs_max_objects, /double)
    object_parameters  =  make_Array(abs_max_objects, /double)
    parameters = make_Array(abs_max_objects,abs_max_parameters, /double)
    Object = {loglike:0.0,logprior:0.0, sigma:0.0, number_of_objects:0,types:object_types, object_parameters:object_parameters, parameters:parameters}
    
    cat = REPLICATE(Object, sample_size)
    for index = 0L, sample_size-1 do begin
      cat[index].loglike = output_data[index, 0]
      cat[index].logprior = output_data[index,1]
      cat[index].sigma = output_data[index,2]
      cat[index].number_of_objects = output_data[index,3]
      cat[index].types = output_data[index,4:13]
      cat[index].object_parameters = output_data[index,14:23]
      cat[index].parameters = transpose(reform(output_data[index,24:173],15,abs_max_objects))
      print, cat[index].loglike, cat[index].logprior,  cat[index].sigma,  cat[index].number_of_objects
     
    endfor
    
    noise = randomn(seed, 100,100)
    zz = new + 1.05*noise

    
    
    image_name = strcompress(dir + '/' + 'mock_data.eps' , /remove_all)
    cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
    ; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
    cgimage,  zz
    cgPS_Close
    
    ; Create a PNG file with a width of 600 pixels.
    cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
    image_name = strcompress(dir + '/' + 'mock_data.png' , /remove_all)
    STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
    spawn, strinG
    
    

    
    

    mock_filename = strcompress('mock_data.fits', /remove_all)
    temp = readfits(dir+'temp2.fits',header)
    writefits,dir+mock_filename,zz, header
    temp = readfits(dir+mock_filename,header)
    
    
   stop 
   

 ;  read in the original image  
  temp = readfits(dir+'temp2.fits')
 
  ; read the PSF 
  psf = readfits(dir+'psfim.fits')

  temp_array = size(output_data ) 
  sample_size = 1 
  
  
  add = '_four_objects_'
  

  object_types      =  make_Array(abs_max_objects, /double)
  object_parameters  =  make_Array(abs_max_objects, /double)
  ;parameters = make_Array(abs_max_objects,abs_max_parameters, /double)
  parameters = transpose(reform(make_array(15*10, value=0, /double),15,10))
  Object = {loglike:0.0,logprior:0.0, sigma:0.0, number_of_objects:0,types:object_types, object_parameters:object_parameters, parameters:parameters}

  cat = REPLICATE(Object, sample_size)
  
   i = 0
   
   
   
   metadata = [100, 100]
    cat[0].sigma =  2 
    cat[0].number_of_objects = 4
    cat[0].types = [0,0,  0,0,  0,0,  0,0,  0,0]
    cat[0].object_parameters = [8,8,  8,8,  0,0,  0,0,  0,0]
    
    ; help, cat.parameters <Expression>    DOUBLE    = Array[10, 15, 1000]IDL>

   seed = !NULL   
 ;  cat[0].parameters[0,*] = [ 79.39-51,    85.49-20,  16.05,  22.27,  1.50, 0.84,  -0.5    ,2*randomu(seed, 1)-1, 0,0,0 ,0 ,0 ,0 ,0]
 ;  cat[0].parameters[1,*] = [109.02 -51,    36.41-20, 19.73 ,  14.79,  0.8, 0.51,   1.0787 ,2*randomu(seed, 1)-1, 0,0,0 ,0 ,0 ,0 ,0]  
 ;  cat[0].parameters[2,*] = [ 55.71-51,    65.14-20,  19.52,  15.03,  0.7, 0.72, -0.866207,2*randomu(seed, 1)-1, 0,0,0 ,0 ,0 ,0 ,0]
 ;  cat[0].parameters[3,*] = [100.87-51,    76.00-20,  19.03,   15.19, 0.82, 0.50,  1.06692 ,2*randomu(seed, 1)-1, 0,0,0 ,0 ,0 ,0 ,0]
  
  
  
  
 ;   cat[0].parameters[0,*] = [ 55.39,    85.49,  17.05,  33.27,  1.50, 0.84,  -0.5    ,1*randomu(seed, 1), 0,0,0 ,0 ,0 ,0 ,0]
  
   cat[0].parameters[0,*] = [ 79.39-51,    85.49-20,  17.05,  23.27,  1.50, 0.84,  -0.5    ,1*randomu(seed, 1), 0,0,0 ,0 ,0 ,0 ,0]
   cat[0].parameters[1,*] = [109.02 -51,    36.41-20, 20.73 ,  13.79,  0.8, 0.51,   1.0787 ,1*randomu(seed, 1), 0,0,0 ,0 ,0 ,0 ,0]
   cat[0].parameters[2,*] = [ 55.71-51,    55.14-20,  20.52,  14.03,  0.7, 0.72, -0.866207,1*randomu(seed, 1), 0,0,0 ,0 ,0 ,0 ,0]
   cat[0].parameters[3,*] = [109.87-51,    76.00-20,  20.03,   14.19, 0.82, 0.50,  8.06692 ,1*randomu(seed, 1), 0,0,0 ,0 ,0 ,0 ,0]
  
  ; 8.0000000       8.0000000       8.0000000       8.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 28.389999       58.019997       4.7099991       49.870003       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 65.489998       16.410000       45.139999       56.000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 16.049999       19.730000       19.520000       19.030001       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 22.270000       14.790000       15.030000       15.190000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 1.5000000      0.80000001      0.69999999      0.81999999       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 0.83999997      0.50999999      0.72000003      0.50000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; -0.50000000       1.0786999     -0.86620700       1.0669200       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000       0.0000000
  ; 0.35162425       1.3140719     -0.92474151     -0.38943553
  
  

;    for index = 0L, 4-1 do begin
;        cat[0].parameters[index, 0] =  randomu(seed, 1)*100
;        cat[0].parameters[index, 1] = randomu(seed, 1)*100
;        cat[0].parameters[index, 2] = randomu(seed, 1)*100
;        cat[0].parameters[index, 3] = randomu(seed, 1)*sqrt(100^2 + 100^2) 
;        cat[0].parameters[index, 4] = randomu(seed, 1)*9.5 + 0.5
;        cat[0].parameters[index, 5] = randomu(seed, 1)
;        cat[0].parameters[index, 6] = randomu(seed, 1)*!pi - 0.5*!pi
;        cat[0].parameters[index, 7] = randomu(seed, 1)
;        
;    endfor
    
    img = create_new_model( cat, metadata)
    zz= convolve(img, psf )
    noise = randomn(seed, 100,100)
    
    zz = zz + noise*1.05*(mean(zz))
    ;         idl_residual = abs(temp - zz)
 
    ; Create a PNG file with a width of 600 pixels.
    
    image_name = strcompress(dir + '/' + 'mock_data_' + add + '.eps' , /remove_all)
    cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
    ; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
    cgimage,zz   
    cgPS_Close
    
    
    ; Create a PNG file with a width of 600 pixels.
    cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
    image_name = strcompress(dir + '/' + 'mock_data_'  +  add +'.png' , /remove_all)
    STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
    spawn, strinG
     
    mock_filename = strcompress('mock_data' + add + '.fits', /remove_all)
    temp = readfits(dir+'temp2.fits',header)
    writefits,dir+mock_filename,zz, header
    temp = readfits(dir+mock_filename,header)
    
    info_filename = strcompress('info_mock_data' + add + '.txt', /remove_all)
    openw, lun2, info_filename, /get_lun    
    for index = 0L, cat[0].number_of_objects -1 do begin
       printf, lun2, transpose(cat[0].parameters[index,*] )
    endfor      
    close, lun2    
    free_lun, lun2
    
    stop
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
;  chi_square = totAl((temp - conv_data)^2); 
;  
;
;
;   print,chi_square
   helpme = size(temp)
   metadata = [helpme[1], helpme[2]]
   helpme = size(output_data)
   sample_size = helpme[1]
   ; recreate image 
   
   new_data = make_Array(helpme[1], helpme[2], /double)
   new_data[*,*] = output_data[*,*] 
   ;output_data = double(output_data, helpme[1], helpme[2])   
   
;   

   image_name = strcompress(dir + '/' + 'log_likelihood_'  + '.eps' , /remove_all)
   cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
   ; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
   plot, output_data[*,0]
   cgPS_Close
   
   
   ; Create a PNG file with a width of 600 pixels.
   cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
   image_name = strcompress(dir + '/' + 'log_likelihood_'  + '.png' , /remove_all)
   STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
      spawn, strinG
  
stop
 
   start = 0
 
   for index = start, sample_size-1 do begin
    
         temp = readfits(dir+'temp2.fits')
     
         ; read the PSF
         psf = readfits(dir+'psfim.fits')
     
         parameters = new_data[index,*]
;         img = create_new_model( parameters[0,*], metadata)         
;         zz= convolve(img, psf )
;         idl_residual = abs(temp - zz)
         
         info_filename = strcompress(dir + '/' + 'file_' + string(index+1) + '.txt' , /remove_all)
         temp_data = READ_ASCII(info_filename)
         final = transpose(temp_data.field001)       
         c_residual = temp - final
         
         image_name = strcompress(dir + '/' + 'multi_img_' + string(index+1) + '.eps' , /remove_all)

         !Y.OMargin = [1, 1]
         !X.OMargin = [8, -2]                      
                         
         ;!Y.OMargin = [0.5, -0.5]  ; onder , boven
        ; !X.OMargin = [2, 2]  ; links, rechts
         !P.Multi=[0,3,1]
         ; !P.Font = 0         ; font things
         cgPS_Open, image_name, /ENCAPSULATED, TT_FONT='Times Italic'
         cgDisplay, 1600, 800
         
         ;cgDisplay, 800, 550
         ;cgPlot, cgDemoData(17), Position=[0.10, 0.10, 0.45, 0.90]
         ;cgPlot, cgDemoData(17), Position=[0.55, 0.40, 0.95, 0.90], /NoErase
         ;cgPlot, cgDemoData(17), Position=[0.55, 0.10, 0.95, 0.30], /NoErase
         
         ;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
         loadct, 33
         ;pos = cgLayout([2,2], OXMargin=[3, 3], OYMargin=[3, 14], XGap=5, YGap=12)
         ;Position=pos[*,0]
         ;p = [0.1, 0.1, 0.95, 0.95] ; [links, onder, rechts, boven]
         ; original image 
         
                    
        ; reconstructed image (IDL)             
;         cgImage, zz, margin=0.7,  Background='white'
;         cgText,  0.1,  -0.15, Alignment=0.5, 'IDL-image', Charsize=cgDefCharsize()*1.0
;                , /NoErase
;         cgImage, idl_residual, margin=0.7,  Background='white'
;         cgText,  3.7,  -0.5, Alignment=0.5, 'IDL Residual', Charsize=cgDefCharsize()*1.0
         ; [left,bottom,right,top]
         cgIMAGE,temp , margin=0.5,  Background='white', Position=[0.05, 0.1, 0.95, 0.95] ; , XRange=[metadata[2],metadata[3]],yrange=[metadata[4],metadata[5]], /Axes, margin=0.5, Background='white',
        ; cgText,    3.7,  -0.5, Alignment=0.5, textoidl('Original Image' ), Charsize=cgDefCharsize()*1.0
        
         cgIMAGE,final, margin=0.5,  Background='white', Position=[0.05, 0.1, 0.95, 0.95]
        ; cgText,  3.7,  -0.32, Alignment=0.5, 'C++-image', Charsize=cgDefCharsize()*1.2
       ;   cgText,    0.0,  8.2, Alignment=0.5, textoidl('Model' ), Charsize=cgDefCharsize()*1.2
                 
         cgIMAGE, c_residual, margin=0.5,  Background='white', Position=[0.05, 0.10, 0.95, 0.95]
        ;  cgText,  3.7,  -0.33, Alignment=0.5, 'C++ Residual', Charsize=cgDefCharsize()*1.2
      ;   cgText,    0.0,  8.2, Alignment=0.5, textoidl('Model' ), Charsize=cgDefCharsize()*1.2
         
         
         cgText,   -15.5,  0.5, Alignment=0.5, textoidl('Original Image' ), Charsize=cgDefCharsize()*1.2
      cgText,  -5.25,  0.5, Alignment=0.5, textoidl('Model Image'), Charsize=cgDefCharsize()*1.2
      cgText,  4.0,  0.5, Alignment=0.5, textoidl('Residual'), Charsize=cgDefCharsize()*1.2
      
      
      cgText,   -15.0,  -0.25, Alignment=0.5, textoidl('Number of Objects:'+string(cat[index].number_of_objects) ), Charsize=cgDefCharsize()*1.2
      
      cgText,   -4.25,  -0.25, Alignment=0.5, textoidl('Log(likelihood):'+string(cat[index].loglike) ), Charsize=cgDefCharsize()*1.2
         cgPS_Close
         
         
         ; Create a PNG file with a width of 600 pixels.
         
         cgPS2Raster, image_name, /PNG, Width=600, /DELETE_PS
         image_name = strcompress(dir + '/' + 'multi_img_' + string(index+1) + '.png' , /remove_all)
;         STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
         ;STRING =  'convert '+ image_name +' -rotate -90 ' + image_name
         ;image_name2 = strcompress(dir + '/' + 'multi_img_' + string(l+1) + '.png' , /remove_all)
         ;string =  'convert -density 300 '+ image_name + ' -rotate -0  ' + image_name2
      ;   spawn, strinG
         
         ;    spawn, 'rm *.eps'
         ;cgPS2Raster, image_name, /DELETE_PS
         !P.Multi = 0
         !p.font=-1
         !p.position=[0.0,0.0,0.0,0.0]
         !Y.OMargin = [0, 0]
         !X.OMargin = [0, 0]
        
   
        
   endfor
   
   
   ;prints the parameters of the first component of the 99th sample
   print, cat[99].parameters[0,*]  


   ;prints all the x-parameters of the 99th sample
   print, cat[99].parameters[*,1]

   stop 
   
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

new_tab = tab[20:151,*]
writefits,dir+'temp.fits',new_tab, header


temp = readfits(dir+'temp.fits')


  end
  
  ;object_types      =  make_Array(abs_max_objects, /double)
  ;object_parameters  =  make_Array(abs_max_objects, /double)
  ;parameters = make_Array(abs_max_objects,abs_max_parameters, /double)
  ;parameters = transpose(reform(make_array(15*10, value=0),15,10))
  ;Object = {loglike:0.0,logprior:0.0, sigma:0.0, number_of_objects:0,types:object_types, object_parameters:object_parameters, parameters:parameters}
  
 function create_new_model, cat, metadata
   xsize = metadata[0]
   ysize = metadata[1]
   
   new_image =  make_Array(xsize, ysize, value=0.0, /double)

   for xindex = 0L, xsize-1 do begin
         for yindex = 0L, ysize-1 do begin
                 
                 for index = 0L, 4 do begin
                  
                   x_pos = round(cat[0].parameters[index,0])
                   y_pos = round(cat[0].parameters[index,1])
                   I_e = cat[0].parameters[index,2];    // Ie
                   R_e = cat[0].parameters[index,3];    // Re
                   sersic_index = cat[0].parameters[index,4];    // sersic-index
                   q =  cat[0].parameters[index,5];    // Axis Ratio (b/a)
                   position_angle = cat[0].parameters[index,6];   // position angle
           
                   c= cat[0].parameters[index,7]
           
                   b_m= 2.0*sersic_index - 0.324;
                   
                   cos_theta = cos(position_angle)
                   sin_theta = sin(position_angle)
                  
                 
                 x_prime =  (xindex - x_pos)*cos_theta + (yindex - y_pos)*sin_theta;
                 y_prime =  -(xindex - x_pos)*sin_theta + (yindex - y_pos)*cos_theta;               
                ; R =    sqrt(     x_prime^2.0 +  (y_prime/q)^(2.0)    ) + 0.00000001;
                  
                 R =    (    abs( x_prime)^(c+2.0) +  abs((y_prime)/q)^(c+2.0)    )^(1.0/(c+2.0) );
                                 
                 e_factor =  q/sqrt(  q*q*cos_theta^2 + sin_theta^2) ;
                ; e_factor =1.000;       
                 factor =  (R/(e_factor*R_e))^(1.0/sersic_index) -1 ;
                 if (-b_m*factor gt  -10 ) then begin
                    new_image[xindex,yindex] = new_image[xindex,yindex] +  I_e*exp(-b_m*factor)
                 endif
                                                                       
                endfor
   
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
 
 
