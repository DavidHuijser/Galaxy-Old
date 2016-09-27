
;function READFITS, filename, header, heap, CHECKSUM=checksum, $
 ; COMPRESS = compress, HBUFFER=hbuf, EXTEN_NO = exten_no, $
 ; NOSCALE = noscale, NSLICE = nslice, $
 ; NO_UNSIGNED = no_unsigned,  NUMROW = numrow, $
 ; POINTLUN = pointlun, SILENT = silent, STARTROW = startrow, $
 ; NaNvalue = NaNvalue, FPACK = fpack, UNIXpipe=unixpipe



function cornerplot, filename, data, dim

 
 cgPS_Open, filename, /ENCAPSULATED, TT_FONT='Times Italic'
  ; plot, output_data[*,0], yrange = [-5.1744e6, -5.1742e6]
  cgDisplay, 1200, 1200
  
  ;cgDisplay, 2000, 2000       ; Load the color table for the display. All zero values will be gray.
  loadct, 33
 
  pos = cgLayout([dim ,dim ])
  FOR j=0,dim -1  DO BEGIN
    for i= 0L, j do begin
      if (i eq  j ) then begin
        cghistoplot, sub_sample[i,*], NoErase=j NE 0, Position=pos[*,i+dim *j], Title='Plot ' + StrTrim(j+1,2)
      endif else begin
        cgplot, sub_sample[i,*], sub_sample[j,*], NoErase=j NE 0, Position=pos[*,i+dim *j], Title='Plot ' + StrTrim(j+1,2), psym=2
      endelse
    endfor
    
    
  ENDFOR


  
  ; Clear the multiple plot variable.
  !P.Multi = 0
  
  
  ;plot, output_data[*,0]
  cgPS_Close
  cgPS2Raster, filename, /PNG, Width=600, /DELETE_PS
  filename = strcompress('triangle.png' , /remove_all)
  STRING =  'convert '+ filename +' -rotate -90 ' + filename
  spawn, strinG


end
