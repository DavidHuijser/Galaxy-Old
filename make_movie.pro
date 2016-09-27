
pro make_movie

    cd ,'/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy'
  dir = '/home/dhui890/Documents/Research/Brendon/DNest4/code/Examples/SingleGalaxy'

  spawn, 'ls multi_img_?.png > list.txt'
  spawn, 'ls multi_img_??.png >> list.txt'
  spawn, 'ls multi_img_???.png >> list.txt'
  spawn, 'ls multi_img_????.png >> list.txt'
  
  string = 'mencoder  "mf://@list.txt" -mf fps=8 -o RJ_movie.avi -ovc lavc  -lavcopts vcodec=msmpeg4v2:vbitrate=800 '
spawn,  string

stop 

spawn, 'ls multi_img_*.png > list.txt'
spawn, 'ls multi_img_*.png >> list.txt'
spawn, 'ls multi_img_*.png >> list.txt'

string = 'mencoder  "mf://@list.txt" -mf fps=1 -o test.avi -ovc lavc  -lavcopts vcodec=msmpeg4v2:vbitrate=400 '
;string = 'mencoder  "mf://@list.txt" -mf fps=2 -o test.avi -ovc lavc  -xy 320  -lavcopts vcodec=msmpeg4v2:vbitrate=1800 '
spawn,  string
;  imdisp, observed, /axis



;imdisp, source_image, /axis
;
stop 
end