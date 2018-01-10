os, shutil

root, dirs, filenames in os.walk("/Users/ljf26/OscLib"):
   for f in filenames: 
       print f
       shutil.copyfile(f,f+"_temp")
