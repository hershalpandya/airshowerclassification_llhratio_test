
import os
import glob

syntax0='svn up'
os.system(syntax0)

nblist = glob.glob('*.ipynb')

svnlist=os.popen('svn list').read()
for nb in nblist:
    #syntax0='svn delete --keep-local %s'%nb
    #print syntax0
    #os.system(syntax0)
    
    pythonname=nb[:-6]+'.py'
    if ' ' in pythonname:
        pythonname=('_').join(pythonname.split(' '))

    syntax1='jupyter nbconvert --to python %s --output %s'%(nb,pythonname)
    print syntax1
    os.system(syntax1)

    
    insvnlist=pythonname in svnlist

    if not insvnlist:
        syntax2='svn add %s'%pythonname
        print syntax2
        os.system(syntax2)

syntax3='svn ci -m "updating python notebooks"'
os.system(syntax3)
