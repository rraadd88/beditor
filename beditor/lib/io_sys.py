import subprocess
import sys

def runbashcmd(cmd,test=False,logf=None):
    from beditor.lib.global_vars import dirs2ps 
    cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
    if test:
        print(cmd)
    err=subprocess.call(cmd,shell=True,stdout=logf,stderr=subprocess.STDOUT)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)

def is_interactive():
    # thanks to https://stackoverflow.com/a/22424821/3521099
    import __main__ as main
    return not hasattr(main, '__file__')
