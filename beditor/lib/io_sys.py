import subprocess
import sys

def runbashcmd(cmd,test=False):
    from beditor.lib.global_vars import dirs2ps 
    cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
    if test:
        print(cmd)
    err=subprocess.call(cmd,shell=True)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)
