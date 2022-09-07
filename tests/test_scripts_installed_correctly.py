import glob
import os
import subprocess as sp


def test_scripts_installed():
    cwd = os.path.dirname(__file__)

    scripts = glob.glob(os.path.join(cwd, '..', 'bin', '*.py'))

    for script in scripts:
        if '__init__.py' in script:
            continue
        name = os.path.basename(script).replace('.py', '')
        print('Testing', name)
        sp.check_call('{} -h'.format(name).split(' '))
