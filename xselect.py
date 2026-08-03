"""
Python wrapper around xselect
"""
import subprocess
import time
import os

class Xselect(object):
    def __init__(self, session_name=None, mission=None):
        self.session_name = 'xsel%d' % int(time.time() % 604800) if session_name is None else session_name
        self.mission = mission
        self.xsel = None

    def start(self):
        args = ['xselect',
                'prefix=%s' % self.session_name]

        self.xsel = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=None, universal_newlines=True,
                                     text=True, bufsize=1)

    def __enter__(self):
        self.start()
        if self.mission is not None:
            self.command('set mission %s' % self.mission)
        return self

    def command(self, cmd):
        self.xsel.stdin.write(cmd + '\n')
        self.xsel.stdin.flush()

    def __exit__(self, exc_type, exc_value, traceback):
        self.command('exit')
        self.command('no')
        self.xsel.stdin.close()
        self.xsel.wait()

    def read_event(self, evl_file):
        evl_dir = os.path.dirname(evl_file)
        if evl_dir == '':
            evl_dir = './'
        evl_filename = os.path.basename(evl_file)

        self.command('read event')
        self.command(evl_dir)
        self.command(evl_filename)

    def filter_time_scc(self, tstart, tend):
        self.command('filter time scc')
        self.command('%g, %g' % (tstart, tend))
        self.command('x')

    def filter_gti(self, gti_file):
        self.command('filter time file %s' % gti_file)
