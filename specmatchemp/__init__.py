import os

HOMEDIR = os.environ['HOME']
SPECMATCHDIR = "{0}/.specmatchemp/".format(HOMEDIR)
LIBPATH = "{0}/.specmatchemp/library.h5".format(HOMEDIR)

SPECMATCH_VERSION = 'v0.2'

SHIFT_REFERENCES = [['nso', 'NSO', 5777, None],
                    ['rj72.718', '123239', 4800, 'nso'],
                    ['rj26.532', '222368', 6200, 'nso'],
                    ['rj59.1926', '216899', 3700, 'rj72.718']]
