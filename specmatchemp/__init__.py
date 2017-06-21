import os

HOMEDIR = os.environ['HOME']
SPECMATCHDIR = "{0}/.specmatchemp/".format(HOMEDIR)
LIBPATH = "{0}/.specmatchemp/library.h5".format(HOMEDIR)

SPECMATCH_VERSION = 'v0.3'

SHIFT_REFERENCES = [['nso', 'NSO', 5777, None],
                    ['j72.718', '123239', 4800, 'nso'],
                    ['j26.532', '222368', 6200, 'nso'],
                    ['j59.1926', '216899', 3700, 'j72.718']]
