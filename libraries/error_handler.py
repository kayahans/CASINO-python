import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    UNDERLINE = '\033[4m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def error(message):
    print (bcolors.FAIL + bcolors.BOLD + str("ERROR: " + message) + bcolors.ENDC)
    sys.exit()

def warning(message):
    print (bcolors.WARNING + bcolors.BOLD + str("WARNING: " + message) + bcolors.ENDC)