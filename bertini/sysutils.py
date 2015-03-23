from __future__ import print_function

from os import chdir
from os.path import dirname
from subprocess import check_output, CalledProcessError

from naglib.startup import TEMPDIR
from naglib.core.base import NAGobject
from naglib.exceptions import BertiniError, NoBertiniException

def __os():
    from sys import platform

    if platform.startswith('win'):
        return 'WINDOWS'
    elif platform.startswith('cygwin'):
        return 'CYGWIN'
    elif platform.startswith('linux'):
        return 'LINUX'
    elif platform.startswith('darwin'):
        return 'OSX'

def __has_bertini():
    platform = __os()
    if platform == 'WINDOWS':
        cmd = 'where.exe'
    else:
        cmd = 'which'

    try:
        bertinipath = check_output([cmd, 'bertini'])
    except CalledProcessError:
        bertinipath = ''

    return bertinipath.strip()

def __proc_err_output(output):
    lines = output.split('\n')
    dex = -1
    for l in lines:
        if l.startswith('ERROR'):
            dex = lines.index(l)
    if dex == -1:
        return output
    else:
        l = lines[dex]
        # strip 'ERROR: '
        lines[dex] = l[l.index(' ')+1:]
        # strip 'Bertini will now exit due to this error'
        return '\n'.join(lines[dex:-1])

BERTINI = __has_bertini()

class BertiniRun(NAGobject):
    #TTEVALSYS = -4
    #TTEVALSYSJAC = -3
    #TTNEWTON = -2
    #TTNEWTJAC = -1
    #TTSOLVE = 0
    #TTNID = 1
    #TTSAMPLE = 2
    #TTMEMBER = 3
    #TTPRINT = 4
    #TTPROJECTION = 5
    #TTSTAB = 6
    #TTREGEN = 7
    def __init__(self, system, start=(), startp=(), finalp=(), config={}, stdin=None):   
        from tempfile import mkdtemp
        self._dirname = mkdtemp(prefix=TEMPDIR)
        self._bertini = BERTINI
        self._system = system
        if start and type(start) not in (tuple, list):
            start = [start]
        self._start = start
        if startp and type(startp) not in (tuple, list):
            startp = [startp]
        self._startparams = startp
        if finalp and type(finalp) not in (tuple, list):
            finalp = [finalp]
        self._finalparams = finalp
        if config and not hasattr(config, 'keys'):
            msg = 'cannot understand keyword arguments'
            raise TypeError(msg)
        self._config = config
        self._stdin = stdin
        self._complete = False
        self.__inputf = []
        
    def __proc_err_output(self, output):
        lines = output.split('\n')
        dex = -1
        for l in lines:
            if l.startswith('ERROR'):
                dex = lines.index(l)
        if dex == -1:
            return output
        else:
            l = lines[dex]
            # strip 'ERROR: '
            lines[dex] = l[l.index(' ')+1:]
            # strip 'Bertini will now exit due to this error'
            return '\n'.join(lines[dex:-1])
        
    def __recover_input(self):
        """
        Reads main_data and recovers the input file
        needed to reproduce a run
        """
        filename = self._dirname + "/main_data"
        key = "*************** input file needed to reproduce this run ***************\n"
        with open(filename, 'r') as fh:
            lines = fh.readlines()
        fh.close()
        try:
            dex = lines.index(key)
        except ValueError:
            msg = 'no main_data file!'
            raise BertiniError(msg)
        inlines = lines[dex+1:]
        
        return inlines
        
    def __write(self, ifile=False):
        from os.path import exists
        from re import sub as resub
        
        dirname = self._dirname
        if not exists(dirname):
            from os import mkdir
            mkdir(dirname)
        
        ### write the system
        system = self._system
        config = self._config
        filename = dirname + '/input'
    
        polynomials  = [p.factor() for p in system.polynomials]
        variables    = system.variables
        parameters   = system.parameters
        homvar       = system.homvar
        num_polys    = system.shape[0]
        
        options = config.keys()
        
        str_poly = [str(p) for p in polynomials]
        str_poly = [resub(string=p, pattern=r'\*\*', repl='^') for p in str_poly]
        str_vars = [str(v) for v in variables]
        str_pars = [str(p) for p in parameters]
        
        poly_names = ['f{0}'.format(i+1) for i in range(num_polys)]
        polys_named = zip(poly_names, str_poly)
        
        poly_list = ','.join([f for f in poly_names])
        vars_list = ','.join([v for v in str_vars])
        pars_list = ','.join([p for p in str_pars])

        fh = open(filename, 'w')
        
        # write the CONFIG section
        print('CONFIG', file=fh)
        for option in options:
            print('{0}:{1};'.format(option, config[option]), file=fh)
        print('END', file=fh)
        
        # write the INPUT section
        print('INPUT', file=fh)
        if parameters:
            print('parameter {0};'.format(pars_list), file=fh)
        if homvar:
            print('hom_variable_group {0};'.format(vars_list), file=fh)
        else:
            print('variable_group {0};'.format(vars_list), file=fh)
        print('function {0};'.format(poly_list), file=fh)
        
        for p in polys_named:
            # p is a key-value pair, e.g., ('f1', 'x^2 - 1')
            print('{0} = {1};'.format(p[0], p[1]), file=fh)
        print('END', file=fh)
        
        # finish up
        fh.close()
        
        ### write the start file, if there's a start point
        start = self._start
        if start:
            from naglib.bertini.fileutils import fprint
            startfile = dirname + '/start'
            fprint(start, startfile)
        else:
            startfile = ''
        
        ### write the start_parameters file, if there are start parameters
        startparams = self._startparams
        if startparams:
            from naglib.bertini.fileutils import fprint
            startpfile = dirname + '/start_parameters'
            fprint(startparams, startpfile)
       
       ### write the final_parameters file, if there are final parameters
        finalparams = self._finalparams
        if finalparams:
            from naglib.bertini.fileutils import fprint
            finalpfile = dirname + '/final_parameters'
            fprint(finalparams, finalpfile)
        
        return filename, startfile
    
    def rerun(self):
        if not self._complete:
            return self.run()
        else:
            ifile = self.__inputf
            fh = open(self._dirname + '/input', 'w')
            for line in ifile:
                fh.write(line)
            fh.close()
            return self.run()
                
    def run(self):
        cmd = self._bertini
        if not cmd:
            raise(NoBertiniException)
        stdin = self._stdin
        dirname = self._dirname
        input_file,startfile = self.__write()
        
        if not startfile:
            arg = [cmd, input_file]
        elif startfile:
            arg = [cmd, input_file, startfile]
            
        chdir(dirname)
        if stdin:
            stdin = open(stdin, 'r')
        try:
            output = check_output(arg, stdin=stdin)
        except CalledProcessError as e:
            msg = self.__proc_err_output(e.output)
            raise BertiniError(msg)

        if stdin:
            stdin.close()
            
        self._complete = True
        self.__inputf = self.__recover_input()
        return output
        
    @property
    def bertini(self):
        return self._bertini
    @bertini.setter
    def bertini(self, bert):
        self._bertini = bert
    @property
    def complete(self):
        return self._complete
    @property
    def dirname(self):
        return self._dirname
    @property
    def inputf(self):
        return self.__inputf

def call_bertini(input_file, start_file='', cmd=BERTINI, stdin=None):
    """
    Call Bertini
    
    Keyword arguments:
    input_file -- string, the path of a Bertini input file
    start_file -- optional string, the path of a Bertini start file
    cmd        -- optional string, the path of the Bertini executable
    redirect   -- optional string, path to stdin redirect
    """
    if not cmd:
        raise(NoBertiniException)
    if not start_file:
        arg = [cmd, input_file]
    elif start_file:
        arg = [cmd, input_file, start_file]
        
    wd = dirname(input_file)
    if wd:
        chdir(wd)
    if stdin:
        stdin = open(stdin, 'r')
    try:
        output = check_output(arg, stdin=stdin)
    except CalledProcessError as e:
        raise(BertiniError(__proc_err_output(e.output)))

    if stdin:
        stdin.close()
    return output