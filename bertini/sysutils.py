from __future__ import print_function

from os import chdir
from os.path import dirname
from subprocess import check_output, CalledProcessError

from naglib.startup import TOL, TEMPDIR as basedir
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
    TEVALP = -4
    TEVALPJ = -3
    TNEWTP = -2
    TNEWTPJ = -1
    TZERODIM = 0
    TPOSDIM = 1
    TSAMPLE = 2
    TMEMTEST = 3
    TPRINTP = 4
    TPROJECT = 5
    TISOSTAB = 6
    TREGENEXT = 7
    def __init__(self, system, tracktype=TZERODIM, config={}, **kwargs):
        kkeys = kwargs.keys()
        ckeys = [k.lower() for k in config.keys()]
        msg = ''
        if tracktype not in range(-4,8):
            msg = 'specify an integer TrackType between -4 and 7 (inclusive)'
            raise ValueError(msg)
        else:
            self._tracktype = tracktype
        if tracktype in (self.TEVALP, self.TEVALPJ, self.TNEWTP, self.TNEWTPJ, self.TMEMTEST, self.TISOSTAB) and 'start' not in kkeys:
            msg = "specify a point or points to evaluate with the keyword argument `start'"
        elif 'start' in kkeys:
            self._start = kwargs['start']
        if tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTP, self.TPROJECT, self.TREGENEXT) and 'component' not in kkeys:
            msg = "specify a component with the keyword argument `component'"
        elif 'component' in kkeys:
            self._component = kwargs['component']
        if tracktype == self.TSAMPLE and 'sample' not in kkeys:
            msg = "specify how many points to sample with the keyword argument `sample'"
        elif tracktype == self.TSAMPLE:
            self._sample = kwargs['sample']
        if tracktype == self.TPROJECT and 'projection' not in kkeys:
            msg = "specify the variables onto which you wish to project with the keyword argument `projection'"
        elif 'projection' in kkeys:
            self._projection = kwargs['projection']
        if msg:
            raise KeyError(msg)
        
        # tolerance
        if 'tol' in kkeys:
            self._tol = kwargs['tol']
        else:
            self._tol = TOL
        
        # parameter homotopy
        if 'parameterhomotopy' in ckeys:
            ckeys2 = config.keys()
            for k in ckeys2:
                if k.lower() == 'parameterhomotopy':
                    phtpy = k
                    break
        msg = ''
        if 'parameterhomotopy' in ckeys and config[phtpy] > 0 and not system.parameters:
            msg = "you have attempted to define a parameter homotopy on a system with no parameters!"
        elif 'parameterhomotopy' in ckeys and config[phtpy] == 0 and system.parameters:
            msg = "a parameterized system requires ParameterHomotopy > 0"
        elif system.parameters and 'parameterhomotopy' not in ckeys:
            msg = "a parameterized system requires ParameterHomotopy > 0"
            
        if msg:
            raise KeyError(msg)
        self._parameterhomotopy = {}
        if tracktype != self.TZERODIM and 'parameterhomotopy' in ckeys:
            # different capitalizations are valid to Bertini
            ckeys2 = config.keys()
            for k in ckeys2:
                if k.lower() == 'parameterhomotopy':
                    phtpy = k
                    break
            if config[phtpy] > 0:
                msg = 'ParameterHomotopy is only appropriate for zero-dimensional solves'
                raise KeyError(msg)
        
        if 'start' in kkeys:
            start = kwargs['start']
            if type(start) not in (tuple, list):
                start = [start]
            self._start = start
        if 'start_parameters' in kkeys:
            startp = kwargs['start_parameters']
            if type(startp) not in (tuple, list):
                startp = [startp]
            self._parameterhomotopy['start parameters'] = startp
        if 'final_parameters' in kkeys:
            finalp = kwargs['final_parameters']
            if type(finalp) not in (tuple, list):
                finalp = [finalp]
            self._parameterhomotopy['final parameters'] = finalp
        
        # user did not specify start or final parameters
        if 'parameterhomotopy' in ckeys and config[phtpy] > 1:
            msg = "specify start and/or final parameters with the keyword arguments `start_parameters' and/or `final_parameters'"
            raise KeyError(msg)
            
        from tempfile import mkdtemp
        self._dirname = mkdtemp(prefix=basedir)
        self._bertini = BERTINI
        self._system = system
        self._config = config
        self._complete = False
        self._inputf = []
        
    def _proc_err_output(self, output):
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
        
    def _recover_data(self):
        """
        recover the information pertinent to a run
        """
        
        if not self._complete:
            return
        
        from naglib.bertini.data import get_components
        from naglib.bertini.fileutils import parse_witness_data, read_points, striplines
        dirname = self._dirname
        system = self._system
        tol = self._tol
        tracktype = self._tracktype
        main_data = dirname + '/main_data'
        fh = open(main_data, 'r')
        self._main_data = striplines(fh.readlines())
        fh.close()
        
        projective = not not system.homvar
        
        if tracktype == self.TEVALP:
            pass
        elif tracktype == self.TEVALPJ:
            pass
        elif tracktype == self.TNEWTP:
            pass
        elif tracktype == self.TNEWTPJ:
            pass
        elif tracktype == self.TZERODIM:
            finites = dirname + '/finite_solutions'
            reals = dirname + '/read_finite_solutions'
            finite_solutions = read_points(finites, tol=tol, projective=projective)
            real_solutions = read_points(reals, tol=tol, projective=projective)
            
            return finite_solutions, real_solutions
        elif tracktype == self.TPOSDIM:
            wdfile = dirname + '/witness_data'
            self._witness_data = parse_witness_data(wdfile)
            components = get_components(system, self._witness_data)
                        
            return components
        elif tracktype == self.TSAMPLE:
            wdfile = dirname + '/witness_data'
            self._witness_data = parse_witness_data(wdfile)
            samplef = dirname + '/sampled'
            sampled = read_points(samplef, tol=tol, projective=projective)
            
            return sampled
        elif tracktype == self.TMEMTEST:
            wdfile = dirname + '/witness_data'
            self._witness_data = parse_witness_data(wdfile)
            inmat = 'incidence_matrix'
            
        elif tracktype == self.TPRINTP:
            pass
        elif tracktype == self.TPROJECT:
            pass
        elif tracktype == self.TISOSTAB:
            pass
        elif tracktype == self.TREGENEXT:
            pass
        
        
    def _recover_input(self):
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
        while inlines[0] == '\n':
            inlines = inlines[1:]
            
        dex = inlines.index('END;\n')
        inlines = inlines[:dex+1]
        
        return inlines
        
    def _write(self, ifile=False):
        from os.path import exists
        
        tracktype = self._tracktype
        dirname = self._dirname
        system = self._system
        if not exists(dirname):
            from os import mkdir
            mkdir(dirname)
        
        ### write the system
        sysconfig = self._config.copy()
        sysconfig.update({'TrackType':tracktype})
        inputf = self._write_system(system, config=sysconfig)
        
        ### write out `start', `start_parameters', `final_parameters'
        if '_start' in dir(self):
            start = self._start
            startfile = dirname + '/start'
            fprint(start, startfile)
        if self._parameterhomotopy:
            from naglib.bertini.fileutils import fprint
            phtpy = self._parameterhomotopy
            pkeys = phtpy.keys()
            if 'start parameters' in pkeys:
                startp = phtpy['start parameters']
                startpfile = dirname + '/start_parameters'
                fprint(startp, startpfile)
            if 'final parameters' in pkeys:
                finalp = phtpy['final parameters']
                finalpfile = dirname + '/final_parameters'
                fprint(finalp, finalpfile)
        
        ### write out component information
        if '_component' in dir(self):
            component = self._component
            cid = component.component_id
            dim = component.dim
            if tracktype == self.TREGENEXT:
                self._write_system(component.system, 'input_old', {'TrackType':1})
                component.write_witness(dirname, 'witness_data_old')
                instructions = ['1', 'input_old', 'witness_data_old', str(dim), str(cid)]
                self._write_instructions(instructions)
            elif tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTP, self.TPROJECT):
                component.write_witness(dirname, 'witness_data')
                if tracktype == self.TSAMPLE:
                    sample = self._sample
                    instructions = [str(dim), str(cid), str(sample), '0', 'sampled']
                    self._write_instructions(instructions)
                elif tracktype == self.TPRINTP:
                    instructions = [str(dim), str(cid), 'points.out', 'sys.out']
                    self._write_instructions(instructions)
                elif tracktype == self.TPROJECT:
                    instructions = [str(dim), str(cid)]
                    self._write_instructions(instructions)
        
                    ### write out projection information
                    projection = self._projection
                    projnum = ['1' if x in projection else '0' for x in system.variables]
                    projnum = [' '.join(projnum)]
                    self._write_instructions(projnum, 'projection')
                    
        return inputf
    
    def _write_instructions(self, lines, filename='instructions'):
        filename = self._dirname + '/' + filename
        fh = open(filename, 'w')
        for line in lines:
            fh.write(line + '\n')
        fh.close()
        
    def _write_system(self, system, inputf='input', config=None):
        from re import sub as resub
        if not config:
            config = self._config
        dirname = self._dirname
        filename = dirname + '/' + inputf
    
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
        
        return filename
    
    def rerun(self):
        if not self._complete:
            return self.run()
        else:
            inputf = self._inputf
            fh = open(self._dirname + '/input', 'w')
            for line in inputf:
                fh.write(line)
            fh.close()
            return self.run()
                
    def run(self):
        cmd = self._bertini
        if not cmd:
            raise NoBertiniException()
        from os.path import exists
        dirname = self._dirname
        if exists(dirname + '/instructions'):
            stdin = dirname + '/instructions'
        else:
            stdin = None
        input_file = self._write()
        
        arg = [cmd, input_file]
            
        chdir(dirname)
        if stdin:
            stdin = open(stdin, 'r')
        try:
            output = check_output(arg, stdin=stdin)
        except CalledProcessError as e:
            msg = self._proc_err_output(e.output)
            raise BertiniError(msg)

        if stdin:
            stdin.close()
            
        self._complete = True
        self._inputf = self._recover_input()
        self._output = output
        data = self._recover_data()
        
        return data
        
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
        return self._inputf
    @property
    def output(self):
        return self._output

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