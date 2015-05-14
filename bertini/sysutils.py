from __future__ import print_function

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
        """
        """
        kkeys = kwargs.keys()
        ckeys = [k.lower() for k in config.keys()]
        if tracktype not in range(-4,8):
            msg = 'specify an integer TrackType between -4 and 7 (inclusive)'
            raise ValueError(msg)
        else:
            self._tracktype = tracktype
            
        # check to see if tracktype jives with kwargs
        msg = ''
        ## start point(s) required
        if tracktype in (self.TEVALP, self.TEVALPJ, self.TNEWTP, self.TNEWTPJ, self.TMEMTEST, self.TISOSTAB) and 'start' not in kkeys:
            msg = "specify a point or points to evaluate with the keyword argument `start'"
        elif 'start' in kkeys:
            self._start = kwargs['start']
        ## component required
        if tracktype in (self.TSAMPLE, self.TMEMTEST, self.TPRINTP, self.TPROJECT, self.TREGENEXT) and 'component' not in kkeys:
            msg = "specify a component with the keyword argument `component'"
        elif 'component' in kkeys:
            self._component = kwargs['component']
        ## sample count required
        if tracktype == self.TSAMPLE and 'sample' not in kkeys:
            msg = "specify how many points to sample with the keyword argument `sample'"
        elif tracktype == self.TSAMPLE:
            self._sample = kwargs['sample']
        ## projection variables required
        if tracktype == self.TPROJECT and 'projection' not in kkeys:
            msg = "specify the variables onto which you wish to project with the keyword argument `projection'"
        elif 'projection' in kkeys:
            self._projection = kwargs['projection']
        if msg:
            raise KeyError(msg)
        
        # tolerance
        # TODO:make this mean something
        if 'tol' in kkeys:
            self._tol = kwargs['tol']
        else:
            self._tol = TOL
        
        # parameter homotopy
        self._parameter_homotopy = {'key':'', 'arg':0}
        if 'parameterhomotopy' in ckeys:
            ckeys2 = config.keys()
            for k in ckeys2:
                if k.lower() == 'parameterhomotopy':
                    self._parameter_homotopy['key'] = k
                    self._parameter_homotopy['arg'] = config[k] # \in (0,1,2)
                    break
            del ckeys2
        
        # ensure the system jives with the call for parameter homotopy
        msg = ''
        if not system.parameters and self._parameter_homotopy['arg'] > 0:
            msg = "you have attempted to define a parameter homotopy on a system with no parameters!"
        elif system.parameters and self._parameter_homotopy['arg'] <= 0:
            msg = "a parameterized system requires ParameterHomotopy > 0"
        elif tracktype != self.TZERODIM and self._parameter_homotopy['arg'] > 0:
            msg = "parameter homotopy only supported for zero-dimensional runs"
            
        if msg:
            raise KeyError(msg)
        
        if 'start' in kkeys:
            start = kwargs['start']
            if type(start) not in (tuple, list):
                start = [start]
            # this doesn't go in self._parameter_homotopy because other kinds of run use start files
            self._start = start
        if 'start_parameters' in kkeys:
            startp = kwargs['start_parameters']
            if type(startp) not in (tuple, list):
                startp = [startp]
            self._parameter_homotopy['start parameters'] = startp
        else:
            startp = None
        if 'final_parameters' in kkeys:
            finalp = kwargs['final_parameters']
            if type(finalp) not in (tuple, list):
                finalp = [finalp]
            self._parameter_homotopy['final parameters'] = finalp
        else:
            finalp = None
        
        # user did not specify start or final parameters
        if 'parameterhomotopy' in ckeys and self._parameter_homotopy['arg'] > 1:
            if not (startp or finalp):
                msg = "specify start and/or final parameters with the keyword arguments `start_parameters' and/or `final_parameters'"
                raise KeyError(msg)
            
        from tempfile import mkdtemp
        self._dirname = mkdtemp(prefix=basedir)
        self._bertini = BERTINI
        self._system = system
        self._config = config
        self._complete = False
        self._inputf = []
        
    def _parse_witness_data(self, filename):
        """
        Parse witness_data file into usable data
        
        Keyword arguments:
        filename -- string, path to witness_data file
        """
        from os.path import isfile
        from sympy import I, Integer, Float, Rational, Matrix
        from naglib.core.misc import dps, striplines
        from naglib.exceptions import UnclassifiedException
        
        if not isfile(filename):
            msg = "{0} does not exist".format(filename)
            raise IOError(msg)
        
        fh = open(filename, 'r')
        lines = striplines(fh.readlines())
        fh.close()

        num_vars, nonempty_codims = int(lines[0]), int(lines[1])
        # previous line includes additional homogenizing variable(s),
        # as appropriate
        lines = lines[2:]

        # following block represents a single codim; repeated for each
        codims = []
        for i in range(nonempty_codims):
            codim = int(lines[0])
            codims.append({'codim':codim})
            num_points = int(lines[1])
            lines = lines[2:]
            # following block represents a single point; repeated for each
            pts = []
            for j in range(num_points):
                prec = int(lines[0])

                lines = lines[1:]
                pt = []
                for k in range(num_vars):
                    real,imag = lines[k].split(' ')
                    pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))
                pt = Matrix(pt)
                lines = lines[num_vars:]
                # the next point is the last approximation of the point
                # on the path before convergence
                prec = int(lines[0])
                lines = lines[1:]
                approx_pt = []
                for k in range(num_vars):
                    real,imag = lines[k].split(' ')
                    approx_pt.append(Float(real, dps(real)) + I*Float(imag, dps(imag)))

                lines = lines[num_vars:]
                condition_number = float(lines[0])
                corank = int(lines[1]) # corank of Jacobian at this point
                smallest_nonzero_singular = float(lines[2])
                largest_zero_singular = float(lines[3])
                pt_type = int(lines[4])
                multiplicity = int(lines[5])
                component_number = int(lines[6])
                if component_number == -1:
                    msg = "components in {0} have unclassified points".format(filename)
                    raise UnclassifiedException(msg)
                deflations = int(lines[7])
                lines = lines[8:]
                pts.append({'coordinates':pt,
                            'corank':corank,
                            'condition number':condition_number,
                            'smallest nonzero':smallest_nonzero_singular,
                            'largest zero':largest_zero_singular,
                            'type':pt_type,
                            'multiplicity':multiplicity,
                            'component number':component_number,
                            'deflations':deflations,
                            'precision':prec,
                            'last approximation':approx_pt})
            codims[-1]['points'] = pts

        # -1 designates the end of witness points
        lines = lines[1:]

        INT = 0
        DOUBLE = 1
        RATIONAL = 2

        # remaining data is related to slices, randomization,
        # homogenization, and patches
        num_format = int(lines[0])
        # previous line describes format for remainder of data
        lines = lines[1:]

        # the following block is repeated for each nonempty codim.
        # first, matrix A used for randomization
        # second, matrix W
        for i in range(nonempty_codims):
            num_rows, num_cols = lines[0].split(' ')
            num_rows = int(num_rows)
            num_cols = int(num_cols)
            AW_size = num_rows*num_cols

            lines = lines[1:]

            if AW_size == 0:
                A = None
                W = None
            else:
                A = lines[:AW_size]
                lines = lines[AW_size:]
                W = lines[:AW_size]
                lines = lines[AW_size:]

                A = [a.split(' ') for a in A] # A is complex-valued
                if num_format == INT:
                    A = [Integer(a[0]) + I*Integer(a[1]) for a in A]
                elif num_format == DOUBLE:
                    for j in range(len(A)):
                        a = A[j]
                        real,imag = a.split(' ')
                        A[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
                elif num_format == RATIONAL:
                    A = [Rational(a[0]) + I*Rational(a[1]) for a in A]
                A = [A[j:j+num_cols] for j in range(0,AW_size,num_cols)]
                A = Matrix(A)

                W = [int(w) for w in W] # W is integer-valued
                W = [W[j:j+num_cols] for j in range(0,AW_size,num_cols)]
                W = Matrix(W)

            # third, a vector H used for homogenization
            # random if projective input
            H_size = int(lines[0])
            lines = lines[1:]
            H = lines[:H_size]
            H = [h.split(' ') for h in H] # H is complex-valued
            if num_format == INT:
                H = [Integer(h[0]) + I*Integer(h[1]) for h in H]
            elif num_format == DOUBLE:
                for j in range(len(H)):
                    h = H[j]
                    real,imag = h.split(' ')
                    H[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                H = [Rational(h[0]) + I*Rational(h[1]) for h in H]

            H = Matrix(H)
            lines = lines[H_size:]

            # fourth, a number homVarConst
            # 0 for affine, random for projective
            hvc = lines[0].split(' ')
            if num_format == INT:
                hvc = Integer(hvc[0]) + I*Integer(hvc[1])
            elif num_format == DOUBLE:
                real,imag = hvc
                hvc = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                hvc = Rational(hvc[0]) + I*Rational(hvc[1])

            lines = lines[1:]

            # fifth, matrix B for linear slice coefficients
            num_rows, num_cols = lines[0].split(' ')
            num_rows, num_cols = int(num_rows), int(num_cols)
            B_size = num_rows*num_cols
            lines = lines[1:]

            if B_size == 0:
                B = None
            else:
                B = lines[:B_size]
                lines = lines[B_size:]

                B = [b.split(' ') for b in B] # B is complex-valued
                if num_format == INT:
                    B = [Integer(b[0]) + I*Integer(b[1]) for b in B]
                elif num_format == DOUBLE:
                    for j in range(len(B)):
                        real,imag = B[j]
                        B[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
                elif num_format == RATIONAL:
                    B = [Rational(b[0]) + I*Rational(b[1]) for b in B]
                B = [B[j:j+num_cols] for j in range(0,B_size,num_cols)]
                B = Matrix(B)

            # sixth and finally, vector p for patch coefficients
            p_size = int(lines[0])
            lines = lines[1:]

            p = lines[:p_size]
            lines = lines[p_size:]
            
            p = [q.split(' ') for q in p]
            if num_format == INT:
                p = [Integer(q[0]) + I*Integer(q[1]) for q in p]
            elif num_format == DOUBLE:
                for j in range(len(p)):
                    real,imag = p[j]
                    p[j] = Float(real, dps(real)) + I*Float(imag, dps(imag))
            elif num_format == RATIONAL:
                p = [Rational(q[0]) + I*Rational(q[1]) for q in p]
            
            p = Matrix(p)
            codims[i]['A'] = A
            codims[i]['W'] = W
            codims[i]['H'] = H
            codims[i]['homVarConst'] = hvc
            codims[i]['slice'] = B
            codims[i]['p'] = p

        return codims
        
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
        
    def _recover_components(self, witness_data):
        """
        """
        from sympy import sympify
        from naglib.core.algebra import LinearSlice
        from naglib.core.base import AffinePoint, ProjectivePoint
        from naglib.core.geometry import IrreducibleComponent
        from naglib.core.witnessdata import WitnessPoint, WitnessSet
        system = self._system
        variables = system.variables
        if system.homvar:
            proj_dim = len(variables)
            #homsys = system
        else:
            proj_dim = len(system.variables) + 1
            #homsys = system.homogenize()
        homvar = sympify('_homvar')
        while homvar in system.variables:
            homvar = sympify('_' + str(homvar))
        homsys  = system.homogenize(homvar)
        homvars = homsys.variables
        
        components = []
        
        for c in witness_data:
            codim       = c['codim']
            homVarConst = c['homVarConst']
            points      = c['points']
            coeffs      = c['slice']
            
            comp_isprojective = homVarConst == 0
            if coeffs:
                if comp_isprojective:
                    slice = LinearSlice(coeffs, homvars, homvar)
                    if not system.homvar:
                        slice = slice.dehomogenize()
                else:
                    slice = LinearSlice(coeffs, variables)
            else:
                slice = None
                
            dim_list = {}
            
            for point in points:
                comp_id = point['component number']
                if comp_isprojective:
                    coord = ProjectivePoint(point['coordinates'])
                    if not system.homvar:
                        coord = coord.dehomogenize()
                else:
                    coord = AffinePoint(point['coordinates'])
                    
                wpoint = WitnessPoint(coord, comp_id)
                
                if not dim_list.has_key(comp_id):
                    dim_list[comp_id] = []
                    
                dim_list[comp_id].append(wpoint)
            
            for comp_id in dim_list.keys():
                ws = WitnessSet(system.copy(), slice, dim_list[comp_id], witness_data)
                component = IrreducibleComponent(ws, codim, comp_id)
                components.append(component)
                
        return components
        
    def _recover_data(self):
        """
        recover the information pertinent to a run
        """
        
        if not self._complete:
            return
        
        from naglib.bertini.fileutils import read_points
        from naglib.core.misc import striplines
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
            startp = dirname + '/start_parameters'
            finite_solutions = read_points(finites, tol=tol, projective=projective)
            
            pkeys = self._parameter_homotopy.keys()
            if 'start parameters' in pkeys and 'start parameters' not in pkeys:
                start_parameters = read_points()
            
            return finite_solutions
        elif tracktype == self.TPOSDIM:
            wdfile = dirname + '/witness_data'
            self._witness_data = self._parse_witness_data(wdfile)
            components = self._recover_components(self._witness_data)
                        
            return components
        elif tracktype == self.TSAMPLE:
            wdfile = dirname + '/witness_data'
            self._witness_data = self._parse_witness_data(wdfile)
            samplef = dirname + '/sampled'
            sampled = read_points(samplef, tol=tol, projective=projective)
            
            return sampled
        elif tracktype == self.TMEMTEST:
            wdfile = dirname + '/witness_data'
            self._witness_data = self._parse_witness_data(wdfile)
            inmat = dirname + '/incidence_matrix'
            fh = open(inmat, 'r')
            lines = [l.strip() for l in fh.readlines()]
            fh.close()
            
            nonempty_codims = int(lines[0])
            codims = {}
            lines = lines[1:]
            for i in range(nonempty_codims):
                codim, numcomponents = lines[0].split(' ')
                codims[int(codim)] = int(numcomponents)
                lines = lines[1:]
                
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
        
    def _write_files(self):
        from os.path import exists
        from naglib.bertini.fileutils import fprint
        
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
        if self._parameter_homotopy:
            phtpy = self._parameter_homotopy
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
                witness_data = component.witness_data
                self._write_witness_data(witness_data, dirname, components=[component])
                if tracktype == self.TSAMPLE:
                    sample = self._sample
                    instructions = [str(dim), '0', str(sample), '0', 'sampled']
                    self._write_instructions(instructions)
                elif tracktype == self.TPRINTP:
                    instructions = [str(dim), '0', 'points.out', 'sys.out']
                    self._write_instructions(instructions)
                elif tracktype == self.TPROJECT:
                    instructions = [str(dim), '0']
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
    
    def _write_witness_data(self, witness_data, dirname, components=[], filename='witness_data'):
        """
        """
        from sympy import Integer, Float
        fh = open(dirname + '/' + filename, 'w')
        
        if components:
            codims = [] # modified witness_data
            compids = [(c.codim, c.component_id) for c in components]
            unique_codims = sorted(list(set([c[0] for c in compids])))
            for i in range(len(witness_data)):
                wd = witness_data[i]
                codim  = wd['codim']
                if codim in unique_codims:
                    cdict = {'codim':codim,
                             'A':wd['A'],
                             'homVarConst':wd['homVarConst'],
                             'p':wd['p'],
                             'slice':wd['slice'],
                             'W':wd['W'],
                             'H':wd['H'],
                             'points':[]}
                    
                    cids   = [c[1] for c in compids if c[0] == codim]
                    # components may not have correct order
                    cidmap = dict([(cids[i], i) for i in range(len(cids))])
                    points = wd['points']
                    for p in points:
                        cn = p['component number']
                        if cn in cids:
                            cdict['points'].append({
                            'largest zero':p['largest zero'],
                            'precision':p['precision'],
                            'last approximation':p['last approximation'],
                            'smallest nonzero':p['smallest nonzero'],
                            'deflations':p['deflations'],
                            'component number':cidmap[cn],
                            'multiplicity':p['multiplicity'],
                            'corank':p['corank'],
                            'coordinates':p['coordinates'],
                            'condition number':p['condition number'],
                            'type':p['type']})
                    
                    codims.append(cdict)
        else:
            codims = witness_data
            
        nonempty_codims = len(codims)
        num_vars = len(codims[0]['points'][0]['coordinates'])
        fh.write('{0}\n'.format(num_vars))
        fh.write('{0}\n'.format(nonempty_codims))
        
        for i in range(nonempty_codims):
            wd_codim = codims[i]
            fh.write('{0}\n'.format(wd_codim['codim']))
            codim_points = [p for p in wd_codim['points']]
            fh.write('{0}\n'.format(len(codim_points)))
            for p in codim_points:
                prec = p['precision']
                fh.write('{0}\n'.format(prec))
                
                coordinates = p['coordinates']
                for c in coordinates:
                    real,imag = c.as_real_imag()
                    fh.write('{0} {1}\n'.format(real, imag))
                
                fh.write('{0}\n'.format(prec))
                approx = p['last approximation']
                for a in approx:
                    real,imag = a.as_real_imag()
                    fh.write('{0} {1}\n'.format(real, imag))
                fh.write('{0}\n'.format(p['condition number']))
                fh.write('{0}\n'.format(p['corank']))
                fh.write('{0}\n'.format(p['smallest nonzero']))
                fh.write('{0}\n'.format(p['largest zero']))
                fh.write('{0}\n'.format(p['type']))
                fh.write('{0}\n'.format(p['multiplicity']))
                fh.write('{0}\n'.format(p['component number']))
                fh.write('{0}\n'.format(p['deflations']))
        fh.write('-1\n\n') # -1 designates the end of witness points
        
        h1 = codims[0]['H'][0]
        if type(h1) == Integer:
            numtype = 0
        elif type(h1) == Float:
            numtype = 1
        else:
            numtype = 2
        fh.write('{0}\n'.format(numtype))
        
        for i in range(nonempty_codims):
            wd_codim = codims[i]
            A   = wd_codim['A']
            W   = wd_codim['W']
            H   = wd_codim['H']
            hvc = wd_codim['homVarConst']
            B   = wd_codim['slice']
            P   = wd_codim['p']
        
            if A: # also W
                num_rows, num_cols = A.shape
                fh.write('{0} {1}\n'.format(num_rows, num_cols))
                for j in range(num_rows):
                    for k in range(num_cols):
                        real, imag = A[j,k].as_real_imag()
                        fh.write('{0} {1}\n'.format(real, imag))
                for j in range(num_rows):
                    for k in range(num_cols):
                        real, imag = W[j,k].as_real_imag()
                        fh.write('{0} {1}\n'.format(real, imag))   
            else:
                fh.write('1 0\n')
                
            fh.write('\n')
            h = len(H)
            fh.write('{0}\n'.format(h))
            for j in range(h):
                real, imag = H[j].as_real_imag()
                fh.write('{0} {1}\n'.format(real, imag))
            
            fh.write('\n')
            real,imag = hvc.as_real_imag()
            fh.write('{0} {1}\n'.format(real, imag))
            if B:
                num_rows, num_cols = B.shape
                fh.write('{0} {1}\n'.format(num_rows, num_cols))
                for j in range(num_rows):
                    for k in range(num_cols):
                        real, imag = B[j,k].as_real_imag()
                        fh.write('{0} {1}\n'.format(real, imag))
            else:
                fh.write('1 0\n')
            
            p = len(P)
            fh.write('{0}\n'.format(p))
            for j in range(p):
                real, imag = P[j].as_real_imag()
                fh.write('{0} {1}\n'.format(real, imag))
                
        fh.close()
    
    def rerun(self, config={}):
        if not self._complete:
            return self.run()
        else:
            inputf = self._inputf
            for k in config.keys():
                inputf.insert(1, '{0}:{1};\n'.format(k, config[k]))
            fh = open(self._dirname + '/input', 'w')
            for line in inputf:
                fh.write(line)
            fh.close()
            return self.run()
                
    def run(self):
        from os import chdir        
        from os.path import exists
        
        cmd = self._bertini
        if not cmd:
            raise NoBertiniException()
        dirname = self._dirname
        
        input_file = self._write_files()
        
        if exists(dirname + '/instructions'):
            stdin = dirname + '/instructions'
        else:
            stdin = None
        
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