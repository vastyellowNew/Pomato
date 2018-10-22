'''************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body
        problem for orbital mechanics applications.

Authors: Wayne Schlei and Xavier Tricoche

Copyright (c) 2013-2018, Purdue University

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
************************************************************************'''


import sys, argparse
import subprocess as sp
import pandas as pd
import os
import numpy as np
import datetime
import socket
import logging as log

now = datetime.datetime.now()

cwd = os.getcwd()

systems = pd.read_csv('systems.csv')

fixpoint_extractor = 'pmateTopologyExtraction'
manifold_extractor = 'pmateManifoldCreator'
fixpoint_merger = 'pmateMergeFPData'
tabulator = 'pmateTabulateFPX'

def tokenize(args):
    single_str = ' '.join(arg for arg in args)
    tokens = single_str.split();
    # print('tokens: ', tokens)
    return tokens

parser = argparse.ArgumentParser(
    description='Convenient command line interface to control POMATO execution')
parser.add_argument('system', metavar='sys', type=str, help='CR3BP system to analyze')
parser.add_argument('Hamiltonian', metavar='C', type=float, help='Jacobi constant')
parser.add_argument('--path', metavar='binary', type=str, default=cwd + '/../build/bin', help='Path to executables')
parser.add_argument('-o', '--output', metavar='out', type=str, default=cwd, help='Output directory')
parser.add_argument('--params', metavar='par', type=str, default=cwd + '/params/', help='Parameter directory')
parser.add_argument('--tabulate', action='store_true', help='Export result statistics in table form')
parser.add_argument('-e', '--epsilon', metavar='eps', type=float, default=1.0e-8, help='Integration precision')
parser.add_argument('-re', '--refinement_epsilon', metavar='reps', type=float, default=1.0e-12, help='Integration precision during refinement phase')
parser.add_argument('-d', '--depth', metavar='depth', type=int, default=3, help='Subdivision depth for map dicretization')
parser.add_argument('--pmax', metavar='pmax', type=int, default=10, help='Max considered period')
parser.add_argument('--pmin', metavar='pmin', type=int, default=1, help='Min considered period')
parser.add_argument('-n', '--iterations', metavar='niter', type=int, default=50, help='Max number of map iterations for winding number computation')
parser.add_argument('-r', '--resolution', metavar='res', type=int, nargs=2, default=[24, 16], help='Initial sampling resolution')
parser.add_argument('--max_angle', metavar='maxa', type=float, default=2.356194490192345, help='Max allowed angle in manifold construction')
parser.add_argument('-b', '--bounds', metavar='bounds', type=float, nargs=4, help='Computation bounds in (x, xdot)')
parser.add_argument('--name', metavar='name', type=str, help='Name given to selected region', default='NoName')
parser.add_argument('--log', metavar='log', type=str, help='Log file name')
parser.add_argument('--lmin', metavar='lmin', type=float, default=1.0e-3, help='Min allowed distance during index computation')
parser.add_argument('--map_param', type=str, help='Map parameter file')
parser.add_argument('--man_param', type=str, help='Manifold construction parameters')
parser.add_argument('--split', metavar='split', type=int, nargs=2, default=[1, 1], help='Domain decomposition for computation')
parser.add_argument('--skip_fps', action='store_true', help='Skip fixed point extraction and perform merge only if necessary')
parser.add_argument('-v', '--verbose', action='store_true', help='toggle on verbose output')

args = parser.parse_args()

if args.log is None:
    args.log = 'log_C_is_{:.8f}.txt'.format(args.Hamiltonian)

# configure logging output
log.basicConfig(filename=args.log, level=log.DEBUG, format='%(levelname)s:%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

log.info('--------------------------------------------------------')
log.info('pomato.py:')
log.info('hostname: ' + socket.gethostname())
log.info('time    : ' + now.strftime('%m/%d/%Y'))
log.info('--------------------------------------------------------')

# verify that we know where the necessary executables are located
if not os.path.isfile(args.path + '/' + fixpoint_extractor) or \
    not os.path.isfile(args.path + '/' + manifold_extractor) or \
    not os.path.isfile(args.path + '/' + fixpoint_merger) or \
    not os.path.isfile(args.path + '/' + tabulator):
    log.debug('ERROR: Provided path (' + args.path + ') does not contain the required executable(s).\nUnable to proceed.')
    sys.exit()

fixpoint_extractor = args.path + '/' + fixpoint_extractor
manifold_extractor = args.path + '/' + manifold_extractor
fixpoint_merger = args.path + '/' + fixpoint_merger
tabulator = args.path + '/' + tabulator

# access system information
datarow = systems.loc[systems['System Name']==args.system]

mu = datarow['mu'].astype(float).values[0]
r1 = datarow['R1'].astype(float).values[0]
r2 = datarow['R2'].astype(float).values[0]
lstar = datarow['lstar'].astype(float).values[0]

if args.bounds is None:
    args.bounds = [
        float(datarow['xmin']),
        float(datarow['xdotmin']),
        float(datarow['xmax']),
        float(datarow['xdotmax'])
    ]

params_txt = [
    'Selected parameters for computation:',
    'System={}'.format(args.system),
    'mu={}'.format(mu),
    'R1={}'.format(r1),
    'R2={}'.format(r2),
    'lstar={}'.format(lstar),
    'C={}'.format(args.Hamiltonian),
    'params path=' + args.params,
    'tabulate={}'.format(args.tabulate),
    'epsilon={}'.format(args.epsilon),
    'refinement epsilon={}'.format(args.refinement_epsilon),
    'refinement depth={}'.format(args.depth),
    'max period={}'.format(args.pmax),
    'min period={}'.format(args.pmin),
    'resolution={}'.format(args.resolution),
    'max angle={}'.format(args.max_angle),
    'bounds={}'.format(args.bounds),
    'min sampling distance (non-dim)={}'.format(args.lmin),
    'domain split={}'.format(args.split),
    'verbose={}'.format(args.verbose),
    'log file={}'.format(args.log)
]

log.info(' '.join(item for item in params_txt))

if args.split[0] < 1:
    args.split[0] = 1
if args.split[1] < 1:
    args.split[1] = 1

basename = datarow['Symbol'].astype(str).values[0] + '_' + 'PMATE' + '_' + args.name + '_' + '{:.8f}'.format(args.Hamiltonian)
prefix = args.output + '/' + basename
xs = np.linspace(args.bounds[0], args.bounds[2], args.split[0]+1, endpoint=True)
ys = np.linspace(args.bounds[1], args.bounds[3], args.split[1]+1)

log.info('bounds={}'.format(args.bounds))

Nstr = '{}'.format(args.split[0]*args.split[1])
ndigits = len('{}'.format(Nstr))
digit_format = '{:0>' + '{}'.format(ndigits) + 'd}'
fpx_files = []
part_res = [
    args.resolution[0]/args.split[0],
    args.resolution[1]/args.split[1]
]

log.info('Starting fixed point extraction over ' + Nstr + ' subdomains')
fpx_start = datetime.datetime.now()
for i in range(0, args.split[0]):
    for j in range(0, args.split[1]):
        k = i*args.split[1]+j
        xmin = xs[i]
        xmax = xs[i+1]
        ymin = ys[j]
        ymax = ys[j+1]
        idx = digit_format.format(int(k+1)) + 'of' +  digit_format.format(int(Nstr))
        filename = prefix + '_' + 'Part_' + idx
        fpx_files.append(filename + '.fpx')

        do_skip = False
        if args.skip_fps:
            if os.path.isfile(filename + '.fpx'):
                do_skip = True
            else:
                log.info(filename + '.fpx was not found')

        if args.map_param is None:
            args.map_param = filename + '.param'

        cmd_line_args = [
                 fixpoint_extractor,
                 '-e {}'.format(args.epsilon),
                 '-re {}'.format(args.refinement_epsilon),
                 '-b {} {} {} {}'.format(xmin, ymin, xmax, ymax),
                 '-r {} {}'.format(part_res[0], part_res[1]),
                 '-m {}'.format(mu),
                 '-C {}'.format(args.Hamiltonian),
                 '-md {}'.format(args.depth),
                 '-xp {}'.format(args.pmax),
                 '-mp {}'.format(args.pmin),
                 '-p {}'.format(args.iterations),
                 '-R1 {}'.format(r1),
                 '-R2 {}'.format(r2),
                 '-l {}'.format(lstar),
                 '-o ' + filename,
                 '-wp 1',
                 '-wg 1',
                 '-of 1',
                 '-og 1',
                 '-lm {}'.format(args.lmin),
                 '-ma {}'.format(args.max_angle) ]
        if args.verbose:
            cmd_line_args.append('-v 1')

        cmd_line_args = tokenize(cmd_line_args)

        start_time = datetime.datetime.now()
        log.info('Starting fixed point extraction in {}th subdomain'.format(k))

        if do_skip:
            log.info('Skipping fixed point extraction: ' + filename + '.fpx already exists')
        else:
            log.info(' '.join(item for item in cmd_line_args))
            completed = sp.run(cmd_line_args)
        end_time = datetime.datetime.now()
        log.info('Computation terminated')
        log.info('Duration: ' + str(end_time-start_time))

fpx_end = datetime.datetime.now()
log.info('Complete fixed point extraction took ' + str(fpx_end-fpx_start))

log.info('Merging fixed points')
cmd_line_args = [
    fixpoint_merger,
    '-i ' + ' '.join(str(name) for name in fpx_files),
    '-o ' + prefix + '.fpx'
]
log.info(' '.join(item for item in cmd_line_args))
cmd_line_args = tokenize(cmd_line_args)

do_skip = False
if args.skip_fps:
    if os.path.isfile(prefix + '.fpx'):
        log.info('Skipping fixed point merging: ' + prefix + '.fpx already exists')
        do_skip = True
    else:
        log.info(prefix + '.fpx was not found')
if not do_skip:
    completed = sp.run(cmd_line_args)

if args.tabulate:
    log.info('Tabulate results:')
    cmd_line_args = [
        tabulator,
        '-i ' + prefix + '.fpx',
        '-o ' + prefix + '.tex'
    ]
    log.info(' '.join(item for item in cmd_line_args))
    do_skip = False
    cmd_line_args = tokenize(cmd_line_args)
    if args.skip_fps:
        if os.path.isfile(prefix + '.tex'):
            log.info('Skipping tabulation: ' + prefix + '.tex already exists')
            do_skip = True
        else:
            log.info(prefix + '.tex was not found')
    if not do_skip:
        completed = sp.run(cmd_line_args)
    log.info('Tabulation complete')

log.info('Computing manifolds')
if args.map_param is None:
    log.warning('No map parameter file provided. Proceeding with default values')
if args.man_param is None:
    log.warning('No manifold construction parameter file provided. Proceeding with default values')

cmd_line_args = [
    'time',
    '-va',
    '-o timestats.txt',
    manifold_extractor,
    '-C {}'.format(args.Hamiltonian),
    '-m {}'.format(mu),
    '-b {} {} {} {}'.format(args.bounds[0], args.bounds[1], args.bounds[2], args.bounds[3]),
    '-ix ' + prefix + '.fpx',
    '-ip ' + args.map_param,
    '-o ' + prefix + '_manifolds',
    '-R1 {}'.format(r1),
    '-R2 {}'.format(r2),
    '-l {}'.format(lstar)
]

if args.verbose:
    cmd_line_args.append('-v 1')

if args.man_param is not None:
    cmd_line_args.append('-is ' + args.man_param)

log.info(' '.join(item for item in cmd_line_args))

cmd_line_args = tokenize(cmd_line_args)

start = datetime.datetime.now()
log.info('Starting manifold construction')

completed = sp.run(cmd_line_args)

end = datetime.datetime.now()
log.info('Manifold construction complete')
log.info('Duration: ' + str(end-start))

print('Computation log was exported to ' + args.log)
