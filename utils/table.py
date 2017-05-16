import os
import PIL.Image
import multiprocessing
import psutil
import shutil
import signal
import subprocess
import sys
import time


s0Settings = {
    'CameraTheta(rad)': 0,
    'CameraPhi(rad)': 0,
    'CameraPsi(rad)': 0,
    'Type': 'Plane',  # Cutaway
}
OrbLet = 'spdfghiklmnoqrtuvwxyz'


def generate_html(opts):
    """
    Generate an html table for the orbitals that would have been generated.
    Enter: opts: a dictionary of options.
    """
    grid = html_grid(opts)
    html = ['<table>', '<tbody>', '<tr>', '<th></th>']
    for colnum in xrange(len(grid[0])):
        header = []
        spec = next(row[colnum]['col'] for row in grid if row[colnum] is not None)
        for entry in spec:
            header.append('-'.join([
                ('<i>%s</i>' % let) for let in entry.keys()[0].split('-')]) +
                '=%s' % (entry.values()[0]))
        html.append('<th>' + ', '.join(header) + '</th>')
    html.append('</tr>')
    for row in grid:
        html.append('<tr>')
        header = []
        spec = next(row[colnum]['row'] for colnum in xrange(len(row))
                    if row[colnum] is not None)
        for entry in spec:
            header.append('-'.join([
                ('<i>%s</i>' % let) for let in entry.keys()[0].split('-')]) +
                '=%s' % (entry.values()[0]))
        html.append('<th>' + '<br/>'.join(header) + '</th>')
        for cell in row:
            value = ''
            if cell is not None:
                name = get_name(cell['n'], cell['l'], cell['m'])
                value = '<img width="%d" height="%d" src="%s.png"></img>' % (
                    opts['size'], opts['size'], name)
                if opts.get('high'):
                    value = '<a href="%s/%s.png">%s</a>' % (
                        opts['high'], name, value)
            html.append('<td>%s</td>' % value)
        html.append('</tr>')
    html.extend(['</tbody>', '</table>'])
    open(opts['html'], 'w').write('\n'.join(html))


def generate_orbital(opts, n, l, m):
    """
    Generate an orbital, creating a PNG file.
    Enter: opts: a dictionary of options.
           n, l, m: orbital parameters.
    """
    starttime = time.time()
    size = opts['size']
    anti = max(1, opts.get('anti', 1))
    # These two equations were determined empirically to ensure that the s0
    # orbitals of different values of n would have all phases visible and the
    # image would frame the entire orbital.
    psi2 = -11.55 * n / (6.874 + n)
    cz = (39.5 * n ** 2.109 + 20) * 1e-10
    param = {
        'n': '%d' % n,
        'l': '%d' % l,
        'm': '%d' % m,
        'Psi^2(log10)': '%4.2f' % psi2,
        'CameraCenterZ(m)': '%20g' % cz,
        'CameraCx': '%d' % (size * anti / 2),
        'Scale(m)': '%20g' % (cz / 25),
        'FixedSize': 'Yes',
        'FixedWidth': '%d' % (size * anti),
        'FixedHeight': '%d' % (size * anti),
        'BackgroundColor': '0xFFFFFF',
    }
    if n > 1 and not l and not m:
        param.update(s0Settings)
    name = get_name(n, l, m)
    dest = os.path.join(opts['dest'], name + '.png')
    tempdest = os.path.join(opts['dest'], 'temp' + name + '.png')
    description = '%10s %4.2f %g' % (name, psi2, cz * 1e10)
    if not opts['multi']:
        sys.stdout.write(description)
        sys.stdout.flush()
    # Create two orbital specification files by modifying the template.  The
    # two only differ by background color; this allows opacity to be extracted,
    # even though it is not explicitly saved.
    temp1 = os.path.join(opts['dest'], 'temp%sw.orb' % name)
    temp2 = os.path.join(opts['dest'], 'temp%sb.orb' % name)
    orb = update_orb_spec(opts['base'], param)
    open(temp1, 'wb').write(orb)
    proc1 = subprocess.Popen([opts['exe'], temp1], stdout=subprocess.PIPE)
    param['BackgroundColor'] = '0x000000'
    orb = update_orb_spec(opts['base'], param)
    open(temp2, 'wb').write(orb)
    proc2 = subprocess.Popen([opts['exe'], temp2], stdout=subprocess.PIPE)
    # We have to parse P3 images, as PIL doesn't support them.
    data1 = {}
    data2 = {}
    white = [read_to_whitespace(proc1.stdout, data1) for x in xrange(4)]
    black = [read_to_whitespace(proc2.stdout, data2) for x in xrange(4)]
    png = []
    ranti2 = 1.0 / (anti * anti)
    for h in xrange(size):
        white = [float(read_to_whitespace(proc1.stdout, data1))
                 for w in xrange(size * anti * anti * 3)]
        black = [float(read_to_whitespace(proc2.stdout, data2))
                 for w in xrange(size * anti * anti * 3)]
        for w in xrange(size):
            r = g = b = a = 0
            for y in xrange(anti):
                for x in xrange(anti):
                    i = ((y) * size * anti + w * anti + x) * 3
                    al = 255 - (
                        white[i] - black[i] +
                        white[i + 1] - black[i + 1] +
                        white[i + 2] - black[i + 2]) / 3
                    if 0 < al < 255:
                        alpha = al / 255
                        r += black[i] / alpha
                        g += black[i + 1] / alpha
                        b += black[i + 2] / alpha
                        a += al
                    else:
                        r += black[i]
                        g += black[i + 1]
                        b += black[i + 2]
                        a += 255 if al else 0
            png.append((
                max(0, min(255, int(r * ranti2))),
                max(0, min(255, int(g * ranti2))),
                max(0, min(255, int(b * ranti2))),
                max(0, min(255, int(a * ranti2)))))
    os.unlink(temp1)
    os.unlink(temp2)
    proc1 = proc2 = None
    img = PIL.Image.new('RGBA', (size, size))
    img.putdata(png)
    img.save(tempdest, 'PNG', optimize=True)
    # Further optimize the png if optipng is available
    try:
        subprocess.call(['optipng', '-clobber', '-quiet', '-o2', '-zs0-1',
                         '-f0,1,5', tempdest])
    except Exception:
        pass
    if os.path.exists(dest):
        try:
            os.unlink(dest)
        except Exception:
            pass
    shutil.move(tempdest, dest)
    if opts['multi']:
        sys.stdout.write(description)
    filesize = os.path.getsize(dest)
    sys.stdout.write(' %d %3.1fs\n' % (filesize, time.time() - starttime))
    sys.stdout.flush()


def generate_table(opts):
    """
    Generate all orbitals in our range.
    Enter: opts: a dictionary of options.  Modified.
    """
    try:
        pool = tasks = None
        if opts['multi']:
            pool = multiprocessing.Pool(
                initializer=worker_init,
                processes=multiprocessing.cpu_count() / 2
                if opts['multi'] is True else opts['multi'])
            priorityLevel = (psutil.BELOW_NORMAL_PRIORITY_CLASS
                             if sys.platform == 'win32' else 10)
            parent = psutil.Process()
            parent.nice(priorityLevel)
            for child in parent.children():
                child.nice(priorityLevel)
            tasks = []
        generate_table_inner(opts, pool, tasks)
        if pool:
            pool.close()
            for task in tasks:
                task.get()
            pool.join()
    except KeyboardInterrupt:
        if pool:
            try:
                pool.terminate()
                pool.join()
            except Exception:
                pass


def generate_table_inner(opts, pool, tasks):
    """
    Generate all orbitals in our range.  This is the inner loop that can handle
     a multiprocessing pool and task list.
    Enter: opts: a dictionary of options.  Modified.
           pool: None or a multiprocessing pool.
           tasks: None or a list to add multiprocessing tasks to.
    """
    for n in xrange(opts['minn'], opts['maxn'] + 1):
        for l in xrange(0, n):
            for m in xrange(-l if opts['allm'] else 0, l + 1):
                if opts.get('skip'):
                    name = get_name(n, l, m)
                    dest = os.path.join(opts['dest'], name + '.png')
                    if os.path.exists(dest):
                        continue
                if pool:
                    tasks.append(pool.apply_async(generate_orbital, args=(opts, n, l, m)))
                else:
                    generate_orbital(opts, n, l, m)


def get_name(n, l, m):
    """
    Return a file name based on n, l, m.
    Enter: n, l, m: values for the file name.
    Exit:  name: the base name without extension.
    """
    return '%d%s%d' % (n, OrbLet[l] if l < len(OrbLet) else '_%d_' % l, m)


def html_grid(opts):
    """
    Generate a grid with cell information for an html table.
    Enter: opts: a dictionary of options.
    Exit:  grid: an array of arrays, each entry of which is either None or has
                n, l, m values, plus 'col' and 'row' with information on what
                was used to generate the column and row.
    """
    rows, columns = opts['table'].split(':')
    rows = [table_to_dict(opts, part) for part in rows.split(',')]
    if len(rows) == 2:
        rows = [(a, b) for a in rows[0] for b in rows[1]]
    else:
        rows = [(a, ) for a in rows[0]]
    columns = [table_to_dict(opts, part) for part in columns.split(',')]
    if len(columns) == 2:
        columns = [(a, b) for a in columns[0] for b in columns[1]]
    else:
        columns = [(a, ) for a in columns[0]]
    grid = []
    for row in rows:
        gridrow = []
        for col in columns:
            nlm = {'col': col, 'row': row}
            [nlm.update(spec) for spec in col]
            [nlm.update(spec) for spec in row]
            if 'l' not in nlm:
                nlm['l'] = nlm['m'] + nlm['l-m']
            if 'n' not in nlm:
                nlm['n'] = nlm['l'] + nlm['n-l']
            if (nlm['n'] < opts['minn'] or nlm['n'] > opts['maxn'] or
                    nlm['l'] < 0 or nlm['l'] >= nlm['n'] or
                    nlm['m'] < -nlm['l'] or nlm['m'] > nlm['l'] or
                    (nlm['m'] < 0 and not opts['allm'])):
                gridrow.append(None)
                continue
            gridrow.append(nlm)
        if gridrow == [None] * len(columns):
            continue
        grid.append(gridrow)
    for colnum in xrange(len(grid[0]) - 1, -1, -1):
        if [row[colnum] for row in grid] == [None] * len(grid):
            for row in grid:
                row.pop(colnum)
    return grid


def prep_table(opts):
    """
    Load the default orbital.
    Enter: opts: a dictionary of options.  Modified.
    """
    opts['base'] = open(opts['orbfile']).read()


def read_to_whitespace(fptr, track):
    """
    Read, skipping over white space.  Then read all non-whitespace until one
    character of white space is consumed.  Return the non-whitespace that was
    read.  If the end of the file is encountered, return any non-whitespace if
    available, or an empty bytes string.
    Enter: fptr: file-like object to read from.
           track: a dictionary used to hold temporary values
    Exit:  val: non-whitespace string
    """
    if 'data' in track:
        if track['next'] < len(track['data']) - 1:
            out = track['data'][track['next']]
            track['next'] += 1
            return out
        data = track['data'][-1][:-1]
    else:
        data = b''
    track['data'] = (data + fptr.read(1024 * 1024) + b'X').split()
    out = track['data'][0]
    track['next'] = 1
    return out


def table_to_dict(opts, spec):
    """
    Convert a table specification to a list of dictionaries of information.
    Enter: opts: a dictionary of options.
           spec: a table specification: one of n, l, m, n-l, l-m.
    Exit:  speclist: a list, each of which has a dictionary of information with
                     values to iterate through in order.
    """
    if spec == 'n':
        return [{'n': n} for n in range(opts['minn'], opts['maxn'] + 1)]
    if spec == 'l':
        return [{'l': l} for l in range(0, opts['maxn'])]
    if spec == 'm':
        return [{'m': m} for m in range(
            -opts['maxn'] - 1 if opts['allm'] else 0, opts['maxn'])]
    if spec == 'n-l':
        return [{'n-l': nml} for nml in range(1, opts['maxn'] + 1)]
    if spec == 'l-m':
        return [{'l-m': lmm} for lmm in range(0, opts['maxn'])]
    raise 'Invalid table specification %s' % spec


def update_orb_spec(spec, param):
    """
    Update an orbital specification.
    Enter: spec: the text specification.
           param: a dictionary of parameters to update.
    """
    lines = [line for line in spec.replace('\r', '\n').split('\n') if line.strip()]
    out = []
    for line in lines:
        parts = line.split('#')[0].strip().split()
        if len(parts) == 2 and parts[0] in param:
            line = line[:line.index(parts[0]) + len(parts[0]) + 1] + str(param[parts[0]])
        out.append(line)
    out.append('')
    return '\n'.join(out)


def worker_init():
    """
    Suppress the ctrl-c signal in the worker processes.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


if __name__ == '__main__':  # noqa
    opts = {
        'allm': False,
        'dest': 'table',
        'exe': 'ansiorb.exe',
        'high': None,
        'html': False,
        'maxn': 10,
        'minn': 1,
        'multi': False,
        'orbfile': 'table.orb',
        'size': 128,
        'skip': False,
        'table': 'm:n,l',
    }
    help = False
    for arg in sys.argv[1:]:
        if arg == '--allm':
            opts['allm'] = True
        elif arg.startswith('--anti='):
            opts['anti'] = int(arg.split('=', 1)[1])
        elif arg.startswith('--dest='):
            opts['dest'] = arg.split('=', 1)[1]
        elif arg.startswith('--exe='):
            opts['exe'] = arg.split('=', 1)[1]
        elif arg.startswith('--high='):
            opts['high'] = arg.split('=', 1)[1]
        elif arg.startswith('--html='):
            opts['html'] = arg.split('=', 1)[1]
        elif arg.startswith('--maxn='):
            opts['maxn'] = int(arg.split('=', 1)[1])
        elif arg.startswith('--minn='):
            opts['minn'] = int(arg.split('=', 1)[1])
        elif arg == '--multi':
            opts['multi'] = True
        elif arg.startswith('--multi='):
            opts['multi'] = int(arg.split('=', 1)[1])
            if opts['multi'] < 2:
                opts['multi'] = False
        elif arg.startswith('--orb='):
            opts['orbfile'] = arg.split('=', 1)[1]
        elif arg == '--posm':
            opts['allm'] = False
        elif arg.startswith('--size='):
            opts['size'] = int(arg.split('=', 1)[1])
        elif arg == '--skip':
            opts['skip'] = True
        elif arg.startswith('--table='):
            opts['table'] = arg.split('=', 1)[1]
        else:
            help = True
    if help:
        print """Generate a table of orbitals.

Syntax: table.py --anti=(antialiasing factor) --dest=(directory)
    --exe=(ansiorb executable) --minn=(n) --maxn=(n) --orb=(orbital file)
    --size=(pixels) --allm|--posm --skip --multi
    --html=(file) --table=(options) --high=(path)

--allm calculates positive and negative p values.
--anti is an integer >= 1.  This is an oversampling, which then scales down
    using a PIL function, which may not look any better than the native
    antialiasing.
--dest specifies a destination directory.
--exe specifies the path to the ansiorb.exe executable.
--high specifies where high-resolution images are relative to an html table.
    If not specified, the original resolution images will not be link-targets
    for high resolution images.
--html outputs an html file to show the table of orbitals rather than
    generating the orbitals.  This is only a table, without any surrounding
    html tags.
--minn is the minimum n value to calculate.
--maxn is the maximum n value to calculate.
--multi uses multiprocessing to parallelize computations.
--orb is the path of a default orbital to use as a basis for calculation.
--posm calculates only non-negative p values (the is the default).
--skip skips existing orbitals.  The default is to regenerate them.
--size is the output size in pixels.
--table specifies options for an html file with a table.  This is of the form
  (n|l|m|n-l|l-m)[,(n|l|m|n-l|l-m)]:(n|l|m|n-l|l-m)[,(n|l|m|n-l|l-m)] where n
  or n-l, l or l-m, and m must each appear once.  To generate a table of n and
  l versus m, specify "n,l:m".  The original Grand Table is a concatenation of
  6 subtables: "n,l:m", "n,m:l", "l,m:n", "n-l,l-m:m", "n-l,m:l-m", and
  "l-m,m:n-l".  The original main page ad a table with "m:n,l".
"""
        sys.exit(0)
    if opts['html']:
        generate_html(opts)
    else:
        prep_table(opts)
        generate_table(opts)
