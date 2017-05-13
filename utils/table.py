import os
import PIL.Image
import multiprocessing
import psutil
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
    name = '%d%s%d' % (n, OrbLet[l] if l < len(OrbLet) else '_%d_' % l, m)
    dest = os.path.join(opts['dest'], name + '.png')
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
    if os.path.exists(dest):
        try:
            os.unlink(dest)
        except Exception:
            pass
    img.save(dest, 'PNG', optimize=True)
    # Further optimize the png if optipng is available
    try:
        subprocess.call(['optipng', '-clobber', '-quiet', dest])
    except Exception:
        pass
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
        pool = None
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
        for n in xrange(opts['minn'], opts['maxn'] + 1):
            for l in xrange(0, n):
                for m in xrange(-l if opts['allm'] else 0, l + 1):
                    if pool:
                        tasks.append(pool.apply_async(generate_orbital, args=(opts, n, l, m)))
                    else:
                        generate_orbital(opts, n, l, m)
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
    Supress the ctrl-c signal in the worker processes.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


if __name__ == '__main__':  # noqa
    opts = {
        'orbfile': 'table.orb',
        'size': 128,
        'dest': 'table',
        'minn': 1,
        'maxn': 10,
        'allm': False,
        'multi': False,
        'exe': 'ansiorb.exe',
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
        else:
            help = True
    if help:
        print """Generate a table of orbitals.

Syntax: table.py --anti=(antialiasing factor) --dest=(directory)
    --exe=(ansiorb executable) --minn=(n) --maxn=(n) --orb=(orbital file)
    --size=(pixels) --allm|--posm --multi

--allm calculates positive and negative p values.
--anti is an integer >= 1.  This is an oversampling, which then scales down
    using a PIL function, which may not look any better than the native
    antialiasing.
--dest specifies a destination directory.
--exe specifies the path to the ansiorb.exe executable.
--minn is the minimum n value to calculate.
--maxn is the maximum n value to calculate.
--multi uses multiprocessing to parallelize computations.
--orb is the path of a default orbital to use as a basis for calculation.
--posm calculates only non-negative p values (the is the default).
--size is the output size in pixels.
"""
        sys.exit(0)
    prep_table(opts)
    generate_table(opts)
