import subprocess
import itertools
import decimal
import sys


#params = [(n, k, round(0.01 * p, 2)) for n, k, p
#          in itertools.product(range(2, 6), range(3, 6), range(0, 55, 5))]

params = [(2, 64, 0.2), (3, 16, 0.2), (4, 8, 0.2)]


for param in params[int(sys.argv[1]):]:
    command = "./main " + " ".join(str(i) for i in param) + " >> compare.csv"
    proc = subprocess.Popen(
        command,
        shell  = True,
        stdin  = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout_data, stderr_data = proc.communicate()
    with open("stderr.log", 'w') as f:
        f.write(stderr_data.decode('UTF-8'))
