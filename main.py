import subprocess
import itertools
import decimal
import sys


params = [(n, k, round(0.01 * p, 2)) for n, k, p
          in itertools.product([3], [16], range(0, 55, 5))]

#params = [(2, 64, 0.2), (3, 16, 0.2), (4, 8, 0.2)]
#params = [(n, k, round(0.01 * p, 2)) for n, k, p
#          in itertools.product([3], [16], range(0, 55, 5))]
#params = [(n, k, p) for n, k, p
#          in itertools.product([2], range(8, 65, 4), [0.5])]

for param in params:
    command = "./main " + " ".join(str(i) for i in param)# + " >> compare.csv"

    #print(command)
    proc = subprocess.Popen(
        command,
        shell  = True,
        stdin  = subprocess.PIPE,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    stdout_data, stderr_data = proc.communicate()
    print(stdout_data.decode('utf-8'), end="")
    with open("stderr.log", 'w') as f:
        f.write(stderr_data.decode('UTF-8'))
