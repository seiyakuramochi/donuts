
lines = '''0.1, 3, 3, 0, 
0.2, 3, 3, 0.001, 
0.3, 3, 3, 0.008, 
0.4, 3, 3, 0.012, 
0.5, 3, 3, 0.06, 
0.1, 3, 4, 0, 
0.2, 3, 4, 0.001, 
0.3, 3, 4, 0.007, 
0.4, 3, 4, 0.026, 
0.5, 3, 4, 0.1, 
0.1, 3, 5, 0, 
0.2, 3, 5, 0.01, 
0.3, 3, 5, 0.024, 
0.4, 3, 5, 0.09, 
0.5, 3, 5, 0.19, 
0.1, 3, 6, 0, 
0.2, 3, 6, 0.008, 
0.3, 3, 6, 0.041, 
0.4, 3, 6, 0.133, 
0.5, 3, 6, 0.287, 
0.1, 3, 7, 0.002, 
0.2, 3, 7, 0.006, 
0.3, 3, 7, 0.06, 
0.4, 3, 7, 0.185, 
0.5, 3, 7, 0.39, 
0.1, 3, 8, 0.001, 
0.2, 3, 8, 0.023, 
0.3, 3, 8, 0.077, 
0.4, 3, 8, 0.205, 
0.5, 3, 8, 0.387, 
0.1, 3, 9, 0.003, 
0.2, 3, 9, 0.02, 
0.3, 3, 9, 0.091, 
0.4, 3, 9, 0.221, 
0.5, 3, 9, 0.482, '''


uniq_p = []
uniq_k = []
for l in lines.split("\n"):
    print(l.strip().split(","))
    p, _, k, unr = tuple(float(v[:4]) for v in l.strip().split(",") if v)
    p = p
    k = int(k)
    if p not in uniq_p:
        uniq_p.append(p)
    if k not in uniq_k:
        uniq_k.append(k)

print("k\p", end="")
for p in uniq_p:
    print(" p="+str(p), end="")
print()

for k in uniq_k:
    print(k, end="")
    for l in lines.split("\n"):
        l = l.strip()
        kk = l.split(",")[-3]
        kv = l.split(",")[-2]
        if int(kk) == k:
            print("", float(kv[:6])*100, end="")
    print()
