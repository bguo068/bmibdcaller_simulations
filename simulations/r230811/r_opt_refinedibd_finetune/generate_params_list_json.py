import json
import math

res = []
for lod in [0.3, 0.7, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 1.9]:
    item = dict(
        lod=lod,
        length=2.0,
        scale=math.sqrt(1000/100),
        minmac=20,
        window=40.0
    )
    res.append(item)

with open("params_list.json", 'w') as f:
    json.dump(res, f, indent=4)