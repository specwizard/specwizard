def ReadPhys():
    import json
    with open('Phys.data', 'r') as f: 
        x = json.load(f)
    return x
