import OpenPNM

'''
quick script to check connectivity options make sense
'''

if __name__ == '__main__':
    scene = OpenPNM.Postprocessing.Graphics.Scene()
    
    for row, l in enumerate([[6, 8, 12], [14, 18, 20], [26]]):
        for col, cval in enumerate(l):
            pn = OpenPNM.Network.Cubic(shape=[3,3,3], connectivity=cval)
            # delete all pairs that do not involve pore #13
            pn['throat.conns'] = [conn for conn in pn['throat.conns'] if 13 in conn]
            offset = [3*col, -3*row, 0]
            wires = OpenPNM.Postprocessing.Graphics.Wires(pn['pore.coords']+offset, pn['throat.conns'])
            scene.add_actors([wires])

    scene.play()
