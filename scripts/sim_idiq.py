import os
import numpy as np

def dq_to_abc(d, q, theta):
    a = d*np.cos(theta) - q*np.sin(theta)
    b = d*np.cos(theta-2*np.pi/3) - q*np.sin(theta-2*np.pi/3)
    c = d*np.cos(theta+2*np.pi/3) - q*np.sin(theta+2*np.pi/3)

    return a, b, c


if __name__ == '__main__':

    os.system('pwd')

    current_angles = np.linspace(0, np.pi/2, 60)
    jds = 10000*np.cos(current_angles)
    jqs = 10000*np.sin(current_angles)

    for index, (jd, jq) in enumerate(zip(jds, jqs)):
        Ja, Jb, Jc = dq_to_abc(jd, jq, 0)

        sources = {'Ja': Ja, '-Ja': -Ja, 'Jb': Jb, '-Jb': -Jb, 'Jc': Jc, '-Jc': -Jc}

        exc_string = '../cmake-build-debug/magas ../cmake-build-debug/motoric_section.json -o motsec-{:d}'.format(index)
        for key, val in sources.items():
            exc_string = exc_string + ' -s{:}={:}'.format(key, val)

        print(exc_string)

        os.system(exc_string)

    print('End of script!')