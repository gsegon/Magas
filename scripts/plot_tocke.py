import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    tocke1 = pd.read_csv('~/Programs/solver/scripts/tocke1.csv')
    tocke2 = pd.read_csv('~/Programs/solver/scripts/tocke2.csv')

    colors = ['r', 'g', 'b', 'm']

    plt.figure()
    for index, (x, y) in enumerate(zip(tocke1['x'], tocke1['y'])):
        plt.plot(x, y, 'o', color=colors[index])

    for index, (x, y) in enumerate(zip(tocke2['x'], tocke2['y'])):
        plt.plot(x, y, 'o', color=colors[index], alpha=0.5)

    plt.gca().set_aspect('equal')
    plt.show()