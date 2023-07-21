import os
import json
import matplotlib.pyplot as plt

if __name__ == '__main__':

    rot_steps = range(180)
    torques = []
    torques_arkkio = []

    for rot_step in rot_steps:
       with open("results- out-{:d}.json".format(rot_step), "r") as jsonFile:
           data = json.load(jsonFile)
           torques.append(data["Torque"])
           torques_arkkio.append(data["Torque_section"])

    plt.plot(rot_steps, torques)
    plt.plot(rot_steps, torques, '*')
    plt.plot(rot_steps, torques_arkkio)
    plt.plot(rot_steps, torques_arkkio, '*')
    plt.show()