import os
import json
import matplotlib.pyplot as plt

if __name__ == '__main__':

    rot_steps = range(60)
    torques = []
    torques_arkkio = []
    flux_a_plus = []
    flux_a_minus = []
    flux_b_plus = []
    flux_b_minus = []
    flux_c_plus = []
    flux_c_minus = []

    for rot_step in rot_steps:
        with open("results-out-{:d}.json".format(rot_step), "r") as jsonFile:
            data = json.load(jsonFile)
            # torques.append(data["Torque"])
            # torques_arkkio.append(data["Torque_section"])
            flux_a_plus.append(data["Flux_a_+"])
            flux_a_minus.append(data["Flux_a_-"])
            flux_b_plus.append(data["Flux_b_+"])
            flux_b_minus.append(data["Flux_b_-"])
            flux_c_plus.append(data["Flux_c_+"])
            flux_c_minus.append(data["Flux_c_-"])

    # plt.figure()
    # plt.plot(rot_steps, torques)
    # plt.plot(rot_steps, torques, '*')
    # plt.plot(rot_steps, torques_arkkio)
    # plt.plot(rot_steps, torques_arkkio, '*')

    flux_a = [a-b for a,b in zip(flux_a_minus, flux_a_plus)]
    flux_b = [a-b for a,b in zip(flux_b_minus, flux_b_plus)]
    flux_c = [a-b for a,b in zip(flux_c_minus, flux_c_plus)]
    plt.figure()

    plt.plot(rot_steps, flux_a)
    plt.plot(rot_steps, flux_a, '*')

    plt.plot(rot_steps, flux_b)
    plt.plot(rot_steps, flux_b, '*')

    plt.plot(rot_steps, flux_c)
    plt.plot(rot_steps, flux_c, '*')

    plt.show()