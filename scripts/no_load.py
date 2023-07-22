import os
import json

if __name__ == '__main__':

    os.system('pwd')
    with open("../examples/motoric_section_airgap/motoric_section_airgap.json", "r") as jsonFile:
        data = json.load(jsonFile)

    for rot_step in range(180):

        data["rotation"]["rot1"] = rot_step

        with open("../examples/motoric_section_airgap/motoric_section_airgap.json", "w") as jsonFile:
            json.dump(data, jsonFile)

        exc_string = '../cmake-build-release/src/magas ../examples/motoric_section_airgap/motoric_section_airgap.json -o out-{:d}.vtu'.format(rot_step)
        os.system(exc_string)

    # for index, (jd, jq) in enumerate(zip(jds, jqs)):
    #     Ja, Jb, Jc = dq_to_abc(jd, jq, 0)
    #
    #     sources = {'Ja': Ja, '-Ja': -Ja, 'Jb': Jb, '-Jb': -Jb, 'Jc': Jc, '-Jc': -Jc}
    #
    #     exc_string = '../cmake-build-debug/magas ../cmake-build-debug/motoric_section.json -o motsec-{:d}'.format(index)
    #     for key, val in sources.items():
    #         exc_string = exc_string + ' -s{:}={:}'.format(key, val)
    #
    #     print(exc_string)
    #
    #     os.system(exc_string)
    #
    # print('End of script!')