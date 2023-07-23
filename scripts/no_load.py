import os
import json

if __name__ == '__main__':

    example = "motoric_section"
    example_path = "../examples/{:}/{:}.json".format(example, example)

    os.system('pwd')
    with open(example_path, "r") as jsonFile:
        data = json.load(jsonFile)

    for rot_step in range(0, 60):

        data["rotation"]["rot1"] = rot_step

        with open(example_path, "w") as jsonFile:
            json.dump(data, jsonFile)

        exc_string = '../cmake-build-release/src/magas {:} -o out-{:d}.vtu'.format(example_path, rot_step)
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