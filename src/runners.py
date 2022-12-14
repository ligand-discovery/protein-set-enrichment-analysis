import os
from tqdm import tqdm
import random
import shutil
import pandas as pd

from enrichers import Enricher, GroupSummarizer, Summarizer, Harmonizer, GroupSummarizer
from util import listdir_util

ROOT = os.path.abspath(os.path.dirname(__file__))

MIN_SET_SIZE = 1

# classes

class BaseRunner(object):
    def __init__(self, entity, randomize=True):
        self.randomize = randomize
        self.entity = entity
        self.universes_folder = os.path.join(ROOT, "../data/universes", entity)
        self.signatures_folder = os.path.join(ROOT, "../data/signatures", entity)
        self.annotations_folder = os.path.join(ROOT, "../data/annotations", entity)
        self.results_folder = os.path.join(ROOT, "../results", entity)
        if not os.path.exists(self.results_folder):
            os.makedirs(self.results_folder)

    def tasks_iter(self, file_to_check="meta.json"):
        tasks = []
        universes_folder = self.universes_folder
        annotations_folder = self.annotations_folder
        signatures_folder = self.signatures_folder
        results_folder = self.results_folder
        for universe_file in listdir_util(universes_folder):
            universe_name = universe_file.split(".")[0]
            universe_file = os.path.abspath(
                os.path.join(universes_folder, universe_file)
            )
            p0 = os.path.join(results_folder, universe_name)
            if not os.path.exists(p0):
                os.mkdir(p0)
            for annotations_file in listdir_util(annotations_folder):
                annotations_name = annotations_file.split(".")[0]
                annotations_file = os.path.abspath(
                    os.path.join(annotations_folder, annotations_file)
                )
                min_set_size = MIN_SET_SIZE
                p1 = os.path.join(p0, annotations_name)
                if not os.path.exists(p1):
                    os.mkdir(p1)
                p2 = os.path.join(p1, "global")
                if not os.path.exists(p2):
                    os.mkdir(p2)
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "global")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    if not os.path.exists(p3):
                        os.makedirs(p3, exist_ok=True)
                    for signature_file in listdir_util(
                        os.path.join(signatures_folder, "global", signatures_subfolder)
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "global",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        output_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        if not os.path.exists(output_folder):
                            os.makedirs(output_folder, exist_ok=True)
                        if os.path.exists(os.path.join(output_folder, file_to_check)):
                            continue
                        task = {
                            "signature_file": signature_file,
                            "annotations_file": annotations_file,
                            "universe_file": universe_file,
                            "output_folder": output_folder,
                            "min_set_size": min_set_size,
                        }
                        tasks += [task]
                p2 = os.path.join(p1, "fragment")
                if not os.path.exists(p2):
                    os.mkdir(p2)
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "fragment")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    if not os.path.exists(p3):
                        os.makedirs(p3, exist_ok=True)
                    for signature_file in listdir_util(
                        os.path.join(
                            signatures_folder, "fragment", signatures_subfolder
                        )
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "fragment",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        output_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        if not os.path.exists(output_folder):
                            os.makedirs(output_folder, exist_ok=True)
                        if os.path.exists(os.path.join(output_folder, file_to_check)):
                            continue
                        task = {
                            "signature_file": signature_file,
                            "annotations_file": annotations_file,
                            "universe_file": universe_file,
                            "output_folder": output_folder,
                            "min_set_size": min_set_size,
                        }
                        tasks += [task]
        if self.randomize:
            random.shuffle(tasks)
        for task in tasks:
            yield task


class PrimaryRunner(BaseRunner):
    def __init__(self, entity, randomize=True):
        BaseRunner.__init__(self, entity=entity, randomize=randomize)

    def run(self):
        c = 0
        for task in tqdm(self.tasks_iter()):
            signature_file = task["signature_file"]
            annotations_file = task["annotations_file"]
            universe_file = task["universe_file"]
            if "hek293t_core" not in universe_file:
                continue  # TODO: remove
            c += 1
            # continue
            print(task)
            print(c)
            min_set_size = task["min_set_size"]
            output_folder = task["output_folder"]
            enricher = Enricher(
                signature_file, annotations_file, universe_file, min_set_size
            )
            enricher.calculate()
            if not os.path.exists(output_folder):
                os.makedirs(output_folder, exist_ok=True)
            enricher.save(output_folder)
        print(c)


class SecondaryRunner(object):
    def __init__(self, entity):
        self.entity = entity
        self.universes_folder = os.path.join(ROOT, "../data/universes", entity)
        self.signatures_folder = os.path.join(ROOT, "../data/signatures", entity)
        self.annotations_folder = os.path.join(ROOT, "../data/annotations", entity)
        self.results_folder = os.path.join(ROOT, "../results", entity)

    def is_ranksum(self, results_path):
        task_name = results_path.rstrip("/").split("/")[-1]
        for x in ["_bin_", "_top_", "_bottom_"]:
            if x in task_name:
                return False
        return True

    def tasks_iter(self):
        tasks = []
        universes_folder = self.universes_folder
        annotations_folder = self.annotations_folder
        signatures_folder = self.signatures_folder
        results_folder = self.results_folder
        for universe_file in listdir_util(universes_folder):
            universe_name = universe_file.split(".")[0]
            if "hek293" not in universe_name:
                continue
            universe_file = os.path.abspath(
                os.path.join(universes_folder, universe_file)
            )
            p0 = os.path.join(results_folder, universe_name)
            for annotations_file in listdir_util(annotations_folder):
                annotations_name = annotations_file.split(".")[0]
                annotations_file = os.path.abspath(
                    os.path.join(annotations_folder, annotations_file)
                )
                p1 = os.path.join(p0, annotations_name)
                p2 = os.path.join(p1, "global")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "global")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(signatures_folder, "global", signatures_subfolder)
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "global",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
                p2 = os.path.join(p1, "fragment")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "fragment")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(
                            signatures_folder, "fragment", signatures_subfolder
                        )
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "fragment",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
        for task in tasks:
            yield task

    def run(self):
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            h = Harmonizer(results_folder=results_folder)
            h.harmonize()
            h.save()

    def clean(self):
        for task in tqdm(self.tasks_iter()):
            p = os.path.join(task["results_folder"], "harmon.tsv")
            if os.path.exists(p):
                os.remove(p)


class TertiaryRunner(object):
    def __init__(self, entity):
        self.entity = entity
        self.universes_folder = os.path.join(ROOT, "../data/universes", entity)
        self.signatures_folder = os.path.join(ROOT, "../data/signatures", entity)
        self.annotations_folder = os.path.join(ROOT, "../data/annotations", entity)
        self.results_folder = os.path.join(ROOT, "../results", entity)

    def tasks_iter(self):
        tasks = []
        universes_folder = self.universes_folder
        annotations_folder = self.annotations_folder
        signatures_folder = self.signatures_folder
        results_folder = self.results_folder
        for universe_file in listdir_util(universes_folder):
            universe_name = universe_file.split(".")[0]
            if "hek293" not in universe_name:
                continue
            universe_file = os.path.abspath(
                os.path.join(universes_folder, universe_file)
            )
            p0 = os.path.join(results_folder, universe_name)
            for annotations_file in listdir_util(annotations_folder):
                annotations_name = annotations_file.split(".")[0]
                annotations_file = os.path.abspath(
                    os.path.join(annotations_folder, annotations_file)
                )
                p1 = os.path.join(p0, annotations_name)
                p2 = os.path.join(p1, "global")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "global")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(signatures_folder, "global", signatures_subfolder)
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "global",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
                p2 = os.path.join(p1, "fragment")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "fragment")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(
                            signatures_folder, "fragment", signatures_subfolder
                        )
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "fragment",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
        for task in tasks:
            yield task

    def run(self):
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            print(results_folder)
            h = Summarizer(results_folder=results_folder)
            h.summarize()
            h.save()

    def clean(self):
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            h = Summarizer(results_folder=results_folder)
            p = h.get_summary_folder()
            if os.path.exists(p):
                shutil.rmtree(p)



class OnlyRanksumRunner(object):
    def __init__(self, entity):
        assert entity == "proteins"
        self.data_folder = os.path.abspath(os.path.join(ROOT, "../data/"))
        self.entity = entity
        self.universes_folder = os.path.join(ROOT, "../data/universes", entity)
        self.signatures_folder = os.path.join(ROOT, "../data/signatures", entity)
        self.annotations_folder = os.path.join(ROOT, "../data/annotations", entity)
        self.results_folder = os.path.join(ROOT, "../results", entity)
        self.group_summarizer = GroupSummarizer(self.data_folder)

    def is_log2fc(self, results_path):
        task_name = results_path.rstrip("/").split("/")[-1]
        if "_log2fc" in task_name:
            return True
        else:
            return False

    def tasks_iter(self):
        tasks = []
        universes_folder = self.universes_folder
        annotations_folder = self.annotations_folder
        signatures_folder = self.signatures_folder
        results_folder = self.results_folder
        for universe_file in listdir_util(universes_folder):
            universe_name = universe_file.split(".")[0]
            if "hek293" not in universe_name:
                continue
            universe_file = os.path.abspath(
                os.path.join(universes_folder, universe_file)
            )
            p0 = os.path.join(results_folder, universe_name)
            for annotations_file in listdir_util(annotations_folder):
                annotations_name = annotations_file.split(".")[0]
                annotations_file = os.path.abspath(
                    os.path.join(annotations_folder, annotations_file)
                )
                p1 = os.path.join(p0, annotations_name)
                p2 = os.path.join(p1, "global")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "global")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(signatures_folder, "global", signatures_subfolder)
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "global",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
                p2 = os.path.join(p1, "fragment")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "fragment")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(
                            signatures_folder, "fragment", signatures_subfolder
                        )
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "fragment",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
        for task in tasks:
            if not self.is_log2fc(task["results_folder"]):
                continue
            yield task

    def run(self):
        c = 0
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            print(results_folder)
            self.group_summarizer.set(results_folder=results_folder)
            self.group_summarizer.summarize()
            self.group_summarizer.save()
            c += 1
        print(c)

    def clean(self):
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            self.group_summarizer.set(results_folder=results_folder)
            self.group_summarizer.clean()




class OnlyRanksumConcatenatorRunner(object):
    def __init__(self, entity):
        assert entity == "proteins"
        self.data_folder = os.path.abspath(os.path.join(ROOT, "../data/"))
        self.entity = entity
        self.universes_folder = os.path.join(ROOT, "../data/universes", entity)
        self.signatures_folder = os.path.join(ROOT, "../data/signatures", entity)
        self.annotations_folder = os.path.join(ROOT, "../data/annotations", entity)
        self.results_folder = os.path.join(ROOT, "../results", entity)

    def is_log2fc(self, results_path):
        task_name = results_path.rstrip("/").split("/")[-1]
        if "_log2fc" in task_name:
            return True
        else:
            return False

    def tasks_iter(self):
        tasks = []
        universes_folder = self.universes_folder
        annotations_folder = self.annotations_folder
        signatures_folder = self.signatures_folder
        results_folder = self.results_folder
        for universe_file in listdir_util(universes_folder):
            universe_name = universe_file.split(".")[0]
            if "hek293" not in universe_name:
                continue
            universe_file = os.path.abspath(
                os.path.join(universes_folder, universe_file)
            )
            p0 = os.path.join(results_folder, universe_name)
            for annotations_file in listdir_util(annotations_folder):
                annotations_name = annotations_file.split(".")[0]
                annotations_file = os.path.abspath(
                    os.path.join(annotations_folder, annotations_file)
                )
                p1 = os.path.join(p0, annotations_name)
                p2 = os.path.join(p1, "global")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "global")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(signatures_folder, "global", signatures_subfolder)
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "global",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
                p2 = os.path.join(p1, "fragment")
                for signatures_subfolder in listdir_util(
                    os.path.join(signatures_folder, "fragment")
                ):
                    p3 = os.path.join(p2, signatures_subfolder)
                    for signature_file in listdir_util(
                        os.path.join(
                            signatures_folder, "fragment", signatures_subfolder
                        )
                    ):
                        signature_name = signature_file.split(".")[0]
                        signature_file = os.path.abspath(
                            os.path.join(
                                signatures_folder,
                                "fragment",
                                signatures_subfolder,
                                signature_file,
                            )
                        )
                        results_folder = os.path.abspath(
                            os.path.join(p3, signature_name)
                        )
                        task = {"results_folder": results_folder}
                        tasks += [task]
        for task in tasks:
            if not self.is_log2fc(task["results_folder"]):
                continue
            yield task

    def run(self):
        all_results_path = os.path.join(ROOT, "../results/all_results.tsv")
        c = 0
        started = False
        for task in tqdm(self.tasks_iter()):
            results_folder = task["results_folder"]
            subf = results_folder.rstrip("/").split("/")
            anno = subf[-4]
            prof = subf[-2]
            df = pd.read_csv(os.path.join(results_folder, "harmon.tsv"), delimiter="\t")
            df["annotation"] = anno
            df["profile"] = prof
            if not started:
                df.to_csv(all_results_path, header=True, index=False, sep="\t")
                started = True
            else:
                df.to_csv(all_results_path, mode="a", header=False, index=False, sep="\t")
            c += 1

    def clean(self):
        pass



if __name__ == "__main__":
    PrimaryRunner("proteins").run()
    SecondaryRunner("proteins").run()
    TertiaryRunner("proteins").run()
    OnlyRanksumRunner("proteins").run()
    OnlyRanksumConcatenatorRunner("proteins").run()
