import os
import subprocess
from DatasetHelper import DatasetHelper

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
DATASET_FILE = os.path.join(GIT_REPO_DIR,"Reader/data/JSON/datasets.json")
TEST_FILE = os.path.join(GIT_REPO_DIR,"Reader/data/JSON/test.json")

def main(dry_run=True):
    ds_helper = DatasetHelper()
    ds_helper.load(DATASET_FILE)

    ds_name = "ttllnunuNoHiggs_SM"
    ds_params = {
        "is_eft": False,
        "loc": "/hadoop/store/user/awightma/postLHE_step/2019_04_19/ttllWithNuNu/v1/mAOD_step_ttllnunuNoHiggs_NoDim6_run0/",
        "central_xsec": 0.244109,
        "on_das": False,
        "dataset": ""
    }

    ds_helper.updateDataset(ds_name,**ds_params)
    if not dry_run:
        ds_helper.save(DATASET_FILE)
    else:
        ds_helper.save(TEST_FILE)

if __name__ == "__main__":
    dry_run = False
    main(dry_run)