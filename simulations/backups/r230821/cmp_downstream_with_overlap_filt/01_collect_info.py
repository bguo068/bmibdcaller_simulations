from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import pickle

res_dir = (
    # "/local/scratch/bing/bmibdcaller_simulations/simulations/r230821/cmp_downstream/res"
    "./res"
)

simulations = [
    dict(genome_set_id="10000", expected_num_peaks=0, label="sp_neu"),
    dict(genome_set_id="10001", expected_num_peaks=14, label="sp_s01"),
    dict(genome_set_id="10002", expected_num_peaks=14, label="sp_s02"),
    dict(genome_set_id="10003", expected_num_peaks=14, label="sp_s03"),
    dict(genome_set_id="10004", expected_num_peaks=14, label="sp_g040"),
    dict(genome_set_id="10005", expected_num_peaks=14, label="sp_g080"),
    dict(genome_set_id="10006", expected_num_peaks=14, label="sp_g120"),
    dict(genome_set_id="10007", expected_num_peaks=14, label="sp_o01"),
    dict(genome_set_id="10008", expected_num_peaks=14, label="sp_o03"),
    dict(genome_set_id="10009", expected_num_peaks=14, label="sp_o27"),
    dict(genome_set_id="20000", expected_num_peaks=0, label="mp_s00"),
    dict(genome_set_id="20001", expected_num_peaks=14, label="mp_s01"),
    dict(genome_set_id="20002", expected_num_peaks=14, label="mp_s02"),
    dict(genome_set_id="20003", expected_num_peaks=14, label="mp_s03"),
    dict(genome_set_id="30000", expected_num_peaks=14, label="sp_rel"),
    dict(genome_set_id="30001", expected_num_peaks=14, label="mp_rel"),
    dict(genome_set_id="30005", expected_num_peaks=0, label="sp_const"),
    dict(genome_set_id="30006", expected_num_peaks=0, label="sp_grow"),
]

ibdcallers = "tskibd hapibd  hmmibd  isorelate  refinedibd  tpbwt".split()

# peaks identified  signal/noise


# ----------------- Signal to noise ratio (per chromosome)   ----------------------
# ----------------- And number of  peaks (expected/observed) ----------------------

# signal/noise ratio results
sn_res = {}
# num of peaks results
np_res = {}
# path to ibd obj
ibdobj_path_res = {}


for simu_ind, simulation in enumerate(simulations):
    # model_ind = 13
    # simulation = simulations[model_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = label.split("_")[0]

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # find the right ibdobj file
        ibdobj_fn = ""
        if model == "sp":
            # e.g. 20003_mp_s03/tskibd/ifm_input/20003_orig.ifm.ibdobj.gz
            ibdobj_fn = next(
                Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ibdne_ibd/").glob(
                    "*_orig.ibdne.ibdobj.gz"
                )
            )
        else:
            # eg: 20003_mp_s03/tskibd/ifm_input/20003_orig.ifm.ibdobj.gz
            ibdobj_fn = next(
                Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_input/").glob(
                    "*_orig.ifm.ibdobj.gz"
                )
            )

        # load the data
        tmp = IBD.pickle_load(ibdobj_fn)

        # calculate signal noise ratio
        sn_ratios = []
        for chrno in range(1, 15):
            cov = tmp._cov_df.loc[
                lambda df: df.Chromosome == chrno, ["Start", "Coverage"]
            ].set_index("Start")
            q5, q95 = np.quantile(cov, [0.05, 0.95])
            noise = cov[(cov >= q5) & (cov <= q95)].mean().to_numpy()[0]
            signal = cov.max().to_numpy()[0]
            sn_ratio = signal / noise
            sn_ratios.append(sn_ratio)
        sn_res[(label, ibdcaller)] = sn_ratios

        num_peaks = tmp._peaks_df.shape[0]
        expected = simulation["expected_num_peaks"]

        np_res[(label, ibdcaller)] = [num_peaks, expected]

        ibdobj_path_res[(label, ibdcaller)] = ibdobj_fn


# ----------------- Ne difference  ------------------------------------------
ne_res = {}

for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = label.split("_")[0]

    # Ne is estimated only for SP model
    if model != "sp":
        continue

    fn = f"{res_dir}/{genome_set_id}_{label}/true_ne/{genome_set_id}_1.true_ne"
    param_popsize = pd.read_csv(fn, sep="\t")
    param_popsize["L95"] = param_popsize["NE"]
    param_popsize["U95"] = param_popsize["NE"]

    ne_res[(label, "param", "param")] = param_popsize

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        fn1 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ne_output").glob(
                "*_orig.ne"
            )
        )
        df = pd.read_csv(fn1, sep="\t")
        df.columns = ["GEN", "NE", "L95", "U95"]
        ne_res[(label, ibdcaller, "orig")] = df

        fn2 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ne_output").glob(
                "*_rmpeaks.ne"
            )
        )
        df = pd.read_csv(fn2, sep="\t")
        df.columns = ["GEN", "NE", "L95", "U95"]
        ne_res[(label, ibdcaller, "rmpeaks")] = df


# Infomap difference
ifm_res = {}
for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = label.split("_")[0]

    # Infomap assignment is estimated only for MP model
    if model != "mp":
        continue

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        # e.g. 20003_mp_s03/tskibd/ifm_output/20003_orig_member.pq
        # print(Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output"))
        fn1 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output").glob(
                "*_orig_member.pq"
            )
        )
        df = pd.read_parquet(fn1)
        ifm_res[(label, ibdcaller, "orig")] = df

        fn2 = next(
            Path(f"{res_dir}/{genome_set_id}_{label}/{ibdcaller}/ifm_output").glob(
                "*_rmpeaks_member.pq"
            )
        )
        df = pd.read_parquet(fn1)
        ifm_res[(label, ibdcaller, "rmpeaks")] = df


# IBD comparison
raw_ibd_file_res = {}

for simu_ind, simulation in enumerate(simulations):
    # simu_ind = 3
    # simulation = simulations[simu_ind]
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    model = label.split("_")[0]

    for caller_ind, ibdcaller in enumerate(ibdcallers):
        # caller_ind = 0
        # ibdcaller = ibdcallers[caller_ind]
        # e.g.  sp_s03/ibd/tskibd/10003_10_tskibd.ibd
        ibd_fn_lst = []
        for chrno in range(1, 15):
            fn = f"{res_dir}/{label}/ibd/{ibdcaller}/{genome_set_id}_{chrno}_{ibdcaller}.ibd"
            # fix inconsistent naming pattern
            if ibdcaller == "tpbwt":
                fn = f"{res_dir}/{label}/ibd/{ibdcaller}ibd/{genome_set_id}_{chrno}_{ibdcaller}ibd.ibd"
            assert Path(fn).exists()
            ibd_fn_lst.append(chrno)
        assert Path(fn).exists()
        raw_ibd_file_res[(label, ibdcaller)] = ibd_fn_lst

# write file
outdir = Path("./analysis_res/info")
outdir.mkdir(parents=True, exist_ok=True)


def pickle_to_file(obj, fn):
    with open(fn, "wb") as f:
        pickle.dump(obj, f)


pickle_to_file(ibdobj_path_res, f"{outdir}/ibdobj_path.pkl")
pickle_to_file(sn_res, f"{outdir}/signal_noise_ratio.pkl")
pickle_to_file(np_res, f"{outdir}/num_peaks.pkl")
pickle_to_file(ne_res, f"{outdir}/ne_res.pkl")
pickle_to_file(ifm_res, f"{outdir}/infomap_res.pkl")
