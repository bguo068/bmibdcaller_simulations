import json
from pathlib import Path

SP_DEFAULTS = dict(
    seqlen=100 * 15000,
    selpos=round(0.33 * 100 * 15000),
    num_origins=1,
    N=10000,
    h=0.5,
    s=0.3,
    g_sel_start=80,
    r=0.01 / 15_000,
    sim_relatedness=0,
    g_ne_change_start=200,
    N0=1000,
    u=1e-8,
    nsam=1000,
)

MP_DEFAULTS = dict(
    seqlen=100 * 15000,
    selpos=round(100 * 15000 * 0.33),
    num_origins=1,
    N=10000,
    h=0.5,
    s=0.3,
    g_sel_start=80,
    r=0.01 / 15_000,
    sim_relatedness=0,
    mig=1e-5,
    sel_mig=0.01,
    npop=5,
    nsam=200,
    Tsplit=500,
    u=1e-8,
)


# for downstream analyses
def downstream():
    parms_dict = dict(  # single  population
        sp_neu=SP_DEFAULTS | dict(s=0.0, genome_set_id=10000),
        sp_s01=SP_DEFAULTS | dict(s=0.1, genome_set_id=10001),
        sp_s02=SP_DEFAULTS | dict(s=0.2, genome_set_id=10002),
        sp_s03=SP_DEFAULTS | dict(s=0.3, genome_set_id=10003),
        # added repeats for genome_set_id = 10000
        sp_neub=SP_DEFAULTS | dict(s=0.0, genome_set_id=10010),
        sp_neuc=SP_DEFAULTS | dict(s=0.0, genome_set_id=10020),
    )
    Path("sp_sets.json").write_text(json.dumps(parms_dict, indent=2))

    parms_dict = dict(  # multiple  population
        mp_neu=MP_DEFAULTS | dict(s=0.0, genome_set_id=30000),
        # added repeats for genome_set_id = 30000
        mp_neub=MP_DEFAULTS | dict(s=0.0, genome_set_id=30010),
        mp_neuc=MP_DEFAULTS | dict(s=0.0, genome_set_id=30020),
    )
    Path("mp_sets.json").write_text(json.dumps(parms_dict, indent=2))


if __name__ == "__main__":
    downstream()
