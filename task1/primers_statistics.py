from Bio import SeqIO
import primer3
import pandas as pd


UIPAC = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "W": ["A", "T"],
    "S": ["C", "G"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "H": ["A", "C", "T"],
    "B": ["C", "G", "T"],
    "V": ["A", "C", "G"],
    "D": ["A", "G", "T"],
    "N": ["G", "A", "T", "C"],
}


def get_non_deg_primers(seq: str) -> list[str]:
    non_deg_primers_from_seq = [seq]
    last_index = 0
    move = True
    pi = 0
    while pi < len(non_deg_primers_from_seq):
        cur_seq = non_deg_primers_from_seq[pi]
        if move:
            move = False
            for ni, nucl in enumerate(cur_seq):
                if nucl in UIPAC:
                    move = True
                    for uipac_nucl in UIPAC[nucl]:
                        new_primer = cur_seq[:ni] + uipac_nucl + cur_seq[ni + 1 :]
                        non_deg_primers_from_seq.append(new_primer)
                    pi += 1
                    break
        else:
            last_index = pi
            break

    return non_deg_primers_from_seq[last_index:]


def calculate_parameters(primer_seq):
    tm = primer3.calc_tm(primer_seq)
    homodimer = primer3.calc_homodimer(primer_seq)
    hairpin = primer3.calc_hairpin(primer_seq)

    return tm, homodimer.dg, homodimer.tm, hairpin.dg, hairpin.tm


def calculate_5th_percentile(x):
    return x.quantile(0.05)


def calculate_95th_percentile(x):
    return x.quantile(0.95)


degenerate_primers = SeqIO.parse("primers_test.fasta", "fasta")

non_degenerate_primers = {}

for primer in degenerate_primers:
    non_degenerate_primers[primer.id] = get_non_deg_primers(str(primer.seq))

primer_data = []

for degenerate_primer, non_degenerate_primers in non_degenerate_primers.items():
    for primer in non_degenerate_primers:
        tm, hom_g, hom_tm, hai_g, hai_tm = calculate_parameters(primer)

        primer_data.append(
            [degenerate_primer, primer, tm, hom_g, hom_tm, hai_g, hai_tm]
        )

df = pd.DataFrame(
    primer_data,
    columns=[
        "degenerate_primer",
        "non_degenerate_primer",
        "tm",
        "homodimer_G",
        "homodimer_tm",
        "hairpin_G",
        "hairpin_tm",
    ],
)

stat_funcs = [
    "mean",
    "median",
    calculate_5th_percentile,
    calculate_95th_percentile,
]

agg_rules = {
    "tm": stat_funcs,
    "homodimer_G": stat_funcs,
    "homodimer_tm": stat_funcs,
    "hairpin_G": stat_funcs,
    "hairpin_tm": stat_funcs,
}

statistics_per_primer = df.groupby("degenerate_primer").agg(agg_rules)

statistics_per_primer.columns = pd.MultiIndex.from_product(
    [
        statistics_per_primer.columns.levels[0],
        ["mean", "median", "5th percentile", "95th percentile"],
    ]
)

overall_stats = df.agg(agg_rules)

overall_stats = overall_stats.rename(
    index={
        "calculate_5th_percentile": "5th percentile",
        "calculate_95th_percentile": "95th percentile",
    }
)

with pd.ExcelWriter("primer_stats.xlsx") as writer:
    statistics_per_primer.to_excel(writer, sheet_name="Primer Stats")
    overall_stats.to_excel(writer, sheet_name="Overall Stats")
