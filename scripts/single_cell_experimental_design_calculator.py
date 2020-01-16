from math import *
from itertools import cycle


LIBRARY_PREP_COST = 2000
ILLUMINA_PER_MILLION = 5
MULTIFACTOR = 1.82
R = 0.5714286
M = 4.5997701e-6


def p(multi):
    return multi / (M * (1 - multi))


def q(nr, multi):
    return - (nr * multi) / (R * M * (1 - multi))


def numCellsLoaded(cells, multi):
    return -0.5 * p(multi) - sqrt(0.25 * p(multi) * p(multi) - q(cells, multi))


def multiplet_rate(ncl):
    return M * ncl


def numCellsRecovered(cells, multi):
    return R * numCellsLoaded(cells, multi)


def singlet_rate(ncl):
    return 1 - multiplet_rate(ncl)


def num_singlet(cells, multi):
    return int(singlet_rate(numCellsLoaded(cells, multi)) * numCellsRecovered(cells, multi))


def num_ident_multiplet(cells, multi):
    numMultiplet = numCellsRecovered(cells, multi) - num_singlet(cells, multi)
    return numMultiplet * (multi - 1) / multi


def num_multiplet(cells, multi):
    return numCellsRecovered(cells, multi) - num_singlet(cells, multi)


def num_nonident_multiplet(cells, multi):
    return int(num_multiplet(cells, multi) - num_ident_multiplet(cells, multi))

def readsX(cells, reads_pc, multi):
    nsing = num_singlet(cells, multi)
    nmult = num_multiplet(cells, multi)
    nidentmulti = num_nonident_multiplet(cells, multi)
    return int(cells * reads_pc) / ((nsing / (nsing + MULTIFACTOR * nmult)) + (nidentmulti / (1/MULTIFACTOR * nsing + nmult)))


def singletAvgReadsX(cells, reads_pc, multi):
    rx = readsX(cells, reads_pc, multi)
    nsing = num_singlet(cells, multi)
    nmulti = num_multiplet(cells, multi)
    return int(rx / (nsing + MULTIFACTOR * nmulti))


def multiAvgReadsX(cells, reads_pc, multi):
    rx = readsX(cells, reads_pc, multi)
    nsing = num_singlet(cells, multi)
    nmulti = num_multiplet(cells, multi)
    return int(rx / ((1/MULTIFACTOR) * nsing + nmulti))



def dichotomy(cells, money, multi, eps=0.00001):
    mini = 1
    maxi = 1000000
    for i in range(20):
        midi = int(0.5 * (mini + maxi))
        midi_reads = readsX(cells, midi, multi)
        money2 = midi_reads * ILLUMINA_PER_MILLION / 1000000
        if abs((money2 - money) * 1.0 / money) < eps:
            return midi
        elif money2 > money:
            maxi = midi
        elif money2 <= money:
            mini = midi
    return midi


def exp_design(budget, lo_cell, hi_cell, lo_p, hi_p, diff_cell=250, multi=8):
    # ASSUMPTION 1: number of persons is divisible by multi(=8)
    # ASSUMPTION 2: number of cells is divisible by diff_cell(=250)
    design = {}
    pers = hi_p
    while pers >= lo_p:
        seq_budget = budget - (pers / multi) * LIBRARY_PREP_COST
        if seq_budget < 0:
            design[pers] = []
            break
        # find budget per sequencing batch
        seq_batch_budget = seq_budget / (pers / multi)
        reads_pp = []
        cn = hi_cell
        while cn >= lo_cell:
            cells_batch = cn * multi
            rpp = dichotomy(cells_batch, seq_batch_budget, multi)
            if rpp > 0:
                singlets_ = num_singlet(cells_batch, multi)
                nonident_multi_ = num_nonident_multiplet(cells_batch, multi)
                singlets_reads = singletAvgReadsX(cells_batch, rpp, multi)
                multiplets_reads = multiAvgReadsX(cells_batch, rpp, multi)
                singlets_pic = int(singlets_ / multi)
                nonident_multi_pic = int(nonident_multi_ / multi)
                reads_pp.append((cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads))
            cn -= diff_cell
        if reads_pp:
            design[pers] = reads_pp
        pers -= multi
    return design



def exp_design_fixed_lane_capacity(budget, lo_cell, hi_cell, lo_p, hi_p, diff_cell=250, diff_person=8, capacity=24000, max_multi=16):
    # ASSUMPTION 1: number of cells per lane is maximized
    # ASSUMPTION 2: maximum number of individuals multiplexed is 16
    # Put greedily cells into lanes taking care to not exceed the maximum lane capacity
    # and not to exceed number of multiplexed persons
    design = {}
    pers = hi_p
    while pers >= lo_p:
        cn = hi_cell
        reads_pp = []
        while cn >= lo_cell:
            number_ind_per_lane = int(capacity * 1.0 / cn)
            number_ind_per_lane = min(number_ind_per_lane, max_multi)
            nr_batches = int(pers * 1.0 / number_ind_per_lane)
            if pers % number_ind_per_lane > 0:
                nr_batches += 1
            ##print nr_batches, number_ind_per_lane, capacity, cn, pers
            number_ind_per_lane_approx = pers / nr_batches
            total_seq_budget = budget - LIBRARY_PREP_COST * nr_batches
            if total_seq_budget <= 0:
                break
            ##print total_seq_budget, "BUDDDD"
            seq_budget_per_person = total_seq_budget / pers
            batch_ind_info = [number_ind_per_lane_approx for i in range(nr_batches)]
            extras = pers - sum(batch_ind_info)
            cyc = cycle(range(nr_batches))
            while extras > 0:
                inc_batch = next(cyc)
                batch_ind_info[inc_batch] += 1
                extras -= 1
            ##print batch_ind_info
            batch_money_info = [seq_budget_per_person * batch_ind_info[i] for i in range(len(batch_ind_info))]
            info = []
            for v, w in zip(batch_money_info, batch_ind_info):
                u = cn * w
                rpp = dichotomy(u, v, w)
                ##print rpp, "RPP"
                if rpp > 0:
                    singlets_ = num_singlet(u, w)
                    nonident_multi_ = num_nonident_multiplet(u, w)
                    singlets_reads = singletAvgReadsX(u, rpp, w)
                    multiplets_reads = multiAvgReadsX(u, rpp, w)
                    singlets_pic = int(singlets_ / w)
                    nonident_multi_pic = int(nonident_multi_ / w)
                    info.append((w, cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads))
            if len(info) == len(batch_ind_info):
                # group by batch cell count
                info_singlet_reads_dict = {}
                info_multiplet_reads_dict = {}
                info_singlets_pic_dict = {}
                info_nonident_multi_pic_dict = {}
                for entry in info:
                    info_singlets_pic_dict[entry[2]] = info_singlets_pic_dict.get(entry[2], 0) + entry[0]
                    info_nonident_multi_pic_dict[entry[3]] = info_nonident_multi_pic_dict.get(entry[3], 0) + entry[0]
                    info_singlet_reads_dict[entry[4]] = info_singlet_reads_dict.get(entry[4], 0) + entry[0]
                    info_multiplet_reads_dict[entry[5]] = info_multiplet_reads_dict.get(entry[5], 0) + entry[0]
                singlets_pic = 0
                for u, v in info_singlets_pic_dict.items():
                    singlets_pic += u * v
                singlets_pic /= float(sum(info_singlets_pic_dict.values()))
                singlets_pic = int(singlets_pic)
                nonident_multi_pic = 0
                for u, v in info_nonident_multi_pic_dict.items():
                    nonident_multi_pic += u * v
                nonident_multi_pic /= float(sum(info_nonident_multi_pic_dict.values()))
                nonident_multi_pic = int(nonident_multi_pic)
                singlets_reads = 0
                for u, v in info_singlet_reads_dict.items():
                    singlets_reads += u * v
                singlets_reads /= float(sum(info_singlet_reads_dict.values()))
                singlets_reads = int(singlets_reads)
                multiplets_reads = 0
                for u, v in info_multiplet_reads_dict.items():
                    multiplets_reads += u * v
                multiplets_reads /= float(sum(info_multiplet_reads_dict.values()))
                multiplets_reads = int(multiplets_reads)
                reads_pp.append((cn, singlets_pic, nonident_multi_pic, singlets_reads, multiplets_reads))
            cn -= diff_cell
            if reads_pp:
                design[pers] = reads_pp
        pers -= diff_person
    return design


