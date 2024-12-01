import numpy as np
import matplotlib.pyplot as plt
import fatpack
import pandas as pd
from fontTools.ttLib.tables.TupleVariation import PRIVATE_POINT_NUMBERS

import utilities as utils

#********************************CODE BEGINS HERE********************************************

def main():
    M = 0.4
    k = 6.5
    N_1 = 500000
    N_2 = 100000
    Sm_kati = 800
    sa_1_kati = 480
    sa_2_kati = sa_1_kati * (N_2 / N_1) ** (-1 / k)
    sa_1_kati_allo = sa_1_kati + M * Sm_kati
    sa_2_kati_allo = sa_2_kati + M * Sm_kati
    b_ = np.log10(sa_1_kati_allo / sa_2_kati_allo) / np.log10(N_1 / N_2)
    sf_ = sa_1_kati_allo / ((2 * N_1) ** b_)

    # print(b_)
    # print(sf_)
    filename = ('Flucht_filtered2', 'Normal_filtered2', 'Ruttelprogram_filtered2')
    N_T = {}
    N_eq_big_dic = {}
    S_big_dic = {}
    Mean_big_dic = {}
    for z in filename:

    #__________________________INPUT DATA SECTION_________________________________________________________
    #*Everything you will need to change are located in this section
    #You will change the filename of you want to read

    #____________________________________________________________________________________________________

    #Read the file of the filtered data
        df_new = utils.reader(z)

    #Exclude the time column
        measurements_points = df_new.keys()[2:]
        if z == 'Flucht_filtered2':
            n_anagwghs = 85 / 4.9
        elif z == 'Normal_filtered2':
            n_anagwghs = 8500 / 10.1
        elif z == 'Ruttelprogram_filtered2':
            n_anagwghs = 20000 / 14.3
    #Damage, fatigue life for every measurement point

    #Damage_arrays = {'MP49' : 0.0, 'MP51' : 0.0, 'MP52' : 0.0, 'MP55' : 0.0, 'MP58' : 0.0}

    #Loop for every measurement point

    # final_array = {}
    # kk=0
    #
    # N_TELIKO = np.zeros(5)
        N_eq_dict = {}
        S_dict = {}
        Mean_dict = {}
        for i in measurements_points:
            #find peaks and valleys using fatpack
            peaks_and_valleys, peaks_and_valleys_index = fatpack.find_reversals(df_new[i].to_numpy())
            #find Δσ and the residue (cycles that are not completed)
            cycles, residue = fatpack.find_rainflow_cycles(peaks_and_valleys)
            processed_residue = fatpack.concatenate_reversals(residue, residue)
            cycles_residue, _ = fatpack.find_rainflow_cycles(processed_residue)
            cycles_total = np.concatenate((cycles, cycles_residue))
            ranges = np.abs(cycles_total[:, 1] - cycles_total[:, 0])
            mean_value = (cycles_total[:, 0] + cycles_total[:, 1])/2
            # Plots
            figsize = np.array([340., 140.]) / 25.
            fig = plt.figure(dpi=96, figsize=figsize)
            # # Plotting signal with reversals.
            # ax_signal = plt.subplot2grid((1, 2), (0, 0))
            # ax_signal.plot(df_new[i].to_numpy())
            # ax_signal.plot(peaks_and_valleys_index, df_new[i].to_numpy()[peaks_and_valleys_index], 'ro', fillstyle = 'none',
            # label = 'Reversals')
            # ax_signal.legend()
            # ax_signal.set(title="Signal", ylabel="y", xlabel="Time", xlim=[400, 700])
            # ax_signal.grid(True)

            # Plotting the cumulative distribution of the cycle count
            ax_cumdist = plt.subplot2grid((1, 2), (0, 0))
            N, S = fatpack.find_range_count(ranges, len(cycles_total))
            Ncum = N.sum() - np.cumsum(N)
            ax_cumdist.semilogx(Ncum, S)
            ax_cumdist.set(title="Cumulative distribution, rainflow ranges",
                           xlabel="Count, N", ylabel="Range, S")
            ax_cumdist.grid(True)
            fig.tight_layout()

            N_eq = np.array(N)*n_anagwghs
            Ncum_eq = N_eq.sum() - np.cumsum(N_eq)
            ax_cumdist = plt.subplot2grid((1, 2), (0, 1))
            ax_cumdist.semilogx(Ncum_eq, S)
            ax_cumdist.set(title="Equivalent cumulative distribution, rainflow ranges",
                           xlabel="Cumulative count, N", ylabel="Range, S")
            ax_cumdist.grid(True)
            fig.tight_layout()
            # plt.show(block=False)
            N_eq_dict[i] = N_eq
            S_dict[i] = S
            Mean_dict[i] = mean_value
        N_eq_big_dic[z] = N_eq_dict
        S_big_dic[z] = S_dict
        Mean_big_dic[z] = Mean_dict

    # Plot the cumulative distribution by combining all trials for each measurement point

    # Plot the cumulative distribution by combining all trials for each measurement point
    figsize = np.array([340., 140.]) / 25.

    measurement_points = ['MP49', 'MP51', 'MP52', 'MP55', 'MP58']
    for mp in measurement_points:
        # Combine data for the specific measurement point across all trials
        N_combined = np.hstack([N_eq_big_dic["Flucht_filtered2"][mp],
                                N_eq_big_dic["Normal_filtered2"][mp],
                                N_eq_big_dic["Ruttelprogram_filtered2"][mp]])
        S_combined = np.hstack([S_big_dic["Flucht_filtered2"][mp],
                                S_big_dic["Normal_filtered2"][mp],
                                S_big_dic["Ruttelprogram_filtered2"][mp]])
        Mean_combined = np.hstack([Mean_big_dic["Flucht_filtered2"][mp],
                                Mean_big_dic["Normal_filtered2"][mp],
                                Mean_big_dic["Ruttelprogram_filtered2"][mp]])

        # Sort by stress range to make the cumulative count correspond to increasing stress
        sorted_indices = np.argsort(S_combined)
        S_combined = S_combined[sorted_indices]
        N_combined = N_combined[sorted_indices]
        Sa_combined = np.array(S_combined)/2

        # Calculate cumulative sum from the end for a cumulative distribution plot
        Ncum_combined = np.cumsum(N_combined[::-1])[::-1]

        # Plotting the combined cumulative distribution
        fig = plt.figure(dpi=96, figsize=figsize)
        ax_cumdist = plt.subplot()
        ax_cumdist.semilogx(Ncum_combined, S_combined, color='red', label=f'Cumulative (P1 + P2 + P3) - {mp}')
        ax_cumdist.set(title=f"Combined Cumulative Distribution for {mp}",
                       xlabel="Cumulative count, N", ylabel="Range, S")
        ax_cumdist.grid(True)
        ax_cumdist.legend()
        fig.tight_layout()
        # plt.show(block=False)
        N = []
        D_m = []
        for w in range(len(N_combined)):
            sa_new = (Sa_combined[w] + M * Mean_combined[w])
            N.append(10 ** (np.log10(sa_new / sf_) / b_ - np.log10(2)))
            # Miner
            D_m.append(N_combined[w] / N[w])

        D_m_sum = sum(D_m)
        N_T[w] = 1 / D_m_sum

    print('------THEORETICAL------')
    print('Cycles to failure for each measurement point:')
    print(N_T)

    fea = pd.read_csv('fea-res.csv')
    n_fea = np.array(fea.iloc[1:18,10]).astype(float)

    mean_fea49 = np.array(fea.iloc[1:18, 21]).astype(float) + 0.05
    semiamp_fea49 = np.array(fea.iloc[1:18,26]).astype(float) - 0.05

    mean_fea51 = np.array(fea.iloc[1:18, 22]).astype(float) + 0.05
    semiamp_fea51 = np.array(fea.iloc[1:18, 27]).astype(float) - 0.05

    mean_fea52 = np.array(fea.iloc[1:18, 23]).astype(float) + 0.05
    semiamp_fea52 = np.array(fea.iloc[1:18, 28]).astype(float) - 0.05

    mean_fea55 = np.array(fea.iloc[1:18, 24]).astype(float) + 0.05
    semiamp_fea55 = np.array(fea.iloc[1:18, 29]).astype(float) - 0.05

    mean_fea58 = np.array(fea.iloc[1:18, 25]).astype(float) + 0.05
    semiamp_fea58 = np.array(fea.iloc[1:18, 30]).astype(float) - 0.05

    n49fea = []
    d49fea = []
    n51fea = []
    d51fea = []
    n52fea = []
    d52fea = []
    n55fea = []
    d55fea = []
    n58fea = []
    d58fea = []
    for e in range(len(n_fea)):
        sa_fea49 = (semiamp_fea49[e] + M * mean_fea49[e])
        n49fea.append(10 ** (np.log10(sa_fea49 / sf_) / b_ - np.log10(2)))

        d49fea.append(n_fea[e] / n49fea[e])

        sa_fea51 = (semiamp_fea51[e] + M * mean_fea51[e])
        n51fea.append(10 ** (np.log10(sa_fea51 / sf_) / b_ - np.log10(2)))

        d51fea.append(n_fea[e] / n51fea[e])

        sa_fea52 = (semiamp_fea52[e] + M * mean_fea52[e])
        n52fea.append(10 ** (np.log10(sa_fea52 / sf_) / b_ - np.log10(2)))

        d52fea.append(n_fea[e] / n52fea[e])

        sa_fea55 = (semiamp_fea55[e] + M * mean_fea55[e])
        n55fea.append(10 ** (np.log10(sa_fea55 / sf_) / b_ - np.log10(2)))

        d55fea.append(n_fea[e] / n55fea[e])

        sa_fea58 = (semiamp_fea58[e] + M * mean_fea58[e])
        n58fea.append(10 ** (np.log10(sa_fea58 / sf_) / b_ - np.log10(2)))

        d58fea.append(n_fea[e] / n58fea[e])



    d49fea_sum = sum(d49fea)
    cycles49_fea = 1 / d49fea_sum
    d51fea_sum = sum(d51fea)
    cycles51_fea = 1 / d51fea_sum
    d52fea_sum = sum(d52fea)
    cycles52_fea = 1 / d52fea_sum
    d55fea_sum = sum(d55fea)
    cycles55_fea = 1 / d55fea_sum
    d58fea_sum = sum(d58fea)
    cycles58_fea = 1 / d58fea_sum

    FEA_failure = np.array([cycles49_fea, cycles51_fea, cycles52_fea, cycles55_fea, cycles58_fea])
    print('-----FEA------')
    print('Cycles to failure for each measurement point:')
    print(FEA_failure/300)

    damage_max = max(d55fea)
    damage_max_idx = d55fea.index(max(d55fea))
    nmax_from_damage = 1 / damage_max


    print('----BLOCK----')
    print("Maximum cycles for Block " + str(damage_max_idx + 1) + " are " + str(n55fea[damage_max_idx]/300000))



if __name__ == "__main__":
    main()
