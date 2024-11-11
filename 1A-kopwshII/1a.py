import numpy as np
import matplotlib.pyplot as plt
import fatpack
import pandas as pd
import utilities as utils

#********************************CODE BEGINS HERE********************************************

def main():

    M = 0.4
    k = 6.5
    N_1 = 500000
    N_2 = 100000
    Sm_kati = 800
    sa_1_kati = 480
    sa_2_kati = sa_1_kati * (N_2/N_1) ** (-1/k)
    sa_1_kati_allo = sa_1_kati + M * Sm_kati
    sa_2_kati_allo = sa_2_kati + M * Sm_kati
    b_ = np.log10(sa_1_kati_allo/sa_2_kati_allo)/np.log10(N_1/N_2)
    sf_ = sa_1_kati_allo/((2*N_1) ** b_)
    filename = ('Flucht_filtered2', 'Normal_filtered2', 'Ruttelprogram_filtered2')
    Counts_T = {}
    Sa_TT = {}
    Sm_TT = {}
    N_T = {}


    #Damage, fatigue life for every measurement point

    Damage_arrays = {'MP49' : 0.0, 'MP51' : 0.0, 'MP52' : 0.0, 'MP55' : 0.0, 'MP58' : 0.0}

    #Loop for every measurement point

    for i in Damage_arrays:
        Counts = []
        Sa_T = []
        Sm_T = []
        C_T = []
        for z in filename:

            # _________________________INPUT DATA SECTION________________________________________________________
            # *Everything you will need to change are located in this section
            # You will change the filename of you want to read

            # __________________________________________________________________________________________________

            # Read the file of the filtered data
            df_new = utils.reader(z)

            # Exclude the time column
            measurements_points = df_new.keys()[2:]
            if z == 'Flucht_filtered2':
                n_anagwghs = 85 / 4.9
            elif z == 'Normal_filtered2':
                n_anagwghs = 8500 / 10.1
            elif z == 'Ruttelprogram_filtered2':
                n_anagwghs = 20000 / 14.3

            #find peaks and valleys using fatpack
            peaks_and_valleys, peaks_and_valleys_index = fatpack.find_reversals(df_new[i].to_numpy())
            #find Δσ and the residue (cycles that are not completed)
            cycles, residue = fatpack.find_rainflow_cycles(peaks_and_valleys)
            processed_residue = fatpack.concatenate_reversals(residue, residue)
            cycles_residue, _ = fatpack.find_rainflow_cycles(processed_residue)
            cycles_total = np.concatenate((cycles, cycles_residue))
            ranges = np.abs(cycles_total[:, 1] - cycles_total[:, 0])
            Sa = ranges/2
            Sm = (cycles_total[:, 0] + cycles_total[:, 1])/2
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

            # Plotting the cumulative distribution of the cycle count ########
            # ax_cumdist = plt.subplot2grid((1, 2), (0, 0))
            C, S = fatpack.find_range_count(Sa, len(cycles_total))

            # mean_stresses = np.array([(cycle[0] + cycle[1]) / 2 for cycle in cycles_total])

            # 6. Calculate the mean stress corresponding to each range (C, S)
            # For each unique range (S), get the mean of corresponding mean_stresses
            # unique_ranges = np.unique(S)
            # mean_stress_per_range = {}
            # Sm = []
            # for s_range in unique_ranges:
            #     indices = np.where(S == s_range)[0]  # Find indices where S matches this range
            #     Sm.append (np.mean(mean_stresses[indices]))
            # print(len(mean_stress_per_range),len(C))
            # mean_stress_per_range now contains the mean stresses for each range (S)
            # print("Mean stress per range (S):", mean_stress_per_range)

            Ccum = C.sum() - np.cumsum(C) ###############
            # ax_cumdist.semilogx(Ccum, S)
            # ax_cumdist.set(title="Cumulative distribution, rainflow ranges",
            #                xlabel="Count, N", ylabel="Range, S")
            # ax_cumdist.grid(True)
            # fig.tight_layout()
            ######### OPEN FOR [cycle count, range, s1, s2] *s1, s2: are s_min and s_max in random order ########
            #print(i)
            #print(np.hstack((np.column_stack((np.array(Ncum).reshape(-1,1), np.array(S).reshape(-1,1))), cycles_total)))

            ######## OPEN FOR [cycle count, range, mean value of cycle] ##############
            # print(i)
            # print(np.column_stack((np.column_stack((np.array(Ncum).reshape(-1, 1), np.array(S).reshape(-1, 1))), mean_value)))

            ######## OPEN FOR .csv file generation -*one for each i(accelerometer ##########
            # final_array[i] = np.column_stack((np.column_stack((np.array(Ncum).reshape(-1, 1), np.array(S).reshape(-1, 1))), mean_value))
            # utils.writer(final_array[i], filename + i + '_rainflow_results.csv')

            C = np.array(C)*n_anagwghs
            Ccum_eq = C.sum() - np.cumsum(C)
            # ax_cumdist = plt.subplot2grid((1, 2), (0, 1)) ##########
            # ax_cumdist.semilogx(Ccum_eq, S)
            # ax_cumdist.set(title="Equivalent cumulative distribution, rainflow ranges",
            #                xlabel="Cumulative count, N", ylabel="Range, S")
            # ax_cumdist.grid(True)
            # fig.tight_layout()
            # plt.show(block=False)
            # Counts = np.hstack((Counts,C))
            Sa_T = np.hstack((Sa_T, S))
            Sm_T = np.hstack((Sm_T, Sm))
            C_T = np.hstack((C_T, C))
        N = []
        D_m =[]
        for j in range(len(C_T)):
            sa_new = (Sa_T[j] + M * Sm_T[j])
            N.append (10 ** (np.log10(sa_new / sf_) / b_ - np.log10(2)))
            # Miner
            D_m.append (C_T[j]/N[j])

        D_m_sum = sum(D_m)
        N_T[i] = 1/D_m_sum

        # Counts_T[i] = Counts
        Sa_TT[i] = Sa_T
        Sm_TT[i] = Sm_T
    print(N_T)

if __name__ == "__main__":
    main()