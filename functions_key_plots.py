def jacknifing_hosts_ids(table, x):
    #takes table that contains host ids and generates and array of jacknifed ids
    #returns: an array of arrays of hostids
    #table = the table containing the hosts you want to jacknife
    #x = the size of the interval you want to jacknife
    host_arr = []
    unit = int(len(table['HOSTID'])/x)
    unit_og = int(len(table['HOSTID'])/x)
    j = 0
    hosts = []
    while j<=unit:
        if j >= len(table['HOSTID']):
            break
        hosts.append(table['HOSTID'][j])
        if j == unit:
            j=unit
            unit = unit+unit_og
            host_arr.append(hosts)
            hosts = [] 
        j = j+1
    return host_arr

def luminosity_distribution_sims(lums, rads, rprojs):
    #this function determines the luminosity distribution of the brightest satellites within different radial cuts
    #returns: a dictionary of arrays
    #lums = the luminosities for all the simulated satellites, a dictinoary of arrays
    #rads = their projected radii, a dictinoary of arrays
    #rprojs = an array of radial cuts
    lum_rad = {}
    for key in lums.keys():
        lum_rad[key] = []
        randnum = np.random.randint(len(lums[key]))
        
        for i in range(len(rprojs)):
            rad_bound = rprojs[i]
            bright = np.min(lums[key][randnum])
            ind = list(lums[key][randnum]).index(np.min(lums[key][randnum])) 
            rad_sat = rads[key][randnum][ind]
            if rad_sat<=rad_bound:
                
                lum_rad[key].append([rad_bound, bright])
            else:
                sat_mag = 0
                for j in range(len(lums[key][randnum])):
                    if rads[key][randnum][j]<=rad_bound and lums[key][randnum][j]<=sat_mag:
                        sat_mag = lums[key][randnum][j]
                lum_rad[key].append([rad_bound, sat_mag])
    return lum_rad

def occ_rate_luminosities_with_radial_cut(index, rad, lum_arr):
    #finds the occurence rate across radial cuts of satellites of a certain luminosity
    #returns: a dictionary
    #index = your index of choice which determines what radial cut you'll be looking at
    #rad = the data structure your trying to determine the occurrence rate of
    #lum_arr = an array of luminosities to compare against
    occs = {}
    for i in range(len(lum_arr)):
        num_true = 0
        for key in rad.keys():
            if rad[key][index][1]<lum_arr[i]:
                num_true+=1
        occs[lum_arr[i]] = num_true/len(rad.keys())
    return occs
def occ_rate_mg_with_radial_cut(index, mg_rad, mg_arr):
    #finds the occurence rate across radial cuts of satellites of a certain luminosity
    #returns: a dictionary
    #index = your index of choice which determines what radial cut you'll be looking at
    #mg_rad = the data structure your trying to determine the occurrence rate of
    #mg_arr = an array of gaps to compare against
    mg_occs = {}
    for i in range(len(mg_arr)):
        num_true = 0
        for key in mg_rad.keys():
            if mg_rad[key][index][1]>mg_arr[i] and mg_rad[key][index][1]!=0:
                num_true+=1
        mg_occs[mg_arr[i]] = num_true/len(mg_rad.keys())
    return mg_occs

def LMCs_sat_host_rad(sorted_dictionary_satellites, host_list, rad_arr):
    #determines the brightest satellite within a radial cut across some set of SAGA hosts
    #returns: a dictionary of arrays
    #sored_dictionary_satellites = dictionary of satellites sorted by hostids containing information on those satellites
    #host_list = list of host ids by which to sort through the satellites
    #rad_arr = array of radii by which to cut
    hosts_LMCs = {}
    for i in range(len(host_list)):
        key = host_list[i]
        hosts_LMCs[key] = []
        if key in sorted_dictionary_satellites.keys():
            for i in range(len(rad_arr)):
                rad_bound = rad_arr[i]
                rad_sat = sorted_dictionary_satellites[key]['Mr'== np.min(sorted_dictionary_satellites[key]["Mr"])]['RHOST_KPC']
                
                if rad_sat<=rad_bound:
                    min_mag = np.sort(sorted_dictionary_satellites[key]["Mr"])[0]
                    hosts_LMCs[key].append([rad_bound,min_mag])
                else:
                    sat_mag = 0
                    for j in range(len(sorted_dictionary_satellites[key])):
                        if sorted_dictionary_satellites[key]['RHOST_KPC'][j] <= rad_arr[i] and sorted_dictionary_satellites[key]['Mr'][j]<sat_mag:
                            sat_mag = sorted_dictionary_satellites[key]['Mr'][j]
                    hosts_LMCs[key].append([rad_bound, sat_mag])
        else: 
            for i in range(len(rad_arr)):
                rad_bound = rad_arr[i]
                hosts_LMCs[key].append([rad_bound, 0])
    return hosts_LMCs

def magnitude_gap_sat_host_rad(sorted_dictionary_satellites, host_list, rad_arr):
    #determines the smallest magnitude gap between host and satellite within a radial cut across some set of SAGA hosts
    #returns: a dictionary of arrays
    #sored_dictionary_satellites = dictionary of satellites sorted by hostids containing information on those satellites
    #host_list = list of host ids by which to sort through the satellites
    #rad_arr = array of radii by which to cut
    mag_gap = {}
    for i in range(len(host_list)):
        key = host_list[i]
        mag_gap[key] = []
        if key in sorted_dictionary_satellites.keys():
            for i in range(len(rad_arr)):
                rad_bound = rad_arr[i]
                rad_sat = sorted_dictionary_satellites[key]['Mr'== np.min(sorted_dictionary_satellites[key]["Mr"])]['RHOST_KPC']
                host_mag = sorted_dictionary_hosts[key]["K_ABS"][0]
                if rad_sat<=rad_bound:
                    min_mag = np.sort(sorted_dictionary_satellites[key]["Mr"])[0]
                    gap = host_mag-min_mag
                    mag_gap[key].append([rad_bound, gap])
                else:
                    sat_mag = 0
                    for j in range(len(sorted_dictionary_satellites[key])):
                        if sorted_dictionary_satellites[key]['RHOST_KPC'][j] <= rad_arr[i] and sorted_dictionary_satellites[key]['Mr'][j]<sat_mag:
                            sat_mag = sorted_dictionary_satellites[key]['Mr'][j]
                    mag_gap[key].append([rad_bound, host_mag-sat_mag])
        else: 
            for i in range(len(rad_arr)):
                rad_bound = rad_arr[i]
                mag_gap[key].append([rad_bound, 0])
    return mag_gap

def brightest_sat_mg_rad_sims_corr(lums, rads, rprojs, hosts_lum):
    #this function determines the magnitude gap distribution of the brightest satellites with their hosts within different radial cuts
    #returns: a dictionary of arrays
    #lums = the luminosities for all the simulated satellites, a dictinoary of arrays
    #rads = their projected radii, a dictionary of arrays
    #rprojs = an array of radial cuts
    #hosts_lum = predicted host luminosities from the simulated data
    mag_gap_rad = {}
    for key in lums.keys():
        mag_gap_rad[key] = []
        randnum = np.random.randint(len(lums[key]))
        host_lum = hosts_lum[key][randnum]
        for i in range(len(rprojs)):
            rad_bound = rprojs[i]
            bright = np.min(lums[key][randnum]) - 0.5
            ind = list(lums[key][randnum]).index(np.min(lums[key][randnum])) 
            rad_sat = rads[key][randnum][ind]
            if rad_sat<=rad_bound:
                gap = host_lum-bright
                mag_gap_rad[key].append([rad_bound, gap])
            else:
                sat_mag = 0
                for j in range(len(lums[key][randnum])):
                    if rads[key][randnum][j]<=rad_bound and lums[key][randnum][j]<=sat_mag:
                        sat_mag = lums[key][randnum][j] -0.5
                mag_gap_rad[key].append([rad_bound, host_lum-sat_mag])
    return mag_gap_rad

def find_spirals_ellipticals(spirals, ellipticals):
    #determines list of SAGA spiral galaxies that are within 0.05 luminosity of each elliptical galaxy
    #returns: dictionary of array, keys are ellipticals, arrays of matched spirals
    #spirals = dataframe of SAGA spiral galaxies
    #ellipticals = dataframe of SAGA elliptical galaxies
    ell_lums = []
    ell_matches = {}
    for i in range(len(ellipticals)):
        ell_matches[ellipticals['HOSTID'][i]] = []
        ell_lums.append(ellipticals['K_ABS'][i])
    for j in range(len(ell_lums)):
        for k in range(len(spirals)):
            if np.abs(spirals['K_ABS'][k] - ell_lums[j])<= 0.05:
                ell_matches[ellipticals['HOSTID'][j]].append(spirals['HOSTID'][k])
    return ell_matches

def select_spirals(matches):
    #takes the dicitonary created by find_spirals_ellipticals and selects one matched spiral for each elliptical
    #returns: a dictionary with ellipticals as keys and one spiral as the value
    single_matched = {}
    for key in matches.keys():
        for i in range(len(matches[key])):
            if test_matches[key][i] not in single_matched.values():
                single_matched[key] = matches[key][i]
                break
    return single_matched

def determining_host_count_arr_rproj(table, sorted_dictionary_satellites, rprojs, Mr):
    #determines number of hosts in the data that have a satellite brighter than some luminosity within a radial cut
    #returns: dictionary of arrays that contain the radial cut and number of hosts with a satellite fo such brithness within the cut
    #table = table of hosts with hostids
    #sorted_dictionary_satellites = dictionary of satellites sorted by hostids containing information on those satellites
    #rprojs = an array of radial cuts
    #Mr = whatever luminosity cutoff you choose for the satellites
    hosts_LMCs = {}
    for key in table['HOSTID']:
        hosts_LMCs[key] = []
        for i in range(len(rprojs)):
            count = 0
            if key in sorted_dictionary_satellites.keys():
                for j in range(len(sorted_dictionary_satellites[key])):
                    if sorted_dictionary_satellites[key]["Mr"][j] <= Mr and sorted_dictionary_satellites[key]["RHOST_KPC"][j] <= rprojs[i]:
                        count = count+1
            arr = [rprojs[i], count]
            hosts_LMCs[key].append(arr) 
    return hosts_LMCs

def determining_host_count_arr_rproj_sims(lums, rads, rprojs, Mr):
    #determines number of hosts in the data that have a satellite brighter than some luminosity within a radial cut
    #returns: dictionary of arrays that contain the radial cut and number of hosts with a satellite for such brightness within the cut
    #lums = the luminosities for all the simulated satellites, a dictinoary of arrays
    #rads = their projected radii, a dictinoary of arrays
    #rprojs = an array of radial cuts
    #Mr = whatever luminosity cutoff you choose for the satellites
    hosts_LMCs = {}
    for key in lums.keys():
        hosts_LMCs[key] = []
        randnum = np.random.randint(len(lums[key]))
        for i in range(len(rprojs)):
            count = 0
            for j in range(len(lums[key][randnum])):
                if lums[key][randnum][j] <= Mr and rads[key][randnum][j] <= rprojs[i]:
                        count = count+1
            arr = [rprojs[i], count]
            hosts_LMCs[key].append(arr) 
    return hosts_LMCs

def ratio_sims_data(sims, data):
    #takes a ratio of simulated values to data values
    #returns: dictionary of that ratio
    #sims = some dictionary of numerical values
    #data = some dictionary of numerical values
    ratio = {}
    for key in sims.keys():
        ratio[key] = sims[key]/data[key]
    return ratio

def determining_mg_count(table, sorted_dictionary_satellites, mgs):
    #determines number of hosts in the data that have a magnitude gap smaller than some gap within a radial cut
    #returns: dictionary of arrays that contain the radial cut and number of hosts with a magnitude gap of such within said cut
    #table = table of hosts with hostids
    #sorted_dictionary_satellites = dictionary of satellites sorted by hostids containing information on those satellites
    #mg = an array of magnitude gap cuts
    mg_count = {}
    for i in range(len(mgs)):
        counts = []
        for key in table['HOSTID']:
            counter = 0
            if key in sorted_dictionary_satellites.keys():
                for j in range(len(sorted_dictionary_satellites[key])):
                    gap = sorted_dictionary_satellites[key]['Mr'][j] - table['K_ABS'][table['HOSTID']==key]
                    if gap<mgs[i]:
                        counter = counter+1
            counts.append(counter)
        arr = [np.mean(counts), np.std(counts)]
        mg_count[mg_arr[i]] = arr
    return mg_count

def mean_std_lum_rad(lmc_dict):
    #determines the means and standard deviations across a dictionary of arrays
    #returns: two dictionaries of arrays
    #lmc_dict = some dictionary of arrays created from data measurements
    means = {}
    stds = {}
    for i in range(len(lmc_dict[list(lmc_dict.keys())[0]])):
        tot = []
        for key in lmc_dict.keys():
            tot.append(lmc_dict[key][i][1])
        means[lmc_dict[key][i][0]] = np.mean(tot)
        stds[lmc_dict[key][i][0]] = np.std(tot)
    return means, stds

def mean_std_sims(occs):
    #determines the means and standard deviations across a dictionary of arrays
    #returns: two dictionaries of arrays
    #occs = some dictionary of arrays created from simulated measurements
    means = {}
    stds = {}
    for i in occs[0].keys():
        tots = []
        for j in range(len(occs)):
            tots.append(occs[j][i])
        means[i] = np.mean(tots)
        stds[i] = np.std(tots)
    return means, stds

def just_mean_std_mg_count(mean_std):
    #separates means and standard deviations from a dctionary containing those means and deviations 
    #returns: two arrays
    #returns: two dicitonarys, one of means and one of standard deviations 
    means = []
    stds = []
    for key in mean_std.keys():
        means.append(mean_std[key][0])
        stds.append(mean_std[key][1])
    return means, stds