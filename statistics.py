
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats
import statsmodels.api as sm

import global_var
from get_ip2as_from_bdrmapit import ConnectToBdrMapItDb, GetIp2ASFromBdrMapItDb,CloseBdrMapItDb, InitBdrCache, \
                                    ConstrBdrCache

def CheckUndo():    
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    #syd-au_20160415
    for vp in global_var.vps:   #vp
        for year in range(2016,2021):
            for month in range(1,13):
                if (year == 2016 and month < 4) or (year == 2020 and month > 4): #2016年4月前peeringdb数据不准，2020年5月后的数据不全
                    continue
                month_str = str(month).zfill(2)              
                date = str(year) + month_str + '15'   #date
                for method in global_var.map_methods:   #method
                    #syd-au_20160415/ribs/record_syd-au.20160415_ribs
                    filepath = vp + '_' + date + '/' + method + '/record_' + vp + '.' + date + '_' + method
                    #print(filepath)
                    if not os.path.isfile(filepath):
                        print('NOTE1!' + vp + date + method)
                        continue
                    with open(filepath, 'r') as f:
                        res = f.read()
                    re_res = re.findall('In ab_filter1, ab num: (\d.*), ab precent: (\d\.\d.*)', res)
                    if not re_res:
                        print('NOTE2!' + vp + date + method)
                    elif float(re_res[0][1]) == 0.0:
                        print('NOTE3!' + vp + date + method)

def CollectMatchIpStat():
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    
    for vp in global_var.vps:   #vp
        for method in global_var.map_methods:   #method
            os.system("cat %s*/%s/record2_* > collect_record2_%s_%s" %(vp, method, vp, method))

def StatNobgp():
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    wr_dir = 'statistics'
    if not os.path.isdir(wr_dir):
        os.makedirs(wr_dir)
    wf = open(wr_dir + '/no_bgp_stat', 'w')
    wf.write("#format: vp_date: no_bgp_percent total_num\n")
    for root,dirs,files in os.walk(path):
        for cur_dir in dirs:
            #if cur_dir == wr_dir:
            if cur_dir != 'nrt-jp_20180415':
                continue
            res = os.popen('wc -l ' + cur_dir + '/ribs/0_no_bgp').readline()
            #print(res) #'89029 nrt-jp_20180415/ribs/0_no_bgp'
            no_bgp_num = int(res.split(' ')[0])
            #print(no_bgp_num)
            res = os.popen("grep -r \"Total\" %s/ribs/record*" %cur_dir).readline()
            #print(res)
            total_num = int(res.split(':')[1].strip('\n').strip(' '))
            #print(total_num)
            wf.write("%s: %f %d\n" %(cur_dir, no_bgp_num / total_num, total_num))
    wf.close()

def StatIxp():
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    wr_dir = 'statistics'
    if not os.path.isdir(wr_dir):
        os.makedirs(wr_dir)
    wf = open(wr_dir + '/ixp', 'w')
    wf.write("#format: vp_date: ixp_percent total_num\n")
    for root,dirs,files in os.walk(path):
        for cur_dir in dirs:
            #if cur_dir == wr_dir:
            if cur_dir != 'nrt-jp_20180415':
                continue
            tmp = cur_dir.replace('_', '.')
            with open("%s/ribs/record_%s_ribs" %(cur_dir, tmp), 'r') as f:
                res = f.read()
            re_res = re.findall("Total valid trace num: (\d.*?)\n", res, re.DOTALL)
            total_num = int(re_res[0])
            re_res = re.findall("ixp num: (\d.*?)\n", res, re.DOTALL)
            ixp_num = int(re_res[0])
            wf.write("%s: %f %d\n" %(cur_dir, ixp_num / total_num, total_num))
            print("%s: %f %d" %(cur_dir, ixp_num / total_num, total_num))
    wf.close()

def StatMultipath():
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    wr_dir = 'statistics'
    if not os.path.isdir(wr_dir):
        os.makedirs(wr_dir)
    wf = open(wr_dir + '/multipath', 'w')
    wf.write("#format: vp_date: multipath_percent total_num\n")
    for root,dirs,files in os.walk('.'):
        for cur_dir in dirs:
            if cur_dir == wr_dir:
            #if cur_dir != 'nrt-jp_20180415':
                continue
            for sub_root,sub_dirs,sub_files in os.walk(cur_dir):
                for cur_method in sub_dirs:
                    tmp = cur_dir.replace('_', '.')
                    if not os.path.isfile("%s/%s/record_%s_%s" %(cur_dir, cur_method, tmp, cur_method)):
                        continue
                    with open("%s/%s/record_%s_%s" %(cur_dir, cur_method, tmp, cur_method), 'r') as f:
                        res = f.read()
                    re_res = re.findall("Total valid trace num: (\d.*?)\n", res, re.DOTALL)
                    total_num = 0
                    if re_res:
                        total_num = int(re_res[0])
                    multipath_num = 0
                    re_res = re.findall("still multi num: (\d.*?)\n", res, re.DOTALL)
                    if re_res:
                        multipath_num = int(re_res[0])
                    else:
                        re_res = re.findall("multi num: (\d.*?)\n", res, re.DOTALL)
                        if re_res:
                            multipath_num = int(re_res[0])
                    if total_num == 0:
                        wf.write("%s_%s: 0 0\n" %(cur_dir, cur_method))
                        print("%s_%s: 0 0" %(cur_dir, cur_method))
                    else:
                        wf.write("%s: %f %d\n" %(cur_dir, multipath_num / total_num, total_num))
                        #print("%s: %f %d" %(cur_dir, multipath_num / total_num, total_num))
    wf.close()

def StatAb():
    path = global_var.par_path + global_var.out_my_anatrace_dir
    os.chdir(path)
    wr_dir = 'statistics'
    if not os.path.isdir(wr_dir):
        os.makedirs(wr_dir)
    wf = open(wr_dir + '/ab', 'w')
    wf.write("#format: vp_date_method: ab_0_percent ab_1_percent total_num\n")
    for root,dirs,files in os.walk('.'):
        for cur_dir in dirs:
            if cur_dir == wr_dir:
            #if cur_dir != 'nrt-jp_20180415':
                continue
            for sub_root,sub_dirs,sub_files in os.walk(cur_dir):
                for cur_method in sub_dirs:
                    tmp = cur_dir.replace('_', '.')
                    if not os.path.isfile("%s/%s/record_%s_%s" %(cur_dir, cur_method, tmp, cur_method)):
                        continue
                    with open("%s/%s/record_%s_%s" %(cur_dir, cur_method, tmp, cur_method), 'r') as f:
                        res = f.read()
                    re_res = re.findall("Total valid trace num: (\d.*?)\n", res, re.DOTALL)
                    total_num = 0
                    if re_res:
                        total_num = int(re_res[0])
                    ab_1 = 0
                    ab_0 = 0
                    re_res = re.findall('In ab_filter1, ab num: (\d.*), ab precent: (0\.\d.*)', res)
                    if re_res:
                        ab_1 = int(re_res[0][0])
                        #print(re_res[0][0])
                        ab_0 = int(int(re_res[0][0]) / float(re_res[0][1]))
                        #print(re_res[0][1])
                    if total_num == 0:
                        wf.write("%s_%s: 0 0 0\n" %(cur_dir, cur_method))
                        print("%s_%s: 0 0 0" %(cur_dir, cur_method))
                    else:
                        wf.write("%s_%s: %f %f %d\n" %(cur_dir, cur_method, ab_0 / total_num, ab_1 / total_num, total_num))
                        print("%s_%s: %f %f %d" %(cur_dir, cur_method, ab_0 / total_num, ab_1 / total_num, total_num))
    wf.close()

def CalDateIndex(date):
    year = int(date[0:4])
    month = int(date[4:6])
    return ((year - 2016) * 12 + month - 4)

#vps = ['nrt-jp', 'per-au', 'syd-au', 'zrh2-ch']
def PlotAbStat():
    #2016.4~2020.4, four years + 1 month, 49 months, thus each vp has 49 time-plot res
    res = dict()
    time_plots = 49
    for vp in global_var.vps:
        res[vp] = dict()#[dict() for i in range(0, 2)] #ab_0, ab_1
        for method in global_var.map_methods:
            res[vp][method] = [] #每个method有ab_0, ab_1两种结果
            for i in range(0, 2): #每个ab_i有49个时间点结果
                res[vp][method].append([0.0 for j in range(0, time_plots)])
    rf = open(global_var.par_path + global_var.out_my_anatrace_dir + '/statistics/ab', 'r')
    curline = rf.readline()
    while curline:
        if curline.startswith('#'):
            curline = rf.readline()
            continue
        #format: vp_date_method: ab_0_percent ab_1_percent total_num
        elems = curline.split(' ')
        sub_elems = elems[0].strip(':').split('_')
        vp = sub_elems[0]
        index = CalDateIndex(sub_elems[1])
        method = sub_elems[2]
        ab_0 = float(elems[1])
        ab_1 = float(elems[2])
        #print(index)
        res[vp][method][0][index] = ab_0
        res[vp][method][1][index] = ab_1
        curline = rf.readline()
    
    color_list = ['#000000', '#0000FF', '#8A2BE2', '#A52A2A', '#DEB887', '#5F9EA0', '#7FFF00', '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C', '#00FFFF', '#00008B', '#008B8B', '#B8860B', '#A9A9A9', '#006400', '#BDB76B', '#8B008B', '#556B2F', '#FF8C00', '#9932CC', '#8B0000', '#E9967A', '#8FBC8F', '#483D8B', '#2F4F4F', '#00CED1', '#9400D3', '#FF1493', '#00BFFF', '#696969', '#1E90FF', '#B22222', '#FFFAF0', '#228B22', '#FF00FF', '#DCDCDC', '#F8F8FF', '#FFD700', '#DAA520', '#808080', '#008000', '#ADFF2F', '#F0FFF0', '#FF69B4', '#CD5C5C', '#4B0082', '#FFFFF0', '#F0E68C', '#E6E6FA', '#FFF0F5', '#7CFC00', '#FFFACD', '#ADD8E6', '#F08080', '#E0FFFF', '#FAFAD2', '#90EE90', '#D3D3D3', '#FFB6C1', '#FFA07A', '#20B2AA', '#87CEFA', '#778899', '#B0C4DE', '#FFFFE0', '#00FF00', '#32CD32', '#FAF0E6', '#FF00FF', '#800000', '#66CDAA', '#0000CD', '#BA55D3', '#9370DB', '#3CB371', '#7B68EE', '#00FA9A', '#48D1CC', '#C71585', '#191970', '#F5FFFA', '#FFE4E1', '#FFE4B5', '#FFDEAD', '#000080', '#FDF5E6', '#808000', '#6B8E23', '#FFA500', '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98', '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F', '#FFC0CB', '#DDA0DD', '#B0E0E6', '#800080', '#FF0000', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#FAA460', '#2E8B57', '#FFF5EE', '#A0522D', '#C0C0C0', '#87CEEB', '#6A5ACD', '#708090', '#FFFAFA', '#00FF7F', '#4682B4', '#D2B48C', '#008080', '#D8BFD8', '#FF6347', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFFFFF', '#F5F5F5', '#FFFF00', '#9ACD32']
    marker_list = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd']
    fig = plt.figure()    
    #设置X轴标签  
    plt.xlabel('date')  
    #设置Y轴标签  
    plt.ylabel('ab_percent') 
    ax = []
    for i in range(0, 4):
        ax.append(fig.add_subplot(221 + i))
    i = 0
    for vp in global_var.vps:
        j = 0
        for method in global_var.map_methods:
            for ab_i in range(0, 2):
                ax[i].scatter(range(0, time_plots), res[vp][method][0], c = color_list[j], marker= marker_list[ab_i], label = method + '_' + str(ab_i))
            j += 1
        i += 1  
    #plt.legend(loc='upper left')
    plt.show()

def CalMatchIpStat():
    path = global_var.par_path + global_var.out_my_anatrace_dir + '/statistics/matchip_stat/'
    os.chdir(path)
    for root,dirs,files in os.walk('.'):
        for filename in files:
            print(filename)
            with open(filename, 'r', encoding='unicode_escape') as f:
                res = f.read()
            re_res = re.findall(", match rate: (.*?), unmatch rate: (.*?), unknown rate: (.*?),", res)
            if not re_res:
                print('Format error 1!')
                continue
            sum_stat = [0 for i in range(3)]
            max_stat = [0 for i in range(3)]
            min_stat = [1 for i in range(3)]
            for elem in re_res:
                for i in range(0, 3):
                    tmp = float(elem[i])
                    sum_stat[i] += tmp
                    if max_stat[i] < tmp:
                        max_stat[i] = tmp
                    if min_stat[i] > tmp:
                        min_stat[i] = tmp
            num1 = len(re_res)
            re_res = re.findall("avg_ip_freq_true: (.*?), avg_ip_freq_false: (.*?), avg_ip_freq_unknown: (.*?)\n", res, re.DOTALL)
            if not re_res:
                print('Format error 2!')
                continue
            sum_freq_stat = [0 for i in range(3)]
            max_freq_stat = [0 for i in range(3)]
            min_freq_stat = [100 for i in range(3)]
            for elem in re_res:
                for i in range(0, 3):
                    tmp = float(elem[i])
                    sum_freq_stat[i] += tmp
                    if max_freq_stat[i] < tmp:
                        max_freq_stat[i] = tmp
                    if min_freq_stat[i] > tmp:
                        min_freq_stat[i] = tmp
            num2 = len(re_res)
            print("\tmatch: %.2f(%.2f, %.2f), %.2f(%.2f, %.2f)" %(sum_stat[0] / num1, max_stat[0], min_stat[0], sum_freq_stat[0] / num2, max_freq_stat[0], min_freq_stat[0]))
            print("\tunmatch: %.2f(%.2f, %.2f), %.2f(%.2f, %.2f)" %(sum_stat[1] / num1, max_stat[1], min_stat[1], sum_freq_stat[1] / num2, max_freq_stat[1], min_freq_stat[1]))
            print("\tunknown: %.2f(%.2f, %.2f), %.2f(%.2f, %.2f)" %(sum_stat[2] / num1, max_stat[2], min_stat[2], sum_freq_stat[2] / num2, max_freq_stat[2], min_freq_stat[2]))

def Tmp():
    record_file_name = '/mountdisk1/ana_c_d_incongruity/out_my_anatrace/nrt-jp_20181015/ribs_midar_bdrmapit/record_nrt-jp.20181015_ribs_midar_bdrmapit'
    with open(record_file_name, 'r') as f:
        res = f.read()
        re_res = re.findall('Total valid trace num: (\d.*)', res)
        if not re_res:
            print('NOTE2!' + vp + date + method)
        else:
            print(int(re_res[0]))

def StatClassify():    
    par_dir = global_var.par_path +  global_var.out_my_anatrace_dir
    os.chdir(par_dir)
    dir_list = os.listdir(par_dir)
    stat_dict = dict()
    classes = ['last_extra', 'first_hop_ab', 'detour', 'bifurc']
    for vp in global_var.vps:
        stat_dict[vp] = dict()
        for cur_class in classes:
            stat_dict[vp][cur_class] = []
        for cur_dir in dir_list:
            if os.path.isdir(os.path.join(par_dir, cur_dir)) and cur_dir.__contains__(vp) and \
            (cur_dir.__contains__('2018') or cur_dir.__contains__('2019')):
                filename = cur_dir + '/ribs_midar_bdrmapit/ana_record_' + cur_dir.replace('_', '.')
                with open(filename, 'r') as rf:
                    data = rf.read()
                re_res = re.findall('last_extra num: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['last_extra'].append(float(re_res[0][1]))
                re_res = re.findall('first_hop_ab num: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['first_hop_ab'].append(float(re_res[0][1]))
                re_res = re.findall('detour num: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['detour'].append(float(re_res[0][1]))
                re_res = re.findall('bifurc num: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['bifurc'].append(float(re_res[0][1]))
    for vp in global_var.vps:
        print(vp)
        for cur_class in classes:
            print("%.2f %.2f %.2f" %(np.mean(stat_dict[vp][cur_class]), np.min(stat_dict[vp][cur_class]), np.max(stat_dict[vp][cur_class])))
     
def StatLastExtraHopPerVP():    
    par_dir = global_var.par_path +  global_var.out_my_anatrace_dir
    os.chdir(par_dir)
    dir_list = os.listdir(par_dir)
    stat_dict = dict()
    classes = ['customer', 'provider', 'peer', 'unknown', 'multi']
    for vp in global_var.vps:
        stat_dict[vp] = dict()
        for cur_class in classes:
            stat_dict[vp][cur_class] = []
        for cur_dir in dir_list:
            if os.path.isdir(os.path.join(par_dir, cur_dir)) and cur_dir.__contains__(vp) and \
            (cur_dir.__contains__('2018') or cur_dir.__contains__('2019')):
                filename = cur_dir + '/ribs_midar_bdrmapit/ana_record_' + cur_dir.replace('_', '.')
                with open(filename, 'r') as rf:
                    data = rf.read()
                re_res = re.findall('reach_dst_last_extra, customer: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['customer'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, provider: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['provider'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, peer: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['peer'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, unknown: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['unknown'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, multi: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict[vp]['multi'].append(float(re_res[0][1]))
    for vp in global_var.vps:
        print(vp)
        for cur_class in classes:
            print("%.2f %.2f %.2f" %(np.mean(stat_dict[vp][cur_class]), np.min(stat_dict[vp][cur_class]), np.max(stat_dict[vp][cur_class])))
                   
def StatLastExtraHop():    
    par_dir = global_var.par_path +  global_var.out_my_anatrace_dir
    os.chdir(par_dir)
    dir_list = os.listdir(par_dir)
    stat_dict = dict()
    classes = ['customer', 'provider', 'peer', 'unknown', 'multi']
    for cur_class in classes:
        stat_dict[cur_class] = []
    for vp in global_var.vps:
        for cur_dir in dir_list:
            if os.path.isdir(os.path.join(par_dir, cur_dir)) and cur_dir.__contains__(vp) and \
            (cur_dir.__contains__('2018') or cur_dir.__contains__('2019')):
                filename = cur_dir + '/ribs_midar_bdrmapit/ana_record_' + cur_dir.replace('_', '.')
                with open(filename, 'r') as rf:
                    data = rf.read()
                re_res = re.findall('reach_dst_last_extra, customer: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict['customer'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, provider: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict['provider'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, peer: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict['peer'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, unknown: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict['unknown'].append(float(re_res[0][1]))
                re_res = re.findall('reach_dst_last_extra, multi: (\d.*), percent: (0\.\d.*)\n', data)
                stat_dict['multi'].append(float(re_res[0][1]))
    for cur_class in classes:
        print("%.2f %.2f %.2f" %(np.mean(stat_dict[cur_class]), np.min(stat_dict[cur_class]), np.max(stat_dict[cur_class])))
         
def StatAbPercent():    
    par_dir = global_var.par_path +  global_var.out_my_anatrace_dir
    os.chdir(par_dir)
    dir_list = os.listdir(par_dir)
    stat_dict = dict()
    classes = ['ab', 'unmap', 'ab_percent', 'unmap_percent']
    for vp in global_var.vps:
        stat_dict[vp] = dict()
        for cur_class in classes:
            stat_dict[vp][cur_class] = []
        for cur_dir in dir_list:
            if os.path.isdir(os.path.join(par_dir, cur_dir)) and cur_dir.__contains__(vp) and \
            (cur_dir.__contains__('2018') or cur_dir.__contains__('2019')):
                filename = cur_dir + '/ribs_midar_bdrmapit/record_' + cur_dir.replace('_', '.') + '_ribs_midar_bdrmapit'
                with open(filename, 'r') as rf:
                    data = rf.read()
                re_res = re.findall('In total, ab num: (\d.*), unmap num: (\d.*), ab precent: (0\.\d.*), unmap precent: (0\.\d.*)\n', data)
                stat_dict[vp]['ab'].append(float(re_res[0][0]))
                stat_dict[vp]['unmap'].append(float(re_res[0][1]))
                stat_dict[vp]['ab_percent'].append(float(re_res[0][2]))
                stat_dict[vp]['unmap_percent'].append(float(re_res[0][3]))
    for vp in global_var.vps:
        print(vp)
        for cur_class in classes:
            print("%.4f %.4f %.4f" %(np.mean(stat_dict[vp][cur_class]), np.min(stat_dict[vp][cur_class]), np.max(stat_dict[vp][cur_class])))

def StatNeigborIp(filename, ip):
    rf = open(filename, 'r')
    curline_trace = rf.readline()
    left_nei = dict()
    right_nei = dict()

    while curline_trace:
        curline_ip = rf.readline()
        #print(curline_ip)
        if curline_ip.__contains__(ip):
            elems = curline_ip.strip('\n').split(']')[1].strip(' ').split(' ')
            index = elems.index(ip)
            if index > 0:
                if elems[index - 1] not in left_nei.keys():
                    left_nei[elems[index - 1]] = [0, '']
                left_nei[elems[index - 1]][0] += 1
            while index < len(elems) and (elems[index] == ip or elems[index].__contains__('<')):
                index += 1
            if index < len(elems):
                if elems[index] not in right_nei.keys():
                    right_nei[elems[index]] = [0, '']
                right_nei[elems[index]][0] += 1
        curline_trace = rf.readline()
    for (key, val) in left_nei.items():
        val[1] = GetIp2ASFromBdrMapItDb(key)
    for (key, val) in right_nei.items():
        val[1] = GetIp2ASFromBdrMapItDb(key)
    print('left neighbor: ')
    sort_list = sorted(left_nei.items(), key=lambda d:d[1][0], reverse=True)
    print(sort_list)
    print('right neighbor: ')
    sort_list = sorted(right_nei.items(), key=lambda d:d[1], reverse=True)
    print(sort_list)

def StatIncompleteTraces():
    os.chdir(global_var.par_path + global_var.traceroute_dir)
    for root,dirs,files in os.walk('.'):
        for filename in files:
            if filename.endswith('warts') or filename == 'process_trace.sh':
                continue
            total_num = 0
            aest_num = 0
            partial_num = 0
            with open(filename, 'r') as rf:
                curline = rf.readline()
                while curline:
                    if not curline.startswith('T'):
                        curline = rf.readline()
                        continue
                    total_num += 1
                    if curline.__contains__('q'):
                        aest_num += 1
                    elems = curline.strip('\n').split('\t')
                    dst_ip = elems[2]
                    if not elems[-1].__contains__(dst_ip):
                        partial_num += 1
                    curline = rf.readline()
                print('Total: %d, aest perc: %.2f, partial_perc: %.2f' %(total_num, aest_num * 100 / total_num, partial_num * 100 / total_num))

#vps = ['nrt-jp', 'per-au', 'syd-au', 'zrh2-ch']
def PlotAbStat_v2():
    #2016.4~2020.4, four years + 1 month, 49 months, thus each vp has 49 time-plot res
    res = dict()
    time_plots = 30
    os.chdir(global_var.par_path + 'tmp_out_my_anatrace/')
    for vp in global_var.vps:
        res[vp] = dict()#[dict() for i in range(0, 2)] #ab_0, ab_1
        for method in global_var.map_methods:
            #res[vp][method] = [0.0 for j in range(0, time_plots)]
            res[vp][method] = []
            for year in range(2018,2021):
                for month in range(1,13):
                    if (year == 2020 and month > 4): #2016年4月前peeringdb数据不准，2020年5月后的数据不全
                        continue
                    date = str(year) + str(month).zfill(2)
                    #offset = (year - 2018) * 12 + month - 1
                    grep_cmd = 'grep \'In total\' ' + vp + '_' + date + '15/' + method + '/' + 'record4_*'
                    output = os.popen(grep_cmd).read()
                    if output:
                        #res[vp][method][offset] = float(re.findall('ab precent: (.+?),', output)[0]) * 100
                        res[vp][method].append(float(re.findall('ab precent: (.+?),', output)[0]) * 100)
    
    #color_list = ['#000000', '#0000FF', '#8A2BE2', '#A52A2A', '#DEB887', '#5F9EA0', '#7FFF00', '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C', '#00FFFF', '#00008B', '#008B8B', '#B8860B', '#A9A9A9', '#006400', '#BDB76B', '#8B008B', '#556B2F', '#FF8C00', '#9932CC', '#8B0000', '#E9967A', '#8FBC8F', '#483D8B', '#2F4F4F', '#00CED1', '#9400D3', '#FF1493', '#00BFFF', '#696969', '#1E90FF', '#B22222', '#FFFAF0', '#228B22', '#FF00FF', '#DCDCDC', '#F8F8FF', '#FFD700', '#DAA520', '#808080', '#008000', '#ADFF2F', '#F0FFF0', '#FF69B4', '#CD5C5C', '#4B0082', '#FFFFF0', '#F0E68C', '#E6E6FA', '#FFF0F5', '#7CFC00', '#FFFACD', '#ADD8E6', '#F08080', '#E0FFFF', '#FAFAD2', '#90EE90', '#D3D3D3', '#FFB6C1', '#FFA07A', '#20B2AA', '#87CEFA', '#778899', '#B0C4DE', '#FFFFE0', '#00FF00', '#32CD32', '#FAF0E6', '#FF00FF', '#800000', '#66CDAA', '#0000CD', '#BA55D3', '#9370DB', '#3CB371', '#7B68EE', '#00FA9A', '#48D1CC', '#C71585', '#191970', '#F5FFFA', '#FFE4E1', '#FFE4B5', '#FFDEAD', '#000080', '#FDF5E6', '#808000', '#6B8E23', '#FFA500', '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98', '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F', '#FFC0CB', '#DDA0DD', '#B0E0E6', '#800080', '#FF0000', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#FAA460', '#2E8B57', '#FFF5EE', '#A0522D', '#C0C0C0', '#87CEEB', '#6A5ACD', '#708090', '#FFFAFA', '#00FF7F', '#4682B4', '#D2B48C', '#008080', '#D8BFD8', '#FF6347', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFFFFF', '#F5F5F5', '#FFFF00', '#9ACD32']
    color_list = ['red', 'yellow', 'blue']
    marker_list = ['o', '<', '*', 'H', '+']
    #fig = plt.figure()  )
    i = 0
    plt.figure()
    for vp in global_var.vps:
        j = 0        
        plt.subplot(511 + i)
        plt.ylim((0, 50))
        plt.xlim((0, 28))
        y_ticks = np.arange(0, 50, 10)
        # plt.ylim((0, 20))
        # y_ticks = np.arange(0, 20, 5)
        plt.tick_params(labelsize=6) 
        plt.yticks(y_ticks)
        frame = plt.gca()
        if i == 2:
            plt.ylabel('Mismatch Percent', fontsize=8)
        if i == 4:
            frame.axes.get_xaxis().set_visible(True)
            plt.xticks([0, 6, 12, 18, 24],[r'$2018.01$', r'$2018.06$', r'$2019.01$', r'$2019.06$', r'$2020.01$'])
            plt.xlabel('Date', fontsize=8)
        else:
            frame.axes.get_xaxis().set_visible(False)
        if i == 0:
            for method in global_var.map_methods:
                    #ax[i].scatter(range(0, time_plots), res[vp][method], c = color_list[j], marker= marker_list[ab_i], label = method + '_' + str(ab_i))
                plt.scatter(range(0, len(res[vp][method])), res[vp][method], c = color_list[j], marker= '*', s = 8., label='%s'%method)
                #plt.text(x, y , '%s' %vp, ha='center', va='center')
                j += 1
            plt.legend(loc='upper right', borderpad=0.5, labelspacing=1, prop={'size': 6}, ncol=3)
        else:
            for method in global_var.map_methods:
                    #ax[i].scatter(range(0, time_plots), res[vp][method], c = color_list[j], marker= marker_list[ab_i], label = method + '_' + str(ab_i))
                plt.scatter(range(0, len(res[vp][method])), res[vp][method], c = color_list[j], marker= marker_list[j], s = 8.)
                j += 1
        ax2 = plt.twinx()  # this is the important function
        #ax2.set_ylim((0, 50))
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        #ax2.set_visible(False)
        ax2.set_ylabel('            %s' %vp, rotation=0, fontsize=8)
        i += 1
    #plt.rcParams['figure.figsize'] = (1.0, 4.0)
    #plt(figsize=(8, 6), dpi=80)
    j = 0
    #'nrt-jp', 'per-au', 'syd-au', 'zrh2-ch', 'sjc2-us'
    plt.tight_layout()
    #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    eps_fig = plt.gcf() # 'get current figure'
    eps_fig.savefig('comparison.eps', format='eps')
    plt.show()

def PlotCdf():
    arr = []
    with open(global_var.par_path +  global_var.out_my_anatrace_dir + '/last_as_in_last_extra', 'r') as rf:
        curline = rf.readline()
        i = 0
        while curline:
            #arr.append(float(curline[curline.index(',')+1:-1]))
            count = int(curline.strip('\n').split(': ')[-1])
            for k in range(0, count):
                arr.append(i)
            i += 1
            curline = rf.readline()
            # if i > 100:
            #     break
    ecdf = sm.distributions.ECDF(arr)
    x = np.linspace(min(arr), max(arr), len(arr))
    y = ecdf(x)
    #print(y)
    y_dot5 = np.linspace(0.5, 0.5, len(arr))
    y_dot9 = np.linspace(0.9, 0.9, len(arr))
    idx_dot5 = np.argwhere(np.diff(np.sign(y - y_dot5))).flatten()
    idx_dot9 = np.argwhere(np.diff(np.sign(y - y_dot9))).flatten()
    plt.axhline(y=0.5, xmin = 0.0, xmax = x[idx_dot5[0]]/x[-1], color="red", linestyle="--")
    plt.axhline(y=0.9, xmin = 0.0, xmax = x[idx_dot9[0]]/x[-1], color="red", linestyle="--")
    plt.axvline(x=x[idx_dot5[0]], ymin = 0.0, ymax = 0.5, color="red", linestyle="--")
    plt.axvline(x=x[idx_dot9[0]], ymin = 0.0, ymax = 0.9, color="red", linestyle="--")
    plt.text(x[idx_dot5[0]]+5, 0.48, '(%d,0.5)'%x[idx_dot5[0]], color="red")
    plt.text(x[idx_dot9[0]]+5, 0.86, '(%d,0.9)'%x[idx_dot9[0]], color="red")
    plt.xlim((0, x[-1]))
    plt.ylim((0,1))
    #plt.axvline(x=22, ymin=0.0, ymax=0.33, color="red", linestyle="--")
    plt.plot(x, y, color='black')
    # cdf = stats.cumfreq(arr)
    # plt.plot(cdf[0])
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    eps_fig = plt.gcf() # 'get current figure'
    eps_fig.savefig('last_as_in_last_extra_cdf.eps', format='eps')
    plt.show()

def PlotCdf_2():
    color_list = ['yellow', 'black', 'blue', 'green', 'pink']
    max_x = 0
    j = 0
    for vp in global_var.vps:
        arr = []        
        with open(global_var.par_path + global_var.out_my_anatrace_dir + '/tmp_bifurc_' + vp, 'r') as rf:
            curline = rf.readline()
            curline = rf.readline()
            i = 0
            while curline:
                #arr.append(float(curline[curline.index(',')+1:-1]))
                count = int(curline[curline.index(':')+1:curline.index(',')])
                for k in range(0, count):
                    arr.append(i)
                i += 1
                curline = rf.readline()
                # if i > 100:
                #     break
        ecdf = sm.distributions.ECDF(arr)
        x = np.linspace(min(arr), max(arr), len(arr))
        y = ecdf(x)
        if x[-1] > max_x:
            max_x = x[-1]
        #print(y)
        # if vp == 'nrt-jp' or vp == 'sjc2-us':
        #     y_dot5 = np.linspace(0.5, 0.5, len(arr))
        #     idx_dot5 = np.argwhere(np.diff(np.sign(y - y_dot5))).flatten()
        #     plt.axhline(y=0.5, xmin = 0.0, xmax = x[idx_dot5[0]]/x[-1], color="red", linestyle="--")
        #     plt.axvline(x=x[idx_dot5[0]], ymin = 0.0, ymax = 0.5, color="red", linestyle="--")
        #     plt.text(x[idx_dot5[0]]+50, 0.48, '(%d,0.5)'%x[idx_dot5[0]], color="red")
        plt.plot(x, y, color=color_list[j], label='%s'%vp)
        j += 1
    plt.xlim((0,max_x))
    plt.ylim((0,1))
    plt.tick_params(labelsize=12)
    plt.legend(loc='lower right', borderpad=0.5, labelspacing=1, prop={'size': 12}, ncol=3)
    plt.tight_layout()
    eps_fig = plt.gcf() # 'get current figure'
    eps_fig.savefig('last_ab_cdf.eps', format='eps')
    plt.show()

if __name__ == '__main__':
    #PlotAbStat_v2()
    PlotCdf()
    #PlotCdf_2()
    #StatIncompleteTraces()
    #StatNobgp()
    #StatIxp()
    #StatAb()
    #StatMultipath()

    #PlotAbStat()

    #CheckUndo()

    #CollectMatchIpStat()
    #CalMatchIpStat()

    #StatClassify()

    #StatLastExtraHop()

    #StatAbPercent()

    # date = sys.argv[1] + '15'
    # vp = 'zrh2-ch'
    # ip = '212.36.135.22'
    # print(date)
    # os.chdir('/mountdisk1/ana_c_d_incongruity/out_my_anatrace/' + vp + '_' + date + '/ribs_midar_bdrmapit/')
    # os.system('cat 1_has_set 2_has_loop 3_single_path 4_multi_path 5_has_ixp_ip > test')
    # ConnectToBdrMapItDb(global_var.par_path + global_var.out_bdrmapit_dir + 'bdrmapit_' + vp + '_' + date + '.db')
    # ConstrBdrCache()
    # print(GetIp2ASFromBdrMapItDb(ip))
    # StatNeigborIp('test', ip)
    #StatNeigborIp('/mountdisk1/ana_c_d_incongruity/traceroute_data/zrh2-ch.' + date, '212.36.135.22')
