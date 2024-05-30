import sys, os
import numpy as np
import pandas as pd
import seaborn as sns


def process_camera(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # sort by adjusted p-value
            high_in_ko.append([line[0], float(line[4])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[4])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_roast(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # sort by FDR.Mixed
            high_in_ko.append([line[0], float(line[8])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[8])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_safe(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[4])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[4])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_gage(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[4])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[4])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_padog(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[8])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[8])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_plage(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[5])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[5])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_zscore(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[5])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[5])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_gsva(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[5])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[5])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_ssgsea(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[5])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[5])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_ora(file_name, result_dir):
    target = file_name.split('.')[0]
    file = open(result_dir + file_name, 'r')
    info = file.readlines()
    tmp_list = []
    high_in_ko = []
    high_in_control = []
    for line in info[1:]:
        line = line.strip().split('\t')
        if line[-2] == '1':
            # adjusted p-value
            high_in_ko.append([line[0], float(line[2])])
        elif line[-2] == '-1':
            high_in_control.append([line[0], float(line[2])])
    # tmp_list.sort(key = lambda x: x[1])
    # the smaller adjusted p-value, the better
    high_in_control.sort(key=lambda x: x[1])
    # wrong direction of enrichment, the less significant the adj p-value the better
    high_in_ko.sort(key=lambda x: x[1])
    high_in_ko.sort()

    tmp_list.extend(high_in_control)
    tmp_list.extend(high_in_ko)

    # tmp_list.reverse()
    rank = 1
    flag = False
    for t in tmp_list:
        if t[0] == target:
            flag = True
            return rank
        else:
            rank += 1
    file.close()
    return 0


def process_egsea_result(result_dir):
    all_result_files = os.listdir(result_dir)
    all_result_files.sort()
    rank_dict = dict()
    for file in all_result_files:
        if not file.endswith('txt'):
            continue
        # egsea_method = file.split('.')[-2]
        target = file.split('.')[0]
        # print(file, result_dir)
        rank = 0
        if file.__contains__('camera'):
            rank = process_camera(file, result_dir)
        elif file.__contains__('roast'):
            rank = process_roast(file, result_dir)
        elif file.__contains__('safe'):
            rank = process_safe(file, result_dir)
        elif file.__contains__('gage'):
            rank = process_gage(file, result_dir)
        elif file.__contains__('padog'):
            rank = process_padog(file, result_dir)
        elif file.__contains__('plage'):
            rank = process_plage(file, result_dir)
        elif file.__contains__('zscore'):
            rank = process_zscore(file, result_dir)
        elif file.__contains__('gsva'):
            rank = process_gsva(file, result_dir)
        elif file.__contains__('ssgsea'):
            rank = process_ssgsea(file, result_dir)
        elif file.__contains__('ora'):
            rank = process_ora(file, result_dir)

        if rank != 0:
            rank_dict[target] = rank
    return rank_dict