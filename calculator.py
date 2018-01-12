import math #数学运算相关库
import os


def least_square_cal(x_array,y_array,z_array):
    '''
    【最小二乘法计算】
    输入：x坐标数组，y坐标数组，z坐标数组
    普通平面计算拟合
    输出：
        拟合平面得出的系数a，b，c，平面方程为ax+by+c=0
    '''
    x_x = 0; y_y = 0; z_z = 0; x_y = 0; x_z = 0; y_z = 0
    sum_x = 0; sum_y = 0; sum_z = 0
    for index in range(len(x_array)):
        x_x += math.pow(x_array[index],2)
        y_y += math.pow(y_array[index],2)
        z_z += math.pow(z_array[index],2)
        x_y += x_array[index]*y_array[index]
        x_z += x_array[index]*z_array[index]
        y_z += y_array[index]*z_array[index]

        sum_x += x_array[index]
        sum_y += y_array[index]
        sum_z += z_array[index]

    s11 = x_x - math.pow(sum_x,2)/len(x_array)
    s12 = x_y - sum_x*sum_y/len(x_array)
    s13 = x_z - sum_x*sum_z/len(x_array)
    s23 = y_z - sum_y * sum_z /len(x_array)

    factor_a = (s12*s23-s13*s11)/(math.pow(s12,2)-math.pow(s11,2))
    factor_b = (s12*s13-s11*s23)/(math.pow(s12,2)-math.pow(s11,2))
    factor_c = (sum_z-factor_a*sum_x-factor_b*sum_y)/len(x_array)
    return factor_a, factor_b, factor_c

def least_square_proj_cal(x_array,y_array,z_array,factor_a, factor_b, factor_c):
    '''
    【最小二乘法投影计算】
    输入：x坐标数组，y坐标数组，z坐标数组, factor_a,b,c: 参考平面ax+by+c的参数
    投影平面计算拟合
    输出：
        拟合空间平面得出的系数a, b, c, d，平面方程为ax+by+cz+d=0
    '''
    x_x = 0; y_y = 0; z_z = 0; x_y = 0; x_z = 0; y_z = 0
    sum_x = 0; sum_y = 0; sum_z = 0
    for index in range(len(x_array)):
        x_x += math.pow(x_array[index],2)
        y_y += math.pow(y_array[index],2)
        z_z += math.pow(z_array[index],2)
        x_y += x_array[index]*y_array[index]
        x_z += x_array[index]*z_array[index]
        y_z += y_array[index]*z_array[index]

        sum_x += x_array[index]
        sum_y += y_array[index]
        sum_z += z_array[index]

    slope_1 = (len(x_array)*x_z-sum_x*sum_z)/(len(x_array)*z_z-math.pow(sum_z,2))
    slope_2 = (len(x_array)*y_z-sum_y*sum_z)/(len(x_array)*z_z-math.pow(sum_z,2))
    intercept_1 = (sum_x*z_z-x_z*sum_z)/(len(x_array)*z_z-math.pow(sum_z,2))
    intercept_2 = (sum_y*z_z-y_z*sum_z)/(len(x_array)*z_z-math.pow(sum_z,2))
        
    factor_aa = slope_2 + factor_b
    factor_bb = -(factor_a+slope_1)
    factor_cc = slope_2*(-factor_bb)-slope_1*factor_aa
    factor_dd = (-factor_bb)*intercept_2 - factor_aa*intercept_1

    return factor_aa, factor_bb, factor_cc, factor_dd

def file_parser(file_name):
    '''
    【文件处理器】
    因为所有数据文件格式相同，将其格式化
    输入：文件名
    输出：文件中数据点的三个坐标数组
    '''
    try:
        data_file = open(file_name, encoding='utf-8')
    except:
        print("open file error")
    data = []
    x_array = []
    y_array = []
    z_array = []
    line_index = 0
    for line in data_file.readlines():
        line_index+=1
        line_strip = line.strip()
        data_line = line_strip.split(' ')
        try:
            assert len(data_line)==3
        except AssertionError:
            print("this line has uncorrect data:",line_index)
        x_array.append(float(data_line[0]))
        y_array.append(float(data_line[1]))
        z_array.append(float(data_line[2]))
    
    return x_array,y_array,z_array


def flatness_cal(flat_file_name):
    '''
    【平面度计算】
    1.读取文件获取单个数据点坐标值，存入数组
    2.调用最小二乘法函数计算平面参数
    3.计算单个点到平面的距离
    4.计算平面度误差
    '''
    x_array, y_array, z_array = file_parser(flat_file_name)

    factor_a, factor_b, factor_c = least_square_cal(x_array, y_array, z_array)

    average = math.sqrt(factor_a*factor_a+factor_b*factor_b+1)
    distance = []

    for index in range(len(x_array)):
        single_dist = (factor_a*x_array[index]+factor_b*y_array[index]-z_array[index]+factor_c)/average
        distance.append(single_dist)
    
    max_distance = max(distance)
    min_distance = min(distance)

    flatness = max_distance - min_distance

    print("平面度：",flatness)

    return flatness

def verticality_cal(vert_file_names):
    '''
    【垂直度计算】
    输入：两个数据文件名，其中第一个为基准平面，第二个为计算对象
    输出：垂直度
    计算：
    1.读取基准平面数据，最小二乘法拟合平面
    2.读取待计算平面数据，计算投影点坐标
    3.计算垂直平面方程
    4.计算平面2数据点到拟合平面的距离，得出垂直度
    '''
    try:
        assert type(vert_file_names) is list
        assert len(vert_file_names)==2
    except AssertionError:
        print("需要两个平面的数据文件")
    x_std_array, y_std_array, z_std_array = file_parser(vert_file_names[0])
    factor_a, factor_b, factor_c = least_square_cal(x_std_array, y_std_array, z_std_array)

    x_array, y_array, z_array = file_parser(vert_file_names[1])
    x_project = []
    y_project = []
    z_project = []
    for index in range(len(x_array)):
        single_x = (((factor_b*factor_b+1)/factor_a)*x_array[index]-factor_b*y_array[index]+z_array[index]-factor_c)/(factor_a+factor_b*factor_b/factor_a+1/factor_a)
        single_y = (factor_b/factor_a)*(single_x-x_array[index])+y_array[index]
        single_z = factor_a*single_x+factor_b*single_y+factor_c

        x_project.append(single_x)
        y_project.append(single_y)
        z_project.append(single_z)

    factor_aa, factor_bb, factor_cc, factor_dd = least_square_proj_cal(x_project, y_project, z_project,factor_a,factor_b,factor_c)

    average = math.sqrt(factor_aa*factor_aa+factor_bb*factor_bb+factor_cc*factor_cc)

    distance = []

    for index in range(len(x_array)):
        single_dist = (factor_aa*x_array[index]+factor_bb*y_array[index]+factor_cc*z_array[index]+factor_dd)/average
        distance.append(single_dist)

    max_distance = max(distance)
    min_distance = min(distance)

    verticality = max_distance-min_distance
    print("垂直度：",verticality)

    return verticality

def cylindricity_cal():
    '''
    【圆柱度计算】
    TODO
    '''
    return

def gradient_cal():
    '''
    【倾斜度计算】
    TODO
    '''
    return

def parallelism_cal():
    '''
    【平行度计算】
    TODO
    '''
    return

def main(file_name, cal_type):
    try:
        assert type(file_name) is list
    except AssertionError:
        print ("file names should be list")
    try:
        assert cal_type>=0 and cal_type<=4
    except AssertionError:
        print ("cal type should be an integer range from 0 to 4")
    
    if cal_type == 0:
        flatness = flatness_cal(file_name[0])
    elif cal_type == 1:
        verticality = verticality_cal(file_name)
    else:
        print ("功能暂未开通")
    return

if __name__ == "__main__":
    file_name = ["文件名","文件名"] #定义输入文件名字，因为有打开多个文件的需求，所以为数组
    cal_type = 0 #定义计算类型，0-平面度，1-垂直度，2-圆柱度，3-倾斜度，4-平行度
    main(file_name, cal_type) #执行主函数