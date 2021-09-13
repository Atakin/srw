import math
import numpy as np
import time

##input_file_name=input("Enter the name of the input file:")
input_file_name = 'input.txt'
input_f = open(input_file_name)
inp_line = input_f.readline()
# выходные файлы
data_C_f = open('Data_C.txt', 'w')
data_H_f = open('Data_H.txt', 'w')
check_f = open('check.txt', 'w')
# Цикл по строкам входного файла
while inp_line:
    time1 = time.time()
    # Получение нужной информации из входного файла
    c = int(inp_line)  # number of C atoms
    h = int(input_f.readline())  # number of H atoms
    print(c)
    print(h)
    N = int(c * 15 + h * 2)  # number of bazis functions
    M = int(c * 3 + h / 2)  # number of MO
    t = int(math.ceil(M / 5))
    out_file_name = input_f.readline()
    input_file_name = input_f.readline()
    inp_line = input_f.readline()
    ##########
    # Запись названий в выходные файлы
    data_C_f.write('\n' + out_file_name + '\n')
    data_H_f.write('\n' + out_file_name + '\n')
    check_f.write(out_file_name + '\n')
    ##########
    # Получение информации из файла .out
    data_f = open('const.txt', 'w')
    out_f = open(out_file_name[:-1])
    line = out_f.readline()
    line_number = 1
    dots_array = []
    while line:
        line = out_f.readline()
        line_number += 1
        if line.find('EIGENVECTORS') != -1:
            line = out_f.readline()
            line_number += 1
            for k in range(0, t):
                for j in range(0, 4):
                    line = out_f.readline()
                    line_number += 1
                for j in range(0, N):
                    line = out_f.readline()
                    line_number += 1
                    line = line[14:]
                    line1 = line.strip()
                    data_f.write(line1 + '\n')

        if line.find('ELECTRON DENSITY INTEGRALS') != -1:
            dots_array.append(line_number)

    length = len(dots_array)

    data_f.close()
    out_f.close()

    if length > 3:
        out_f = open(out_file_name[:-1])
        dots_f = open('dots.txt', 'w')
        i = 1
        j = 1
        line = out_f.readline()
        while j < 4:
            i += 1
            line = out_f.readline()
            if (i == dots_array[j] - 2):
                j += 1
                dots_f.write(line)
        out_f.close()
        dots_f.close()
        dots = np.loadtxt('dots.txt')

    data = np.loadtxt('const.txt')
    const = np.zeros((N, M), float)
    for k in range(0, t):
        for i in range(0, 5):
            if ((k * 5 + i) == M):
                break
            for j in range(0, N):
                const[j][k * 5 + i] = data[k * N + j][i]

    ###########
    # Получение информации из файла .inp
    if (input_file_name[-1] == '\n'):
        input_file_name = input_file_name[:-1]
    inp_f = open(input_file_name)
    inp_f1 = open('inp.txt', 'w')
    line = inp_f.readline()
    i = 1
    while line:
        i += 1
        line = inp_f.readline()
        if line.find('$DATA') != -1:
            line = inp_f.readline()
            line = inp_f.readline()
            for i in range(0, c + h):
                line = inp_f.readline()
                line1 = line.strip()
                j = line1.find(' ')
                line2 = line1[j:]
                line3 = line2.strip()
                j = line3.find(' ')
                line4 = line3[j:]
                line5 = line4.strip()
                inp_f1.write(line5 + '\n')
    inp_f.close()
    inp_f1.close()
    inp = np.loadtxt('inp.txt')
    inp = inp * 1.889725989

    ################
    # Массив констант базиса 6-31G(сначала заполняется для C атомов потом для H
    fi = np.zeros((7, N * 2), float)
    # f1
    fi[0][0] = int(0)
    fi[1][0] = math.pow(2 * 3047.5248800 / math.pi, 3 / 4) * 0.001834737132
    fi[1][1] = 3047.5248800
    fi[2][0] = math.pow(2 * 457.3695180 / math.pi, 3 / 4) * 0.014037322813
    fi[2][1] = 457.3695180
    fi[3][0] = math.pow(2 * 103.9486850 / math.pi, 3 / 4) * 0.068842622264
    fi[3][1] = 103.9486850
    fi[4][0] = math.pow(2 * 29.2101553 / math.pi, 3 / 4) * 0.232184443216
    fi[4][1] = 29.2101553
    fi[5][0] = math.pow(2 * 9.2866630 / math.pi, 3 / 4) * 0.467941348435
    fi[5][1] = 9.2866630
    fi[6][0] = math.pow(2 * 3.1639270 / math.pi, 3 / 4) * 0.362311985337
    fi[6][1] = 3.1639270
    # f2
    fi[0][2] = int(0)
    fi[1][2] = -0.119332419775 * math.pow(2 * 7.8682723 / math.pi, 3 / 4)
    fi[1][3] = 7.8682723
    fi[2][2] = -math.pow(2 * 1.8812885 / math.pi, 3 / 4) * 0.160854151696
    fi[2][3] = 1.8812885
    fi[3][2] = math.pow(2 * 0.5442493 / math.pi, 3 / 4) * 1.143456437840
    fi[3][3] = 0.5442493

    # f3
    fi[0][4] = int(1)
    fi[1][4] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591)
    fi[1][5] = 7.8682723
    fi[2][4] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(1.8812885, 5 / 4) * 0.316423960957
    fi[2][5] = 1.8812885
    fi[3][4] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.5442493, 5 / 4) * 0.744308290898
    fi[3][5] = 0.5442493

    # f4
    fi[0][6] = int(2)
    fi[1][6] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591)
    fi[1][7] = 7.8682723
    fi[2][6] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(1.8812885, 5 / 4) * 0.316423960957
    fi[2][7] = 1.8812885
    fi[3][6] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.5442493, 5 / 4) * 0.744308290898
    fi[3][7] = 0.5442493

    # f5
    fi[0][8] = int(3)
    fi[1][8] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591)
    fi[1][9] = 7.8682723
    fi[2][8] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(1.8812885, 5 / 4) * 0.316423960957
    fi[2][9] = 1.8812885
    fi[3][8] = 2 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.5442493, 5 / 4) * 0.744308290898
    fi[3][9] = 0.5442493

    # f6
    fi[0][10] = int(0)
    fi[1][10] = math.pow(2 * 0.1687145 / math.pi, 3 / 4)
    fi[1][11] = 0.1687145

    # f7
    fi[0][12] = int(1)
    fi[1][12] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(0.1687145, 5 / 4))
    fi[1][13] = 0.1687145

    # f8
    fi[0][14] = int(2)
    fi[1][14] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(0.1687145, 5 / 4))
    fi[1][15] = 0.1687145

    # f9
    fi[0][16] = int(3)
    fi[1][16] = 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(0.1687145, 5 / 4))
    fi[1][17] = 0.1687145

    fi[0][18] = int(1)
    fi[0][19] = int(1)
    fi[1][18] = 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4)
    fi[1][19] = 0.8

    fi[0][20] = int(2)
    fi[0][21] = int(2)
    fi[1][20] = 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4)
    fi[1][21] = 0.8

    fi[0][22] = int(3)
    fi[0][23] = int(3)
    fi[1][22] = 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4)
    fi[1][23] = 0.8

    fi[0][24] = int(1)
    fi[0][25] = int(2)
    fi[1][24] = 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4)
    fi[1][25] = 0.8

    fi[0][26] = int(1)
    fi[0][27] = int(3)
    fi[1][26] = 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4)
    fi[1][27] = 0.8

    fi[0][28] = int(2)
    fi[0][29] = int(3)
    fi[1][28] = 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4)
    fi[1][29] = 0.8

    for i in range(30, c * 15 * 2):
        for j in range(0, 7):
            fi[j][i] = fi[j][i - 30]

    if (h > 0):
        # h1
        fi[0][c * 15 * 2] = int(0)
        fi[1][c * 15 * 2] = 0.033494604338 * math.pow(2 * 18.7311370 / math.pi, 3 / 4)
        fi[1][c * 15 * 2 + 1] = 18.7311370
        fi[2][c * 15 * 2] = math.pow(2 * 2.8253944 / math.pi, 3 / 4) * 0.234726953484
        fi[2][c * 15 * 2 + 1] = 2.8253944
        fi[3][c * 15 * 2] = math.pow(2 * 0.6401217 / math.pi, 3 / 4) * 0.813757326146
        fi[3][c * 15 * 2 + 1] = 0.6401217
        # h2
        fi[0][c * 15 * 2 + 2] = int(0)
        fi[1][c * 15 * 2 + 2] = math.pow(2 * 0.1612778 / math.pi, 3 / 4)
        fi[1][c * 15 * 2 + 3] = 0.1612778
        for i in range(c * 15 * 2 + 4, N * 2):
            for j in range(0, 7):
                fi[j][i] = fi[j][i - 4]


    # Массив функций в явном виде(ускоряет прогрмаму)
    def cfunc(x, y, z):
        """returns array of function values for carbon at a point (x,y,z)
        """

        def c1(x, y, z):
            return (math.pow(2 * 3047.5248800 / math.pi, 3 / 4) * 0.001834737132 * math.exp(
                -(x * x + y * y + z * z) * 3047.5248800) + math.pow(2 * 457.3695180 / math.pi,
                                                                    3 / 4) * 0.014037322813 * math.exp(
                -(x * x + y * y + z * z) * 457.3695180) + math.pow(2 * 103.9486850 / math.pi,
                                                                   3 / 4) * 0.068842622264 * math.exp(
                -(x * x + y * y + z * z) * 103.9486850) + math.pow(2 * 29.2101553 / math.pi,
                                                                   3 / 4) * 0.232184443216 * math.exp(
                -(x * x + y * y + z * z) * 29.2101553) + math.pow(2 * 9.2866630 / math.pi,
                                                                  3 / 4) * 0.467941348435 * math.exp(
                -(x * x + y * y + z * z) * 9.2866630) + math.pow(2 * 3.1639270 / math.pi,
                                                                 3 / 4) * 0.362311985337 * math.exp(
                -(x * x + y * y + z * z) * 3.1639270))

        def c2(x, y, z):
            return (-0.119332419775 * math.pow(2 * 7.8682723 / math.pi, 3 / 4) * math.exp(
                -(x * x + y * y + z * z) * 7.8682723) - math.pow(2 * 1.8812885 / math.pi,
                                                                 3 / 4) * 0.160854151696 * math.exp(
                -(x * x + y * y + z * z) * 1.8812885) + math.pow(2 * 0.5442493 / math.pi,
                                                                 3 / 4) * 1.143456437840 * math.exp(
                -(x * x + y * y + z * z) * 0.5442493))

        def c6(x, y, z):
            return math.pow(2 * 0.1687145 / math.pi, 3 / 4) * math.exp(-(x * x + y * y + z * z) * 0.1687145)

        def c3(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591 * x * math.exp(
                -(x * x + y * y + z * z) * 7.8682723) + math.pow(1.8812885, 5 / 4) * 0.316423960957 * x * math.exp(
                -(x * x + y * y + z * z) * 1.8812885) + math.pow(0.5442493, 5 / 4) * 0.744308290898 * x * math.exp(
                -(x * x + y * y + z * z) * 0.5442493))

        def c7(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (
                    math.pow(0.1687145, 5 / 4) * x * math.exp(-(x * x + y * y + z * z) * 0.1687145))

        def c4(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591 * y * math.exp(
                -(x * x + y * y + z * z) * 7.8682723) + math.pow(1.8812885, 5 / 4) * 0.316423960957 * y * math.exp(
                -(x * x + y * y + z * z) * 1.8812885) + math.pow(0.5442493, 5 / 4) * 0.744308290898 * y * math.exp(
                -(x * x + y * y + z * z) * 0.5442493))

        def c8(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (
                    math.pow(0.1687145, 5 / 4) * y * math.exp(-(x * x + y * y + z * z) * 0.1687145))

        def c5(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (math.pow(7.8682723, 5 / 4) * 0.068999066591 * z * math.exp(
                -(x * x + y * y + z * z) * 7.8682723) + math.pow(1.8812885, 5 / 4) * 0.316423960957 * z * math.exp(
                -(x * x + y * y + z * z) * 1.8812885) + math.pow(0.5442493, 5 / 4) * 0.744308290898 * z * math.exp(
                -(x * x + y * y + z * z) * 0.5442493))

        def c9(x, y, z):
            return 2 * math.pow(2 / math.pi, 3 / 4) * (
                    math.pow(0.1687145, 5 / 4) * z * math.exp(-(x * x + y * y + z * z) * 0.1687145))

        def c10(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4) * x * x * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        def c11(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4) * y * y * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        def c12(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) / math.sqrt(3) * math.pow(0.8, 7 / 4) * z * z * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        def c13(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4) * x * y * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        def c14(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4) * x * z * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        def c15(x, y, z):
            return 4 * math.pow(2 / math.pi, 3 / 4) * math.pow(0.8, 7 / 4) * y * z * math.exp(
                -(x * x + y * y + z * z) * 0.8)

        cfunc = np.zeros((15, 1), float)
        cfunc[0][0] = c1(x, y, z)
        cfunc[1][0] = c2(x, y, z)
        cfunc[2][0] = c3(x, y, z)
        cfunc[3][0] = c4(x, y, z)
        cfunc[4][0] = c5(x, y, z)
        cfunc[5][0] = c6(x, y, z)
        cfunc[6][0] = c7(x, y, z)
        cfunc[7][0] = c8(x, y, z)
        cfunc[8][0] = c9(x, y, z)
        cfunc[9][0] = c10(x, y, z)
        cfunc[10][0] = c11(x, y, z)
        cfunc[11][0] = c12(x, y, z)
        cfunc[12][0] = c13(x, y, z)
        cfunc[13][0] = c14(x, y, z)
        cfunc[14][0] = c15(x, y, z)
        return cfunc


    ##hydragen
    def hfunc(x, y, z):
        """returns array of function values for hydrogen at a point (x,y,z)
        """

        def h1(x, y, z):
            return (0.033494604338 * math.pow(2 * 18.7311370 / math.pi, 3 / 4) * math.exp(
                -(x * x + y * y + z * z) * 18.7311370) - math.pow(2 * 2.8253944 / math.pi,
                                                                  3 / 4) * 0.234726953484 * math.exp(
                -(x * x + y * y + z * z) * 2.8253944) + math.pow(2 * 0.6401217 / math.pi,
                                                                 3 / 4) * 0.813757326146 * math.exp(
                -(x * x + y * y + z * z) * 0.6401217))

        def h2(x, y, z):
            return math.pow(2 * 0.1612778 / math.pi, 3 / 4) * math.exp(-(x * x + y * y + z * z) * 0.1612778)

        hfunc = np.zeros((2, 1), float)
        hfunc[0][0] = h1(x, y, z)
        hfunc[1][0] = h2(x, y, z)
        return hfunc


    def func(x, y, z):
        """returns a column of function values for molecule at a point (x,y,z)"""
        mass = np.zeros((N, 1), float)
        for j in range(0, c):
            m = cfunc(x - inp[j][0], y - inp[j][1], z - inp[j][2])
            for i in range(0, 15):
                mass[j * 15 + i][0] = m[i][0]
        for j in range(c, c + h):
            if (h == 0):
                break
            m = hfunc(x - inp[j][0], y - inp[j][1], z - inp[j][2])
            for i in range(0, 2):
                mass[c * 15 + (j - c) * 2 + i][0] = m[i][0]
        return mass


    # Функция возвращает элеткронную плотность в точке(x,y,z)
    def n(x, y, z):
        """returns an electron density  at a point (x,y,z)
        """
        m = np.matmul(np.transpose(func(x, y, z)), const)
        return float(2 * np.linalg.norm(m) ** 2)


    # Фукция для численного интегрирования корень из плотность умножить на базисную функцию

    def integralA(N):
        """numerical integration"""
        q = 0
        A = np.zeros((N, 1), float)
        Np = 20  # чило точек на оси при интегрировании
        p1x = np.min(inp[:, 0])
        p1x -= 3
        p2x = np.max(inp[:, 0])
        p2x += 3
        ex = (p2x - p1x) / Np
        p1y = np.min(inp[:, 1])
        p1y -= 3
        p2y = np.max(inp[:, 1])
        p2y += 3
        ey = (p2y - p1y) / Np
        p1z = np.min(inp[:, 2])
        p1z -= 3
        p2z = np.max(inp[:, 2])
        p2z += 3
        ez = (p2z - p1z) / Np
        for x in np.arange(p1x, p2x, ex):
            q += 1
            print(q, 'from', Np)
            for y in np.arange(p1y, p2y, ey):
                for z in np.arange(p1z, p2z, ez):
                    f = func(x, y, z)
                    m = np.matmul(np.transpose(f), const)
                    n = 2 * np.linalg.norm(m) ** 2
                    for i in range(0, N):
                        A[i][0] += (math.sqrt(n) * f[i]) * ex * ey * ez
        return A


    # Фукция для точного интегрирования произведения базисных функций
    def integralB(i, k):
        """precise integration"""
        s = 0
        numb1 = int(i // 15)
        if (numb1 > (c - 1)):
            numb1 = int(c + ((i - c * 15) // 2))
        numb2 = int(k // 15)
        if (numb2 > (c - 1)):
            numb2 = int(c + (k - c * 15) // 2)
        for m in range(1, 7):
            for n in range(1, 7):
                if (fi[m][2 * i] * fi[n][2 * k] != 0):
                    a = fi[m][2 * i + 1] + fi[n][2 * k + 1]
                    b = np.zeros(3, float)
                    b[0] = fi[m][2 * i + 1] * inp[numb1][0] + fi[n][2 * k + 1] * inp[numb2][0]
                    b[1] = fi[m][2 * i + 1] * inp[numb1][1] + fi[n][2 * k + 1] * inp[numb2][1]
                    b[2] = fi[m][2 * i + 1] * inp[numb1][2] + fi[n][2 * k + 1] * inp[numb2][2]
                    cx = fi[m][2 * i + 1] * math.pow(inp[numb1][0], 2) + fi[n][2 * k + 1] * math.pow(inp[numb2][0], 2)
                    cy = fi[m][2 * i + 1] * math.pow(inp[numb1][1], 2) + fi[n][2 * k + 1] * math.pow(inp[numb2][1], 2)
                    cz = fi[m][2 * i + 1] * math.pow(inp[numb1][2], 2) + fi[n][2 * k + 1] * math.pow(inp[numb2][2], 2)
                    piexp = math.pow((math.pi / a), 3 / 2) * math.exp(-((a * cx - b[0] * b[0]) / a)) * math.exp(
                        -((a * cy - b[1] * b[1]) / a)) * math.exp(-((a * cz - b[2] * b[2]) / a))
                    pos1 = int(fi[0][2 * i] - 1)
                    pos2 = int(fi[0][2 * k] - 1)
                    pos12 = int(fi[0][2 * i + 1] - 1)
                    pos22 = int(fi[0][2 * k + 1] - 1)
                    if ((fi[0][2 * i + 1] < 0.5) and (fi[0][2 * k + 1] < 0.5)):
                        if ((fi[0][2 * i] < 0.5) and (fi[0][2 * k] < 0.5)):
                            s += fi[m][2 * i] * fi[n][2 * k] * piexp
                        elif ((fi[0][2 * i] < 0.5) and (fi[0][2 * k] > 0.5)):
                            s += fi[m][2 * i] * fi[n][2 * k] * piexp * (b[pos2] / a - inp[numb2][pos2])
                        elif ((fi[0][2 * i] > 0.5) and (fi[0][2 * k] < 0.5)):
                            s += fi[m][2 * i] * fi[n][2 * k] * piexp * (b[pos1] / a - inp[numb1][pos1])
                        elif (abs((fi[0][2 * i]) - (fi[0][2 * k])) > 0.5):
                            s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                    b[pos1] * b[pos2] / math.pow(a, 2) - inp[numb1][pos1] * b[pos2] / a -
                                    inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos1] * inp[numb2][pos2])
                        elif (abs((fi[0][2 * i]) - (fi[0][2 * k])) < 0.5):
                            s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                    (2 * b[pos1] * b[pos2] + a) / (2 * math.pow(a, 2)) - inp[numb1][pos1] * b[
                                pos2] / a - inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos1] * inp[numb2][pos2])
                        else:
                            print("ERROR")
                    elif ((fi[0][2 * i + 1] > 0.5) and (fi[0][2 * k + 1] < 0.5)):
                        if (abs(fi[0][2 * i + 1] - fi[0][2 * i]) > 0.5):
                            if (fi[0][2 * k] < 0.5):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        b[pos1] * b[pos12] / math.pow(a, 2) - inp[numb1][pos1] * b[pos12] / a -
                                        inp[numb1][pos12] * b[pos1] / a + inp[numb1][pos1] * inp[numb1][pos12])
                            elif fi[0][2 * k] > 0.5:
                                if ((abs(fi[0][2 * k] - fi[0][2 * i]) > 0.5) and (
                                        abs(fi[0][2 * k] - fi[0][2 * i + 1]) > 0.5)):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            b[pos1] * b[pos12] * b[pos2] / a ** 3 - inp[numb2][pos2] * b[pos1] * b[
                                        pos12] / a ** 2 - inp[numb1][pos12] * b[pos1] * b[pos2] / a ** 2 -
                                            inp[numb1][pos1] * b[pos12] * b[pos2] / a ** 2 + inp[numb1][pos12] *
                                            inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos1] * inp[numb2][pos2] *
                                            b[pos12] / a + inp[numb1][pos12] * inp[numb1][pos1] * b[pos2] / a -
                                            inp[numb1][pos12] * inp[numb1][pos1] * inp[numb2][pos2])
                                elif (abs(fi[0][2 * k] - fi[0][2 * i]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos1] * b[pos2] + a) * b[pos12] / (2 * math.pow(a, 3)) -
                                            inp[numb2][pos2] * b[pos1] * b[pos12] / a ** 2 - inp[numb1][pos1] * b[
                                                pos12] * b[pos2] / a ** 2 - inp[numb1][pos12] * (
                                                    2 * b[pos1] * b[pos2] + a) / (2 * math.pow(a, 2)) +
                                            inp[numb1][pos12] * inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos1] *
                                            inp[numb2][pos2] * b[pos12] / a + inp[numb1][pos12] * inp[numb1][pos1] *
                                            b[pos2] / a - inp[numb1][pos12] * inp[numb1][pos1] * inp[numb2][pos2])
                                elif (abs(fi[0][2 * k] - fi[0][2 * i + 1]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos12] * b[pos2] + a) * b[pos1] / (2 * math.pow(a, 3)) -
                                            inp[numb2][pos2] * b[pos1] * b[pos12] / a ** 2 - inp[numb1][pos12] * b[
                                                pos1] * b[pos2] / a ** 2 - inp[numb1][pos1] * (
                                                    2 * b[pos12] * b[pos2] + a) / (2 * math.pow(a, 2)) +
                                            inp[numb1][pos12] * inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos1] *
                                            inp[numb2][pos2] * b[pos12] / a + inp[numb1][pos12] * inp[numb1][pos1] *
                                            b[pos2] / a - inp[numb1][pos12] * inp[numb1][pos1] * inp[numb2][pos2])
                                else:
                                    print("ERROR")
                            else:
                                print("ERROR")
                        elif (abs(fi[0][2 * i + 1] - fi[0][2 * i]) < 0.5):
                            if (fi[0][2 * k] < 0.5):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos1] * b[pos12] + a) / (2 * math.pow(a, 2)) - 2 * b[pos1] *
                                        inp[numb1][pos12] / a + inp[numb1][pos1] * inp[numb1][pos12])
                            elif (fi[0][2 * k] > 0.5):
                                if (abs(fi[0][2 * k] - fi[0][2 * i]) > 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos1] * b[pos12] + a) * b[pos2] / (2 * math.pow(a, 3)) - 2 *
                                            inp[numb1][pos12] * b[pos1] * b[pos2] / a ** 2 - inp[numb2][pos2] * (
                                                    2 * b[pos1] * b[pos12] + a) / (2 * math.pow(a, 2)) + 2 *
                                            inp[numb1][pos12] * inp[numb2][pos2] * b[pos1] / a + inp[numb1][pos12] *
                                            inp[numb1][pos1] * b[pos2] / a - inp[numb1][pos12] * inp[numb1][pos1] *
                                            inp[numb2][pos2])
                                elif (abs(fi[0][2 * k] - fi[0][2 * i]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos1] ** 2 + 3 * a) * b[pos1] / (2 * a ** 3) - 2 * inp[numb1][
                                        pos1] * (2 * b[pos1] ** 2 + a) / (2 * a ** 2) - inp[numb2][pos2] * (
                                                    2 * b[pos1] ** 2 + a) / (2 * a ** 2) + 2 * inp[numb2][
                                                pos2] * inp[numb1][pos1] * b[pos1] / a + inp[numb1][pos1] ** 2 * b[
                                                pos1] / a - inp[numb1][pos1] ** 2 * inp[numb2][pos2])
                                else:
                                    print("ERROR")
                            else:
                                print("ERROR")
                        else:
                            print("ERROR")
                    elif ((fi[0][2 * i + 1] < 0.5) and (fi[0][2 * k + 1] > 0.5)):
                        if (abs(fi[0][2 * k + 1] - fi[0][2 * k]) > 0.5):
                            if (fi[0][2 * i] < 0.5):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        b[pos2] * b[pos22] / math.pow(a, 2) - inp[numb2][pos2] * b[pos22] / a -
                                        inp[numb2][pos22] * b[pos2] / a + inp[numb2][pos2] * inp[numb2][pos22])
                            elif (fi[0][2 * i] > 0.5):
                                if ((abs(fi[0][2 * i] - fi[0][2 * k]) > 0.5) and (
                                        abs(fi[0][2 * i] - fi[0][2 * k + 1]) > 0.5)):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            b[pos2] * b[pos22] * b[pos1] / a ** 3 - inp[numb1][pos1] * b[pos2] * b[
                                        pos22] / a ** 2 - inp[numb2][pos22] * b[pos2] * b[pos1] / a ** 2 -
                                            inp[numb2][pos2] * b[pos22] * b[pos1] / a ** 2 + inp[numb2][pos22] *
                                            inp[numb1][pos1] * b[pos2] / a + inp[numb2][pos2] * inp[numb1][pos1] *
                                            b[pos22] / a + inp[numb2][pos22] * inp[numb2][pos2] * b[pos1] / a -
                                            inp[numb2][pos22] * inp[numb2][pos2] * inp[numb1][pos1])
                                elif (abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos2] * b[pos1] + a) * b[pos22] / (2 * math.pow(a, 3)) -
                                            inp[numb1][pos1] * b[pos2] * b[pos22] / a ** 2 - inp[numb2][pos2] * b[
                                                pos22] * b[pos1] / a ** 2 - inp[numb2][pos22] * (
                                                    2 * b[pos2] * b[pos1] + a) / (2 * math.pow(a, 2)) +
                                            inp[numb2][pos22] * inp[numb2][pos2] * b[pos1] / a + inp[numb2][pos2] *
                                            inp[numb1][pos1] * b[pos22] / a + inp[numb2][pos22] * inp[numb2][pos2] *
                                            b[pos1] / a - inp[numb2][pos22] * inp[numb1][pos1] * inp[numb2][pos2])
                                elif (abs(fi[0][2 * i] - fi[0][2 * k + 1]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos22] * b[pos1] + a) * b[pos1] / (2 * math.pow(a, 3)) -
                                            inp[numb1][pos1] * b[pos2] * b[pos22] / a ** 2 - inp[numb2][pos22] * b[
                                                pos2] * b[pos1] / a ** 2 - inp[numb2][pos2] * (
                                                    2 * b[pos22] * b[pos1] + a) / (2 * math.pow(a, 2)) +
                                            inp[numb2][pos22] * inp[numb2][pos2] * b[pos1] / a + inp[numb2][pos2] *
                                            inp[numb1][pos1] * b[pos22] / a + inp[numb2][pos22] * inp[numb2][pos2] *
                                            b[pos1] / a - inp[numb2][pos22] * inp[numb1][pos1] * inp[numb2][pos2])
                                else:
                                    print("ERROR")
                            else:
                                print("ERROR")
                        elif (abs(fi[0][2 * k + 1] - fi[0][2 * k]) < 0.5):
                            if (fi[0][2 * i] < 0.5):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos2] * b[pos22] + a) / (2 * math.pow(a, 2)) - 2 * b[pos2] *
                                        inp[numb2][pos22] / a + inp[numb2][pos2] * inp[numb2][pos22])
                            elif (fi[0][2 * i] > 0.5):
                                if (abs(fi[0][2 * i] - fi[0][2 * k]) > 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos2] * b[pos22] + a) * b[pos1] / (2 * math.pow(a, 3)) - 2 *
                                            inp[numb2][pos22] * b[pos2] * b[pos1] / a ** 2 - inp[numb1][pos1] * (
                                                    2 * b[pos2] * b[pos22] + a) / (2 * math.pow(a, 2)) + 2 *
                                            inp[numb2][pos22] * inp[numb1][pos1] * b[pos2] / a + inp[numb2][pos22] *
                                            inp[numb2][pos2] * b[pos1] / a - inp[numb2][pos22] * inp[numb2][pos2] *
                                            inp[numb1][pos1])
                                elif (abs(fi[0][2 * k] - fi[0][2 * i]) < 0.5):
                                    s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                            (2 * b[pos2] ** 2 + 3 * a) * b[pos2] / (2 * a ** 3) - 2 * inp[numb2][
                                        pos2] * (2 * b[pos2] ** 2 + a) / (2 * a ** 2) - inp[numb1][pos1] * (
                                                    2 * b[pos2] ** 2 + a) / (2 * a ** 2) + 2 * inp[numb1][
                                                pos1] * inp[numb2][pos2] * b[pos2] / a + inp[numb2][pos2] ** 2 * b[
                                                pos2] / a - inp[numb2][pos2] ** 2 * inp[numb1][pos1])
                                else:
                                    print("ERROR")
                            else:
                                print("ERROR")
                        else:
                            print("ERROR")
                    elif ((fi[0][2 * i + 1] > 0.5) and (fi[0][2 * k + 1] > 0.5)):
                        if ((abs(fi[0][2 * i + 1] - fi[0][2 * i]) < 0.5) and (
                                abs(fi[0][2 * k + 1] - fi[0][2 * k]) < 0.5)):
                            if (abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5):

                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (4 * b[pos1] ** 4 + 12 * a * b[pos1] ** 2 + 3 * a ** 2) / (
                                        4 * a ** 4) - 2 * (inp[numb1][pos1] + inp[numb2][pos2]) * b[pos1] * (
                                                2 * b[pos1] ** 2 + 3 * a) / (2 * a ** 3) + (
                                                inp[numb1][pos1] ** 2 + 4 * inp[numb1][pos1] * inp[numb2][
                                            pos2] + inp[numb2][pos2] ** 2) * (2 * b[pos1] ** 2 + a) / (
                                                2 * a ** 2) - 2 * inp[numb1][pos1] * inp[numb2][pos2] * (
                                                inp[numb1][pos1] + inp[numb2][pos2]) * b[pos1] / a + inp[numb1][
                                            pos1] ** 2 * inp[numb2][pos2] ** 2)
                            elif (abs(fi[0][2 * i] - fi[0][2 * k]) > 0.5):

                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos1] + a) * (2 * b[pos2] + a) / (4 * a ** 4) - 2 * inp[numb2][
                                    pos2] * (2 * b[pos1] + a) * b[pos2] / (2 * a ** 3) - 2 * inp[numb1][pos1] * (
                                                2 * b[pos2] + a) * b[pos1] / (2 * a ** 3) + inp[numb2][
                                            pos2] ** 2 * (2 * b[pos1] + a) / (2 * a ** 2) + inp[numb1][
                                            pos1] ** 2 * (2 * b[pos2] + a) / (2 * a ** 2) + 4 * inp[numb1][pos1] *
                                        inp[numb2][pos2] * b[pos1] * b[pos2] / a ** 2 - 2 * inp[numb1][pos1] ** 2 *
                                        inp[numb2][pos2] * b[pos2] / a - 2 * inp[numb1][pos1] * inp[numb2][
                                            pos2] ** 2 * b[pos1] / a + inp[numb1][pos1] ** 2 * inp[numb2][
                                            pos2] ** 2)
                            else:
                                print("ERROR")
                        elif ((abs(fi[0][2 * i + 1] - fi[0][2 * i]) < 0.5) and (
                                abs(fi[0][2 * k + 1] - fi[0][2 * k]) > 0.5)):

                            if ((abs(fi[0][2 * i + 1] - fi[0][2 * k]) > 0.5) and (
                                    abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) > 0.5)):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos1] + a) * b[pos2] * b[pos22] / (2 * a ** 4) - 2 * inp[numb1][
                                    pos1] * b[pos1] * b[pos2] * b[pos22] / a ** 3 + inp[numb1][pos1] ** 2 * b[
                                            pos2] * b[pos22] / a ** 2 - inp[numb2][pos22] * b[pos2] * (
                                                2 * b[pos1] + a) / (2 * a ** 3) + 2 * inp[numb2][pos22] *
                                        inp[numb1][pos1] * b[pos2] * b[pos1] / a ** 2 - inp[numb1][pos1] ** 2 *
                                        inp[numb2][pos22] * b[pos2] / a - inp[numb2][pos2] * (2 * b[pos1] + a) * b[
                                            pos22] / (2 * a ** 3) + 2 * inp[numb2][pos2] * inp[numb1][pos1] * b[
                                            pos1] * b[pos22] / a ** 2 - inp[numb1][pos1] ** 2 * inp[numb2][pos2] *
                                        b[pos22] / a + inp[numb2][pos2] * inp[numb2][pos22] * (2 * b[pos1] + a) / (
                                                2 * a ** 2) - 2 * inp[numb1][pos1] * inp[numb2][pos2] *
                                        inp[numb2][pos22] * b[pos1] / a + inp[numb1][pos1] ** 2 * inp[numb2][pos2] *
                                        inp[numb2][pos22])
                            elif ((abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5) or (
                                    abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) < 0.5)):
                                if (abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5):
                                    same = int(fi[0][2 * k] - 1)
                                    dif = int(fi[0][2 * k + 1] - 1)
                                elif (abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) < 0.5):
                                    same = int(fi[0][2 * k + 1] - 1)
                                    dif = int(fi[0][2 * k] - 1)
                                else:
                                    print("ERROR")
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        b[same] * (2 * b[same] ** 2 + 3 * a) * b[dif] / (2 * a ** 4) - b[same] *
                                        inp[numb2][dif] * (2 * b[same] ** 2 + 3 * a) / (2 * a ** 3) - 2 *
                                        inp[numb1][same] * (2 * b[same] ** 2 + a) * b[dif] / (2 * a ** 3) -
                                        inp[numb2][same] * (2 * b[same] ** 2 + a) * b[dif] / (2 * a ** 3) + 2 *
                                        inp[numb1][same] * inp[numb2][dif] * (2 * b[same] ** 2 + a) / (2 * a ** 2) +
                                        inp[numb2][same] * inp[numb2][dif] * (2 * b[same] ** 2 + a) / (2 * a ** 2) +
                                        inp[numb1][same] ** 2 * b[same] * b[dif] / a ** 2 + 2 * inp[numb1][same] *
                                        inp[numb2][same] * b[same] * b[dif] / a ** 2 - inp[numb1][same] ** 2 *
                                        inp[numb2][dif] * b[same] / a - 2 * inp[numb1][same] * inp[numb2][dif] *
                                        inp[numb2][same] * b[same] / a - inp[numb1][same] ** 2 * inp[numb2][same] *
                                        b[dif] / a + inp[numb1][same] ** 2 * inp[numb2][same] * inp[numb2][dif])
                            else:
                                print("ERROR")
                        elif ((abs(fi[0][2 * i + 1] - fi[0][2 * i]) > 0.5) and (
                                abs(fi[0][2 * k + 1] - fi[0][2 * k]) < 0.5)):
                            if ((abs(fi[0][2 * i + 1] - fi[0][2 * k]) > 0.5) and (
                                    abs(fi[0][2 * i] - fi[0][2 * k]) > 0.5)):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos2] + a) * b[pos1] * b[pos12] / (2 * a ** 4) - 2 * inp[numb2][
                                    pos2] * b[pos2] * b[pos1] * b[pos12] / a ** 3 + inp[numb2][pos2] ** 2 * b[
                                            pos1] * b[pos12] / a ** 2 - inp[numb1][pos12] * b[pos1] * (
                                                2 * b[pos2] + a) / (2 * a ** 3) + 2 * inp[numb1][pos12] *
                                        inp[numb2][pos2] * b[pos1] * b[pos2] / a ** 2 - inp[numb2][pos2] ** 2 *
                                        inp[numb1][pos12] * b[pos1] / a - inp[numb1][pos1] * (2 * b[pos2] + a) * b[
                                            pos12] / (2 * a ** 3) + 2 * inp[numb1][pos1] * inp[numb2][pos2] * b[
                                            pos2] * b[pos12] / a ** 2 - inp[numb2][pos2] ** 2 * inp[numb1][pos1] *
                                        b[pos12] / a + inp[numb1][pos1] * inp[numb1][pos12] * (2 * b[pos2] + a) / (
                                                2 * a ** 2) - 2 * inp[numb2][pos2] * inp[numb1][pos1] *
                                        inp[numb1][pos12] * b[pos2] / a + inp[numb2][pos2] ** 2 * inp[numb1][pos1] *
                                        inp[numb1][pos12])
                            elif ((abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5) or (
                                    abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5)):
                                if abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5:
                                    same = int(fi[0][2 * i + 1] - 1)
                                    dif = int(fi[0][2 * i] - 1)
                                elif abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5:
                                    same = int(fi[0][2 * i] - 1)
                                    dif = int(fi[0][2 * i + 1] - 1)
                                else:
                                    print("ERROR")
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        b[same] * (2 * b[same] ** 2 + 3 * a) * b[dif] / (2 * a ** 4) - b[same] *
                                        inp[numb1][dif] * (2 * b[same] ** 2 + 3 * a) / (2 * a ** 3) - 2 *
                                        inp[numb2][same] * (2 * b[same] ** 2 + a) * b[dif] / (2 * a ** 3) -
                                        inp[numb1][same] * (2 * b[same] ** 2 + a) * b[dif] / (2 * a ** 3) + 2 *
                                        inp[numb2][same] * inp[numb1][dif] * (2 * b[same] ** 2 + a) / (2 * a ** 2) +
                                        inp[numb1][same] * inp[numb1][dif] * (2 * b[same] ** 2 + a) / (2 * a ** 2) +
                                        inp[numb2][same] ** 2 * b[same] * b[dif] / a ** 2 + 2 * inp[numb2][same] *
                                        inp[numb1][same] * b[same] * b[dif] / a ** 2 - inp[numb2][same] ** 2 *
                                        inp[numb1][dif] * b[same] / a - 2 * inp[numb2][same] * inp[numb1][dif] *
                                        inp[numb1][same] * b[same] / a - inp[numb2][same] ** 2 * inp[numb1][same] *
                                        b[dif] / a + inp[numb2][same] ** 2 * inp[numb1][same] * inp[numb1][dif])
                            else:
                                print("ERROR")
                        elif ((abs(fi[0][2 * i + 1] - fi[0][2 * i]) > 0.5) and (
                                abs(fi[0][2 * k + 1] - fi[0][2 * k]) > 0.5)):
                            if (((abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5) and (
                                    abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) < 0.5)) or (
                                    (abs(fi[0][2 * i] - fi[0][2 * k + 1]) < 0.5) and (
                                    abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5))):
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[pos1] ** 2 + a) * (2 * b[pos12] ** 2 + a) / (4 * a ** 4) - (
                                        2 * b[pos1] ** 2 + a) * b[pos12] / (2 * a ** 3) * (
                                                inp[numb2][pos12] + inp[numb1][pos12]) - (
                                                2 * b[pos12] ** 2 + a) * b[pos1] / (2 * a ** 3) * (
                                                inp[numb2][pos1] + inp[numb1][pos1]) + (
                                                2 * b[pos1] ** 2 + a) / (2 * a ** 2) * inp[numb1][pos12] *
                                        inp[numb2][pos12] + (2 * b[pos12] ** 2 + a) / (2 * a ** 2) * inp[numb1][
                                            pos1] * inp[numb2][pos1] + b[pos1] * b[pos2] / a ** 2 * (
                                                inp[numb1][pos1] * inp[numb1][pos12] + inp[numb1][pos1] *
                                                inp[numb2][pos12] + inp[numb2][pos1] * inp[numb1][pos12] +
                                                inp[numb2][pos1] * inp[numb2][pos12]) - b[pos1] / a * (
                                                inp[numb1][pos1] * inp[numb1][pos12] * inp[numb2][pos12] +
                                                inp[numb2][pos1] * inp[numb1][pos12] * inp[numb2][pos12]) - b[
                                            pos2] / a * (inp[numb1][pos1] * inp[numb2][pos1] * inp[numb1][pos12] +
                                                         inp[numb1][pos1] * inp[numb2][pos1] * inp[numb2][pos12]) +
                                        inp[numb1][pos1] * inp[numb2][pos1] * inp[numb1][pos12] * inp[numb2][pos12])
                            else:
                                if ((abs(fi[0][2 * i] - fi[0][2 * k]) < 0.5) and (
                                        abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) > 0.5)):
                                    same = int(fi[0][2 * i] - 1)
                                    dif1 = int(fi[0][2 * i + 1] - 1)
                                    dif2 = int(fi[0][2 * k + 1] - 1)
                                elif ((abs(fi[0][2 * i] - fi[0][2 * k + 1]) < 0.5) and (
                                        abs(fi[0][2 * i + 1] - fi[0][2 * k]) > 0.5)):
                                    same = int(fi[0][2 * i] - 1)
                                    dif1 = int(fi[0][2 * i + 1] - 1)
                                    dif2 = int(fi[0][2 * k] - 1)
                                elif ((abs(fi[0][2 * i + 1] - fi[0][2 * k]) < 0.5) and (
                                        abs(fi[0][2 * i] - fi[0][2 * k + 1]) > 0.5)):
                                    same = int(fi[0][2 * i + 1] - 1)
                                    dif1 = int(fi[0][2 * i] - 1)
                                    dif2 = int(fi[0][2 * k + 1] - 1)
                                elif ((abs(fi[0][2 * i + 1] - fi[0][2 * k + 1]) < 0.5) and (
                                        abs(fi[0][2 * i] - fi[0][2 * k]) > 0.5)):
                                    same = int(fi[0][2 * i + 1] - 1)
                                    dif1 = int(fi[0][2 * i] - 1)
                                    dif2 = int(fi[0][2 * k] - 1)
                                else:
                                    print("ERROR")
                                s += fi[m][2 * i] * fi[n][2 * k] * piexp * (
                                        (2 * b[same] ** 2 + a) * b[dif1] * b[dif2] / (2 * a ** 4) - (
                                        2 * b[same] ** 2 + a) * b[dif2] / (2 * a ** 3) * inp[numb1][dif1] - (
                                                2 * b[same] ** 2 + a) * b[dif1] / (2 * a ** 3) * inp[numb2][
                                            dif2] - b[dif1] * b[dif2] * b[same] / a ** 3 * (
                                                inp[numb1][same] + inp[numb2][same]) + b[same] * b[
                                            dif2] / a ** 2 * (
                                                inp[numb1][same] * inp[numb1][dif1] + inp[numb2][same] *
                                                inp[numb1][dif1]) + b[same] * b[dif1] / a ** 2 * (
                                                inp[numb1][same] * inp[numb2][dif2] + inp[numb2][same] *
                                                inp[numb2][dif2]) + (2 * b[same] ** 2 + a) / (2 * a ** 2) *
                                        inp[numb1][dif1] * inp[numb2][dif2] + b[dif2] * b[dif1] / a ** 2 *
                                        inp[numb1][same] * inp[numb2][same] - b[same] / a * (
                                                inp[numb1][same] * inp[numb2][dif2] * inp[numb1][dif1] +
                                                inp[numb2][same] * inp[numb2][dif2] * inp[numb1][dif1]) - b[
                                            dif2] / a * inp[numb1][same] * inp[numb2][same] * inp[numb1][dif1] - b[
                                            dif1] / a * inp[numb1][same] * inp[numb2][same] * inp[numb2][dif2] +
                                        inp[numb1][same] * inp[numb2][same] * inp[numb2][dif2] * inp[numb1][dif1])


                        else:
                            print("ERROR")

        return s


    ##    #А-столбец численных интегралов
    A = integralA(N)

    # В-матрица точных итегралов
    B = np.zeros((N, N), float)

    for i in range(0, N):
        for j in range(0, i + 1):
            B[i][j] = integralB(i, j)
    ##            print('B',i,j,B[i][j])
    for i in range(0, N - 1):
        for j in range(i + 1, N):
            B[i][j] = B[j][i]

    # C-стобец констант разложения в базис
    C = np.linalg.solve(B, A)

    # Заполнение выходных файлов
    inp01 = inp / 1.889725989
    for i in range(0, c):
        neib = 0
        lineC = 'C' + str(i + 1)
        for j in range(0, c):
            if j != i:
                if (math.sqrt((inp01[j][0] - inp01[i][0]) ** 2 + (inp01[j][1] - inp01[i][1]) ** 2 + (
                        inp01[j][2] - inp01[i][2]) ** 2) < 1.7):
                    lineC += '   C' + ' ' + str(inp01[i][0] - inp01[j][0]) + ' ' + str(
                        inp01[i][1] - inp01[j][1]) + ' ' + str(inp01[i][2] - inp01[j][2])
                    neib += 1
        for j in range(c, c + h):
            if (math.sqrt((inp01[j][0] - inp01[i][0]) ** 2 + (inp01[j][1] - inp01[i][1]) ** 2 + (
                    inp01[j][2] - inp01[i][2]) ** 2) < 1.2):
                lineC += '   H' + ' ' + str(inp01[i][0] - inp01[j][0]) + ' ' + str(
                    inp01[i][1] - inp01[j][1]) + ' ' + str(inp01[i][2] - inp01[j][2])
                neib += 1
        while neib < 4:
            lineC += '   X 0 0 0'
            neib += 1
        lineC += '\n'
        for k in range(i * 15, i * 15 + 15):
            lineC += str(C[k][0]) + ' '
        lineC += '\n'
        data_C_f.write(lineC)
    for i in range(c, c + h):
        lineH = 'H' + str(i - c + 1)
        for j in range(0, c):
            if (math.sqrt((inp01[j][0] - inp01[i][0]) ** 2 + (inp01[j][1] - inp01[i][1]) ** 2 + (
                    inp01[j][2] - inp01[i][2]) ** 2) < 1.2):
                lineH += '   C' + ' ' + str(inp01[i][0] - inp01[j][0]) + ' ' + str(
                    inp01[i][1] - inp01[j][1]) + ' ' + str(inp01[i][2] - inp01[j][2])
        lineH += '\n'
        for k in range(c * 15 + (i - c) * 2, c * 15 + (i - c) * 2 + 2):
            lineH += str(C[k][0]) + ' '
        lineH += '\n'
        data_H_f.write(lineH)

    if length > 3:
        ch_line = "Плотность в точках  "
        for i in range(0, 3):
            ch_line += '(' + str(dots[i][1]) + ',' + str(dots[i][2]) + ',' + str(dots[i][3]) + ') '
        ch_line += ', посчитанная GAMESS:\n'
        for i in range(0, 3):
            ch_line += str(dots[i][4]) + '  '
        check_f.write(ch_line)
        ch_line = "\nПлотность в точках  "
        for i in range(0, 3):
            ch_line += '(' + str(dots[i][1]) + ',' + str(dots[i][2]) + ',' + str(dots[i][3]) + ') '
        ch_line += ', посчитанная программой:\n'
        for i in range(0, 3):
            ch_line += str(n(dots[i][1], dots[i][2], dots[i][3])) + '  '
        check_f.write(ch_line)
        ch_line = "\nКорень из плотности в точках  "
        for i in range(0, 3):
            ch_line += '(' + str(dots[i][1]) + ',' + str(dots[i][2]) + ',' + str(dots[i][3]) + ') '
        ch_line += ', посчитанный программой:\n'
        for i in range(0, 3):
            ch_line += str(math.sqrt(n(dots[i][1], dots[i][2], dots[i][3]))) + '  '
        check_f.write(ch_line)
        ch_line = "\nКорень из плотности в точках  "
        for i in range(0, 3):
            ch_line += '(' + str(dots[i][1]) + ',' + str(dots[i][2]) + ',' + str(dots[i][3]) + ') '
        ch_line += ', посчитанный с помощью разложения в базис:\n'
        for i in range(0, 3):
            ch_line += str(float(np.matmul(np.transpose(func(dots[i][1], dots[i][2], dots[i][3])), C))) + '  '
        check_f.write(ch_line)
        check_f.write('\n\n')
    else:
        check_f.write('Недостаточно точек для сравнения с GAMESS\n')
        ch_line = "Плотность в точках (0,0,0) (0.5,0.5,0.5) (1,1,1):\n" + str(n(0, 0, 0)) + str(n(0.5, 0.5, 0.5)) + str(
            n(1, 1, 1)) + '\n'
        check_f.write(ch_line)
        ch_line = "Корень из плотности в точках (0,0,0) (0.5,0.5,0.5) (1,1,1):\n" + str(
            math.sqrt(n(0, 0, 0))) + '  ' + str(math.sqrt(n(0.5, 0.5, 0.5))) + '  ' + str(math.sqrt(n(1, 1, 1))) + '\n'
        check_f.write(ch_line)
        ch_line = "Корень плотности в точках (0,0,0) (0.5,0.5,0.5) (1,1,1) базиса:\n " + str(
            float(np.matmul(np.transpose(func(0, 0, 0)), C))) + '  ' + str(
            float(np.matmul(np.transpose(func(0.5, 0.5, 0.5)), C))) + '  ' + str(
            float(np.matmul(np.transpose(func(1, 1, 1)), C))) + '\n'
        check_f.write(ch_line)

    print("Плотность в точках (0,0,0) (0.5,0.5,0.5) (1,1,1):\n", n(0, 0, 0), n(0.5, 0.5, 0.5), n(1, 1, 1))
    print("Корень плотности в точках (0,0,0) (0.5,0.5,0.5) (1,1,1):\n", math.sqrt(n(0, 0, 0)),
          math.sqrt(n(0.5, 0.5, 0.5)), math.sqrt(n(1, 1, 1)))
    print("Корень плотности в точках (0,0,0) (0.5,0.5,0.5) (1,1,1) базиса:\n ",
          np.matmul(np.transpose(func(0, 0, 0)), C), np.matmul(np.transpose(func(0.5, 0.5, 0.5)), C),
          np.matmul(np.transpose(func(1, 1, 1)), C))
    time2 = time.time() - time1
    time2 = time2 / 60
    print("Время выполнения(мин):", time2)
print('end')
input_f.close()
data_C_f.close()
data_H_f.close()
check_f.close()
