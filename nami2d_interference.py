# wave diffraction and interference, free end refrections
# 波の回折と干渉 自由端反射 2次元版

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# param memo x52,y52 137で壁に届く
# param memo x50,y100 170で角に届く
# param memo x150,y300 step20 510で壁に届く

X_SU = 150+2
Y_SU = 300+1
T_STEP = 20
SLIT_POS = 10
WALL_POS = 10
SRC_POW = 0.2

fig = plt.figure(figsize = (12.8, 7.2))
ax = fig.add_subplot(111, projection="3d")

plt.subplots_adjust(left=0, bottom=0, right=1, top=1)

x, y = np.meshgrid(range(0,X_SU),range(0,Y_SU))
z = 0.0 * x * y
edge = np.full((Y_SU,X_SU),False)

wall = int(X_SU/WALL_POS)
slit1 = int(Y_SU/2+Y_SU/SLIT_POS+0.5)-1
slit2 = int(Y_SU/2-Y_SU/SLIT_POS+0.5)-1

#edge[0,:] = True
#edge[-1,:] = True
edge[:,0] = True
edge[:,-1] = True
edge[:,wall] = True
edge[slit1,wall] = False
edge[slit2,wall] = False

for i in range(Y_SU):
    for j in range(X_SU):
        if edge[i][j]:
            print('o',end='')
        else:
            print('.',end='')
    print(':')

z1 = z.copy()
k = z.copy()

def update(i):
    global z,z1,k
    print(i) #,np.max(z[:,wall+1:]))

    # 上下の壁を吸収端にする
    z[0,:] = z[1,:] - k[1,:] * 2.32475437
    z[-1,:] = z[-2,:] - k[-2,:] * 2.32475437

    z[int(Y_SU/2+0.5)-1,:wall] = SRC_POW*np.sin(2*np.pi*i/T_STEP)

    for ix in range(1,X_SU-1):
        for iy in range(1,Y_SU-1):
            if not edge[iy][ix]:
                z1[iy,ix] = np.sum(z[iy-1:iy+2,ix]) + z[iy,ix-1] + z[iy,ix+1]
                if edge[iy-1,ix]:
                    z1[iy,ix] += z[iy,ix]
                if edge[iy+1,ix]:
                    z1[iy,ix] += z[iy,ix]
                if edge[iy,ix-1]:
                    z1[iy,ix] += z[iy,ix]
                if edge[iy,ix+1]:
                    z1[iy,ix] += z[iy,ix]
                z1[iy,ix] = z1[iy,ix]/5 + k[iy,ix]
    k = (z1 - z)
    ax.clear()
    vminmax = np.max(z[:,wall+int((X_SU-wall)/5):])
    vminmax = np.max([vminmax,0.01])
    pl = ax.plot_surface(x[1:-1,wall+1:-1], y[1:-1,wall+1:-1], z[1:-1,wall+1:-1], cmap = "viridis", vmin=-vminmax, vmax=vminmax, ccount=X_SU-2, rcount=Y_SU-1)
    pl = ax.plot_surface(x[1:-1,wall:wall+2], y[1:-1,wall:wall+2], z[1:-1,wall:wall+2], cmap = "hot", vmin=-vminmax, vmax=vminmax, ccount=X_SU-2, rcount=Y_SU-1)
    vminmax = np.max(z)
    vminmax = np.max([vminmax,SRC_POW])
    pl = ax.plot_surface(x[1:-1,1:wall+1], y[1:-1,1:wall+1], z[1:-1,1:wall+1], cmap = "viridis", vmin=-vminmax, vmax=vminmax, ccount=X_SU-2, rcount=Y_SU-1)
    ax.set_zlim(-1.0,1.0)
    XY_SU=max([X_SU,Y_SU])
    XYHAMI = int(XY_SU*0.2)
    ax.set_xlim(XYHAMI-XY_SU/4,XY_SU*3/4-XYHAMI)
    ax.set_ylim(XYHAMI,XY_SU-XYHAMI)

    z = z1.copy()

ani = animation.FuncAnimation(fig, update, interval=120, blit=False, save_count=500)
#ani.save("nami2d_interference.mp4")
plt.show()
