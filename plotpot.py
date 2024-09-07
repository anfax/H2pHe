import matplotlib.pyplot as plt
import numpy as np 
import os 
import scienceplots
plt.style.use('science')

fil = 'pot_He_vs_H_2.dat'
if os.path.exists(fil):
    data = np.loadtxt(fil)
    R = data[:,0]
    V = data[:,1]
    
    fig, ax = plt.subplots()
    ax.plot(R, V, 'k-')
    ax.set_xlabel('$R$ (a.u.)')
    ax.set_ylabel('$V(R)$ (a.u.)')
    ax.set_title('Potential curve of He + $\mathrm{H_2}^+$')
    
    # 添加放大图
    axins = ax.inset_axes([0.3, 0.4, 0.6, 0.5])  # 放大图的位置和大小
    axins.plot(R, V, 'k-')
    axins.set_xlim(2.0, 8.0)
    axins.set_ylim(-0.02, 0.01)
    # axins.set_xticklabels([])
    # axins.set_yticklabels([])
    
    ax.indicate_inset_zoom(axins)
    
    plt.tight_layout()
    plt.savefig(fil[:-4]+'.jpeg', dpi=666)
    plt.show()
    plt.close()
else:
    print(fil + ' is not found!')