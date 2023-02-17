## Integration
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import griddata


Name_file_serie = ["study7"]

linetype = ["r--", "b--", "y--", "g--", "m--"]

for idx, job in enumerate(Name_file_serie):
    legend = [
        "$n=0$",
        "$n=2$",
        "$n=4$",
        "$n=6$",
        "$\\theta=45\degree$",
        "$\\theta=60\degree$",
    ]
    name_file = job
    Int_name = name_file + "Int.npz"
    Dis_name = name_file + "Dis.npz"
    For_name = name_file + "For.npz"
    # phi_name = name_file+'phi.npz'
    data = np.load(Int_name)["arr_0"]
    New_data = []
    Dis_data = np.load(Dis_name)["arr_0"] * 1000
    For_data = np.load(For_name)["arr_0"]
    # phi_data=np.load(phi_name)['arr_0']
    for steps in data:
        # sort y data
        step_info = np.array(steps)
        step_info = step_info[np.lexsort((step_info[:, 2],))]
        New_data.append(step_info)
    Integration_val = []
    # Obtain  current vs time_steps
    # Obtain  Force vs displacement

    for steps in New_data:

        X_coord = steps[:, 1]
        Y_coord = steps[:, 2]
        concateXY = np.concatenate((X_coord, Y_coord), axis=0)
        points = np.transpose(concateXY.reshape(2, int(len(concateXY) / 2)))
        values = steps[:, 3]
        grid_z1 = lambda x, y: griddata(points, values, (x, y), method="linear")
        x = np.linspace(np.min(X_coord), np.max(X_coord), 100)
        y = np.linspace(np.min(Y_coord), np.max(Y_coord), 100)
        z = grid_z1(x[:, None], y)

        int_surf = -integrate.trapz(integrate.trapz(z, y), x)

        print(int_surf)
        Integration_val.append(int_surf)
    # plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(
        Dis_data, 1000 * np.array(Integration_val), linetype[idx], label=legend[idx]
    )
    plt.xlabel("Displacement [mm]")
    plt.ylabel("Current [mA]")
    plt.ylim([0.44, 0.55])
    plt.xlim([0.3, 0.39])
    # plt.figure(2)
    # plt.plot(Dis_data,For_data)
    # plt.xlabel('Displacement [mm]')
    # plt.ylabel('Force [N]')
    plt.legend()




for idx, job in enumerate(Name_file_serie):
    legend = [
        "$n=0$",
        "$n=2$",
        "$n=4$",
        "$n=6$",
        "$\\theta=45\degree$",
        "$\\theta=60\degree$",
    ]
    name_file = job
    Int_name = name_file + "Int.npz"
    Dis_name = name_file + "Dis.npz"
    For_name = name_file + "For.npz"
    data = np.load(Int_name)["arr_0"]
    New_data = []
    Dis_data = np.load(Dis_name)["arr_0"] * 1000
    For_data = np.load(For_name)["arr_0"]
    for steps in data:
        # sort y data
        step_info = np.array(steps)
        step_info = step_info[np.lexsort((step_info[:, 2],))]
        New_data.append(step_info)
    Integration_val = []
    # Obtain  current vs time_steps
    # Obtain  Force vs displacement

    for steps in New_data:

        X_coord = steps[:, 1]
        Y_coord = steps[:, 2]
        concateXY = np.concatenate((X_coord, Y_coord), axis=0)
        points = np.transpose(concateXY.reshape(2, int(len(concateXY) / 2)))
        values = steps[:, 3]
        grid_z1 = lambda x, y: griddata(points, values, (x, y), method="linear")
        x = np.linspace(np.min(X_coord), np.max(X_coord), 100)
        y = np.linspace(np.min(Y_coord), np.max(Y_coord), 100)
        z = grid_z1(x[:, None], y)

        int_surf = -integrate.trapz(integrate.trapz(z, y), x)
        Integration_val.append(int_surf)
    # plt.figure(1)
    plt.subplot(2, 1, 2)
    # plt.plot(Dis_data,np.array(Integration_val),label=legend[idx])
    # plt.xlabel('Displacement [mm]')
    # plt.ylabel('Current [A]')
    # plt.ylim([0,0.016])
    # plt.figure(2)
    plt.ylim([0, 7.5])
    plt.xlim([0, 0.399])
    plt.plot(Dis_data, For_data, linetype[idx], label=legend[idx])
    plt.xlabel("Displacement [mm]")
    plt.ylabel("Force [N]")
    plt.legend()


plt.tight_layout()
plt.savefig("plot_figure.png", dpi=600)
