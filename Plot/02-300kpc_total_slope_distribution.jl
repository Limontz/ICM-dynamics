using CSV
using DataFrames
using Statistics
using PyCall
using PyPlot

main = "/home/marco/Scrivania/Ettori_project/"

# Specifica il percorso della tua cartella
folder_path_box_3D = string(main,"300_kpc_3D_box_filtering/")
folder_path_box_2D = string(main,"300_kpc_2D_box_filtering/")
folder_path_radial_3D = string(main,"02_3D_radial_filtering/")
folder_path_radial_2D = string(main,"02_2D_radial_filtering/")

# Crea un DataFrame vuoto
df3d_box = DataFrame()

# 100 kpc 3D box filtering reading
for file in readdir(folder_path_box_3D)
    if endswith(file, ".txt")
        df = CSV.File(joinpath(folder_path_box_3D, file), delim=',', header=true) |> DataFrame
        for col_name in names(df)[1:2]
            df3d_box[!, col_name] = df[!, col_name]
        end
    end
end

# 100 kpc 2D box filtering reading
df2d_box = DataFrame()

# Itera sui file e leggi i dati
for file in readdir(folder_path_box_2D)
    if endswith(file, ".txt")
        df = CSV.File(joinpath(folder_path_box_2D, file), delim=',', header=true) |> DataFrame
        for col_name in names(df)[1:2]
            new_col_name = Symbol(col_name * "_" * file[1])  # Aggiungi il nome del file per rendere il nome univoco
            df2d_box[!, new_col_name] = df[!, col_name]
        end
    end
end

# 01 shell 3D radial filtering reading
df3d_radial = DataFrame()

for file in readdir(folder_path_radial_3D)
    if endswith(file, ".txt")
        df = CSV.File(joinpath(folder_path_radial_3D, file), delim=',', header=true) |> DataFrame
        for col_name in names(df)[1:2]
            df3d_radial[!, col_name] = df[!, col_name]
        end
    end
end

# 100 kpc 2D box filtering reading
df2d_radial = DataFrame()

# Itera sui file e leggi i dati
for file in readdir(folder_path_radial_2D)
    if endswith(file, ".txt")
        df = CSV.File(joinpath(folder_path_radial_2D, file), delim=',', header=true) |> DataFrame
        for col_name in names(df)[1:2]
            new_col_name = Symbol(col_name * "_" * file[1])  # Aggiungi il nome del file per rendere il nome univoco
            df2d_radial[!, new_col_name] = df[!, col_name]
        end
    end
end

dtslope_mean_box = mean.(eachrow(df2d_box[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_mean_box = mean.(eachrow(df2d_box[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_mean_box = mean.(eachrow(df2d_box[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))

dtslope_max_box = maximum.(eachrow(df2d_box[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_max_box = maximum.(eachrow(df2d_box[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_max_box = maximum.(eachrow(df2d_box[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))

dtslope_min_box = minimum.(eachrow(df2d_box[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_min_box = minimum.(eachrow(df2d_box[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_min_box = minimum.(eachrow(df2d_box[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))



dtslope_mean_radial = mean.(eachrow(df2d_radial[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_mean_radial = mean.(eachrow(df2d_radial[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_mean_radial = mean.(eachrow(df2d_radial[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))

dtslope_max_radial = maximum.(eachrow(df2d_radial[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_max_radial = maximum.(eachrow(df2d_radial[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_max_radial = maximum.(eachrow(df2d_radial[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))

dtslope_min_radial = minimum.(eachrow(df2d_radial[:, [:dtslope_x, :dtslope_y, :dtslope_z]]))
dvslope_min_radial = minimum.(eachrow(df2d_radial[:, [:dvslope_x, :dvslope_y, :dvslope_z]]))
tvslope_min_radial = minimum.(eachrow(df2d_radial[:, [:tvslope_x, :tvslope_y, :tvslope_z]]))

fig, ax = subplots(figsize=(6,6))
fig.subplots_adjust(top=0.99,bottom=0.125,left=0.165, right=0.99)

upper_error_box = abs.(dtslope_max_box - dtslope_mean_box)
lower_error_box =  abs.(dtslope_min_box - dtslope_mean_box)
upper_error_radial = abs.(dtslope_max_radial - dtslope_mean_radial)
lower_error_radial =  abs.(dtslope_min_radial - dtslope_mean_radial)

errorbar(dtslope_mean_box, df3d_box.dtslope, xerr=(lower_error_box, upper_error_box), linestyle="", marker=".", markersize = 16, color="blue", label="box filtering", markeredgecolor="black")
errorbar(dtslope_mean_radial, df3d_radial.dtslope, xerr=(lower_error_radial, upper_error_radial), linestyle="", marker=".", markersize = 16, color="red", label="radial_filtering", markeredgecolor="black")
xlabel(L"(\delta T/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta T/\delta \rho)_{3D}",fontsize=20)
xticks([-0.4, -0.2, 0, 0.2],fontsize=18)
yticks([-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2],fontsize=18)
xlim(-0.4, 0.3)
ylim(-0.4, 0.2)
legend(fontsize=15, loc="upper right")
show()
readline()
clf()

upper_error_box = abs.(dvslope_max_box - dvslope_mean_box)
lower_error_box =  abs.(dvslope_min_box - dvslope_mean_box)
upper_error_radial = abs.(dvslope_max_radial - dvslope_mean_radial)
lower_error_radial =  abs.(dvslope_min_radial - dvslope_mean_radial)

errorbar(dvslope_mean_box, df3d_box.dvslope, xerr=(lower_error_box, upper_error_box), linestyle="", marker="s", markersize = 7, color="blue", label="box filtering", markeredgecolor="black")
errorbar(dvslope_mean_radial, df3d_radial.dvslope, xerr=(lower_error_radial, upper_error_radial), linestyle="", marker="s", markersize = 7, color="red", label="radial_filtering", markeredgecolor="black")
xlabel(L"(\delta v/\delta \rho)_{2D}",fontsize=20)
ylabel(L"(\delta v/\delta \rho)_{3D}",fontsize=20)
xticks(fontsize=18)
yticks(fontsize=18)
#xlim(-0.4, 0.3)
#ylim(-0.4, 0.2)
legend(fontsize=15, loc="lower right")
show()
readline()
clf()



upper_error_box = abs.(tvslope_max_box - tvslope_mean_box)
lower_error_box =  abs.(tvslope_min_box - tvslope_mean_box)
upper_error_radial = abs.(tvslope_max_radial - tvslope_mean_radial)
lower_error_radial =  abs.(tvslope_min_radial - tvslope_mean_radial)

errorbar(tvslope_mean_box, df3d_box.tvslope, xerr=(lower_error_box, upper_error_box), linestyle="", marker="D", markersize = 7, color="blue", label="box filtering", markeredgecolor="black")
errorbar(tvslope_mean_radial, df3d_radial.tvslope, xerr=(lower_error_radial, upper_error_radial), linestyle="", marker="D", markersize = 7, color="red", label="radial_filtering", markeredgecolor="black")
xlabel(L"(\delta v/\delta T)_{2D}",fontsize=20)
ylabel(L"(\delta v/\delta T)_{3D}",fontsize=20)
xticks([0.6, 1.0, 1.4, 1.8, 2.2],fontsize=18)
yticks([1.4, 1.6, 1.8, 2.0, 2.2],fontsize=18)
#xlim(-0.4, 0.3)
#ylim(-0.4, 0.2)
legend(fontsize=15, loc="upper left")
show()
