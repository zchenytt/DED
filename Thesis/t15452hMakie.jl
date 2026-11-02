# phase1
rgap1 = [1.499, 1.29, 2.331, 2.52, 0.617, 0.542, 1.315, 2.58, 1.096, 0.98, 0.995, 1.725, 1.344, 1.251, 1.239, 2.516, 1.264, 0.46, 1.373, 1.79, 1.184, 2.057, 0.604, 0.66, 1.402, 1.301, 1.091, 1.893, 0.489, 1.29];
ubratio1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
bottime1 = [109, 115, 114, 118, 118, 115, 109, 112, 113, 116, 110, 112, 110, 108, 117, 118, 110, 116, 115, 109, 114, 118, 112, 117, 113, 115, 118, 117, 108, 118];
ncut1 = [12314, 13128, 13046, 13120, 12250, 13038, 13066, 12689, 12866, 12912, 12134, 12445, 12153, 12196, 12303, 12939, 12287, 12483, 12434, 13158, 12876, 12816, 12800, 12303, 13043, 12794, 12235, 13023, 12380, 12939];
# phase2
rgap2 = [0.39733392, 4.029e-5, 0.04468459, 0.07029165, 0.00230526, 0.41685617, 0.00259881, 0.23666351, 0.37471668, 0.18886952, 0.37256756, 0.30538008, 3.598e-5, 0.0628511, 2.075e-5, 0.29745276, 0.02307396, 0.07542717, 0.39816345, 0.00493117, 0.34042679, 0.40544532, 1.361e-5, 0.00225484, 0.00511729, 0.29451405, 0.00146045, 0.00320929, 0.00140226, 0.16388265];
ubratio2 = [0.837, 0.688, 0.693, 0.819, 0.251, 0.4, 0.709, 0.701, 0.305, 0.805, 0.459, 1.0, 1.0, 0.546, 0.826, 0.639, 0.862, 0.655, 0.68, 0.422, 0.346, 1.0, 0.601, 0.838, 0.837, 0.614, 0.605, 0.242, 0.668, 0.672];
bottime2 = [222, 237, 232, 344, 238, 286, 246, 228, 345, 239, 342, 268, 252, 345, 247, 298, 282, 220, 237, 330, 257, 266, 342, 265, 269, 337, 260, 324, 344, 296];
ncut2 = [21324, 19527, 19508, 18728, 17891, 21094, 20460, 20837, 21351, 21503, 21350, 21223, 19772, 21011, 19728, 20269, 21350, 18279, 18856, 19342, 20614, 21037, 20586, 20570, 19460, 20844, 20681, 21331, 18645, 18306];
# phase3
rgap3 = [0.00301884, 0.00214521, 0.00082895, 0.00021935, 8.31e-6, 9.299e-5, 0.0014006, 3.58e-6, 0.00037968, 0.00017904, 0.00280017, 0.00071942, 0.00181051, 0.00266093, 0.00154007, 0.00267062, 0.0002106, 0.0017804, 0.00226072, 0.00229041, 0.00251048, 0.00177011, 0.0025602, 0.0007204, 0.00010069, 4.048e-5, 0.00148017, 0.00182003, 0.00021014, 0.0029006];
ubratio3 =  [0.221, 0.416, 0.574, 0.663, 0.437, 0.738, 0.922, 0.406, 0.646, 0.8, 0.359, 0.426, 0.556, 0.439, 0.784, 0.744, 0.671, 0.384, 0.816, 0.406, 0.507, 0.807, 0.473, 0.299, 0.596, 0.787, 0.482, 0.748, 0.725, 0.419];
bottime3 = [391, 424, 368, 356, 350, 418, 369, 377, 434, 435, 402, 345, 347, 319, 324, 395, 300, 360, 379, 333, 428, 379, 390, 407, 438, 420, 401, 375, 441, 379];
ncut3 = [23463, 23179, 22323, 21169, 24231, 23867, 21904, 21002, 23607, 23849, 24326, 24790, 23766, 23182, 22052, 22329, 24775, 23071, 23354, 22695, 22513, 20984, 24504, 22790, 23992, 21637, 21141, 21443, 23970, 23080];

using CairoMakie
Makie.update_theme!(
    fonts = (
        regular = "Times New Roman",
        bold = "Times New Roman Bold",
        italic = "Times New Roman Italic",
        bold_italic = "Times New Roman Bold Italic",
    )
)
fig = Figure(; figure_padding = 4.9);
ssf(s) = rich(s;font="SimSun")
tnf(s) = rich(s;font="Times New Roman")
ax1_1 = Axis(fig[1,1];ylabel=ssf("相对间隙"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax1_2 = Axis(fig[1,2];ylabel=ssf("上界对其初始值的比值"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax1_3 = Axis(fig[2,1];xlabel=tnf("15min"),ylabel=rich("瓶颈时长/",tnf("s");font="SimSun"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax1_4 = Axis(fig[2,2];xlabel=tnf("15min"),ylabel=ssf("割平面数"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax2_1 = Axis(fig[1,3];ylabel=ssf("相对间隙"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax2_2 = Axis(fig[1,4];ylabel=ssf("上界对其初始值的比值"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax2_3 = Axis(fig[2,3];xlabel=tnf("45min"),ylabel=rich("瓶颈时长/",tnf("s");font="SimSun"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax2_4 = Axis(fig[2,4];xlabel=tnf("45min"),ylabel=ssf("割平面数"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax3_1 = Axis(fig[1,5];ylabel=ssf("相对间隙"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax3_2 = Axis(fig[1,6];ylabel=ssf("上界对其初始值的比值"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax3_3 = Axis(fig[2,5];xlabel=tnf("2h"),ylabel=rich("瓶颈时长/",tnf("s");font="SimSun"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
ax3_4 = Axis(fig[2,6];xlabel=tnf("2h"),ylabel=ssf("割平面数"),xticksvisible = false,xticklabelsvisible = false,xgridvisible = false,ygridvisible = false);
boxplot!(ax1_1, ones(length(rgap1)), rgap1);
boxplot!(ax2_1, ones(length(rgap2)), rgap2);
boxplot!(ax3_1, ones(length(rgap3)), rgap3);
boxplot!(ax1_2, ones(length(ubratio1)), ubratio1, color = :darkgoldenrod1);
boxplot!(ax2_2, ones(length(ubratio2)), ubratio2, color = :darkgoldenrod1);
boxplot!(ax3_2, ones(length(ubratio3)), ubratio3, color = :darkgoldenrod1);
boxplot!(ax1_3, ones(length(bottime1)), bottime1, color = :seagreen3);
boxplot!(ax2_3, ones(length(bottime2)), bottime2, color = :seagreen3);
boxplot!(ax3_3, ones(length(bottime3)), bottime3, color = :seagreen3);
boxplot!(ax1_4, ones(length(ncut1)), ncut1, color = :coral2);
boxplot!(ax2_4, ones(length(ncut2)), ncut2, color = :coral2);
boxplot!(ax3_4, ones(length(ncut3)), ncut3, color = :coral2);
fig
